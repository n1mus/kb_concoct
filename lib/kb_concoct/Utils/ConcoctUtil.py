import errno
import json
import os
import subprocess
import sys
import time
import uuid
import zipfile

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.WorkspaceClient import Workspace


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class ConcoctUtil:
    CONCOCT_BASE_PATH = '/kb/deployment/bin/CONCOCT'
    CONCOCT_RESULT_DIRECTORY = 'concoct_output_dir'
    CONCOCT_BIN_RESULT_DIR = 'final_bins'
    MAPPING_THREADS = 16
    BBMAP_MEM = '30g'

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.shock_url = config['shock-url']
        self.ws_url = config['workspace-url']
        self.dfu = DataFileUtil(self.callback_url)
        self.ru = ReadsUtils(self.callback_url)
        self.au = AssemblyUtil(self.callback_url)
        self.mgu = MetagenomeUtils(self.callback_url)

    # Function to return object ID based on object name
    def get_obj_id(object_name):
        """
        Example: get_obj_id("MAG-QC_Archaea.SAGs.Prokka"), returns '32342/1/3'
        """
        ws = Workspace(self.ws_url)
        ws_name = os.environ['KB_WORKSPACE_ID']
        try:
            object_id = ws.get_object_info3({'objects': [{'workspace': ws_name, 'name': object_name}]})['paths'][0][0]
            return object_id
        except:
            return False

    # Function to return object data based on object ID
    def get_object_data(object_id):
        """
        Fetch data from the workspace.
        Example1: get_object_data(get_obj_id("MAG-QC_Archaea.SAGs.RAST"))
        Example2: get_object_data(u'43402/2132/1')
        """
        from biokbase.workspace.client import Workspace
        import os

        ws = Workspace(self.ws_url)
        return ws.get_objects([{"ref":object_id}])[0]


    # Function to return object type based on object ID
    def get_object_type(object_name):
        """
        Fetch data type from the workspace, for example "PairedEndLibrary".
        Example1: get_object_type('lib2.oldstyle.fastq_reads')
        """
        ws = Workspace(self.ws_url)
        return get_object_data(get_object_id(object_name)).values()[0][2]

    def validate_run_concoct_params(self, task_params):
        """
        validate_run_concoct_params:
                validates params passed to run_concoct method
        """
        log('Start validating run_concoct params')

        # check for required parameters
        for p in ['assembly_ref', 'binned_contig_name', 'workspace_name', 'reads_list', 'read_mapping_tool']:
            if p not in task_params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def mkdir_p(self, path):
        """
        mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def run_command(self, command):
        """
        run_command: run command and print result
        """
        os.chdir(self.scratch)
        log('Start executing command:\n{}'.format(command))
        log('Command is running from:\n{}'.format(self.scratch))
        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        output, stderr = pipe.communicate()
        exitCode = pipe.returncode

        if (exitCode == 0):
            log('Executed command:\n{}\n'.format(command) +
                'Exit Code: {}\n'.format(exitCode))
        else:
            error_msg = 'Error running command:\n{}\n'.format(command)
            error_msg += 'Exit Code: {}\nOutput:\n{}\nStderr:\n{}'.format(exitCode, output, stderr)
            raise ValueError(error_msg)
            sys.exit(1)
        return (output, stderr)

    def stage_reads_list_file(self, reads_list):
        """
        stage_reads_list_file: download fastq file associated to reads to scratch area
                          and return result_file_path
        """

        log('Processing reads object list: {}'.format(reads_list))

        result_file_path = []

        # getting from workspace and writing to scratch. The 'reads' dictionary now has file paths to scratch.
        reads = self.ru.download_reads({'read_libraries': reads_list, 'interleaved': None})['files']

        print("\n\n\n\n: reads variable: {}".format(reads))
        print("\n\n\n\n: reads variable: {}".format(list(reads.values())[0]))
        print("\n\n\n\n: reads variable: {}".format(list(reads.values())[13][5]))


        # reads_list is the list of file paths on workspace? (i.e. 12804/1/1).
        # "reads" is the hash of hashes where key is "12804/1/1" or in this case, read_obj and
        # "files" is the secondary key. The tertiary keys are "fwd" and "rev", as well as others.
        for read_obj in reads_list:
            files = reads[read_obj]['files']    # 'files' is dictionary where 'fwd' is key of file path on scratch.
            result_file_path.append(files['fwd'])
            if 'rev' in files and files['rev'] is not None:
                result_file_path.append(files['rev'])

        return result_file_path

    def get_contig_file(self, assembly_ref):
        """
        get_contig_file: get contig file from GenomeAssembly object
        """
        contig_file = self.au.get_assembly_as_fasta({'ref': assembly_ref}).get('path')

        sys.stdout.flush()
        contig_file = self.dfu.unpack_file({'file_path': contig_file})['file_path']

        return contig_file

    def retrieve_and_clean_assembly(self, task_params):
        if os.path.exists(task_params['contig_file_path']):
            assembly = task_params['contig_file_path']
            print("FOUND ASSEMBLY ON LOCAL SCRATCH")
        else:
            # we are on njsw so lets copy it over to scratch
            assembly = self.get_contig_file(task_params['assembly_ref'])

        # remove spaces from fasta headers because that breaks bedtools
        assembly_clean = os.path.abspath(assembly).split('.fa')[0] + "_clean.fa"

        command = '/bin/bash reformat.sh in={} out={} addunderscore'.format(assembly, assembly_clean)

        log('running reformat command: {}'.format(command))
        out, err = self.run_command(command)

        return assembly_clean

    def generate_stats_for_genome_bins(self, task_params, genome_bin_fna_file, bbstats_output_file):
        """
        generate_command: bbtools stats.sh command
        """
        log("\n\nRunning generate_stats_for_genome_bins on {}".format(genome_bin_fna_file))
        result_directory = os.path.join(self.scratch, self.CONCOCT_RESULT_DIRECTORY)
        genome_bin_fna_file = os.path.join(self.scratch, self.CONCOCT_RESULT_DIRECTORY, genome_bin_fna_file)
        command = '/bin/bash stats.sh in={} format=3 > {}'.format(genome_bin_fna_file, bbstats_output_file)
        self.run_command(command)
        bbstats_output = open(bbstats_output_file, 'r').readlines()[1]
        n_scaffolds = bbstats_output.split('\t')[0]
        n_contigs = bbstats_output.split('\t')[1]
        scaf_bp = bbstats_output.split('\t')[2]
        contig_bp = bbstats_output.split('\t')[3]
        gap_pct = bbstats_output.split('\t')[4]
        scaf_N50 = bbstats_output.split('\t')[5]
        scaf_L50 = bbstats_output.split('\t')[6]
        ctg_N50 = bbstats_output.split('\t')[7]
        ctg_L50 = bbstats_output.split('\t')[8]
        scaf_N90 = bbstats_output.split('\t')[9]
        scaf_L90 = bbstats_output.split('\t')[10]
        ctg_N90 = bbstats_output.split('\t')[11]
        ctg_L90 = bbstats_output.split('\t')[12]
        scaf_max = bbstats_output.split('\t')[13]
        ctg_max = bbstats_output.split('\t')[14]
        scaf_n_gt50K = bbstats_output.split('\t')[15]
        scaf_pct_gt50K = bbstats_output.split('\t')[16]
        gc_avg = float(bbstats_output.split('\t')[17]) * 100 # need to figure out if correct
        gc_std = float(bbstats_output.split('\t')[18]) * 100 # need to figure out if correct

        log('Generated generate_stats_for_genome_bins command: {}'.format(command))

        return {'n_scaffolds': n_scaffolds, 'n_contigs': n_contigs, 'scaf_bp': scaf_bp, 'contig_bp': contig_bp, 'gap_pct': gap_pct, 'scaf_N50': scaf_N50, 'scaf_L50': scaf_L50, 'ctg_N50': ctg_N50, 'scaf_N90': scaf_N90, 'ctg_N90': ctg_N90, 'ctg_L90': ctg_L90, 'scaf_max': scaf_max, 'ctg_max': ctg_max, 'scaf_n_gt50K': scaf_n_gt50K, 'gc_avg': gc_avg, 'gc_std':gc_std}



    def run_read_mapping_unpaired_mode(self, task_params, assembly_clean, fastq, sam):
        if task_params['read_mapping_tool'] == 'bbmap':
            command = 'bbmap.sh -Xmx{} fast threads={} ref={} in={} out={} interleaved=false mappedonly nodisk overwrite'.format(self.BBMAP_MEM, self.MAPPING_THREADS, assembly_clean, fastq, sam) # BBMap is deterministic without the deterministic flag if using single-ended reads
        elif task_params['read_mapping_tool'] == 'bwa':
            command = 'bwa index {} && bwa mem -t {} {} {} > {}'.format(assembly_clean, self.MAPPING_THREADS, assembly_clean, fastq, sam)
        elif task_params['read_mapping_tool'] == 'bowtie2':
            bt2index = os.path.basename(assembly_clean) + '.bt2'
            command = 'bowtie2-build -f {} --threads {} {} && bowtie2 -x {} -U {} --threads {} -S {}'.format(assembly_clean, self.MAPPING_THREADS, bt2index, bt2index, fastq, self.MAPPING_THREADS, sam)
        elif task_params['read_mapping_tool'] == 'minimap2':
            command = 'minimap2 -ax sr -t {} {} {} > {}'.format(self.MAPPING_THREADS, assembly_clean, fastq, sam)
        log('running alignment command: {}'.format(command))
        out, err = self.run_command(command)


    def deinterlace_raw_reads(self, fastq):
        fastq_forward = fastq.split('.fastq')[0] + "_forward.fastq"
        fastq_reverse = fastq.split('.fastq')[0] + "_reverse.fastq"
        command = 'reformat.sh in={} out1={} out2={}'.format(fastq, fastq_forward, fastq_reverse)
        self.run_command(command)
        return (fastq_forward, fastq_reverse)


    def run_read_mapping_interleaved_pairs_mode(self, task_params, assembly_clean, fastq, sam):
        if task_params['read_mapping_tool'] == 'bbmap':
            command = 'bbmap.sh -Xmx{} fast threads={} ref={} in={} out={} interleaved=true mappedonly nodisk overwrite'.format(self.BBMAP_MEM, self.MAPPING_THREADS, assembly_clean, fastq, sam)
        elif task_params['read_mapping_tool'] == 'bwa':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            command = 'bwa index {} && bwa mem -t {} {} {} {} > {}'.format(assembly_clean, self.MAPPING_THREADS, assembly_clean, fastq_forward, fastq_reverse, sam)
        elif task_params['read_mapping_tool'] == 'bowtie2':
            bt2index = os.path.basename(assembly_clean) + '.bt2'
            command = 'bowtie2-build -f {} --threads {} {} && bowtie2 -x {} --interleaved {} --threads {} -S {}'.format(assembly_clean, self.MAPPING_THREADS, bt2index, bt2index, fastq, self.MAPPING_THREADS, sam)
        elif task_params['read_mapping_tool'] == 'minimap2':
            #need to split interleaved read file to run in paired mode
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            command = 'minimap2 -ax sr -t {} {} {} {} > {}'.format(self.MAPPING_THREADS, assembly_clean, fastq_foward, fastq_reverse, sam)
        log('running alignment command: {}'.format(command))
        out, err = self.run_command(command)


    def convert_sam_to_sorted_and_indexed_bam(self, sam, sorted_bam):
        # create bam files from sam files
        command = 'samtools view -F 0x04 -uS {} | samtools sort - -o {}'.format(sam, sorted_bam)

        log('running samtools command to generate sorted bam: {}'.format(command))
        self.run_command(command)

        # verify we got bams
        if not os.path.exists(sorted_bam):
            log('Failed to find bam file\n{}'.format(sorted_bam))
            sys.exit(1)
        elif(os.stat(sorted_bam).st_size == 0):
            log('Bam file is empty\n{}'.format(sorted_bam))
            sys.exit(1)

        # index the bam file
        command = 'samtools index {}'.format(sorted_bam)

        log('running samtools command to index sorted bam: {}'.format(command))
        self.run_command(command)

    def generate_alignment_bams(self, task_params):
        """
            This function runs the bbmap.sh alignments and creates the
            bam files from sam files using samtools.
        """
        # first, clean the assembly file so that there are no spaces in the fasta headers
        assembly_clean = self.retrieve_and_clean_assembly(task_params)

        # even though read_name is just one file, make read_name a
        # list because that is what stage_reads_list_file expects
        reads_list = task_params['reads_list']

        # we need to copy the fastq file from workspace to local scratch.
        # 'fastq' should be path to one file (if we had read.fwd.fq and read.rev.fq then they
        # should have been interleaved already). However, I should interleave them in case
        # this function gets called other than the narrative.
        # this should return a list of 1 fastq file path on scratch
        read_scratch_path = self.stage_reads_list_file(reads_list)

        print("read_scratch_path {}".format(read_scratch_path))
        print("read_scratch_path {}".format(read_scratch_path[0]))

        # list of reads files, can be 1 or more. assuming reads are either type unpaired or interleaved
        # will not handle unpaired forward and reverse reads input as seperate (non-interleaved) files
        task_params['read_scratch_file'] = read_scratch_path

        for i in range(len(read_scratch_path)):
            fastq = read_scratch_path[i]
            #fastq = read_scratch_path[0] #needs to be a loop
            #reads_type = self.get_obj_id(str(os.path.basename(fastq)[1]))
            print("OS.PATH.BASENAME0 {}".format([os.path.basename(fastq)]))
            print("OS.PATH.BASENAME1 {}".format(type([os.path.basename(fastq)])))

            #print("OS.PATH.BASENAME" + str(self.get_obj_id(os.path.basename(fastq))))
            # reads_type = self.get_object_type(str(os.path.basename(fastq)))
            # reads_type = self.get_object_type(self.get_obj_id([os.path.basename(fastq)])).values()[0][2]
            # print("read_type is " + reads_type)

            #result_directory = task_params['result_directory']
            sam = os.path.basename(fastq).split('.fastq')[0] + ".sam"
            sam = os.path.join(self.CONCOCT_RESULT_DIRECTORY, sam)

            self.run_read_mapping_unpaired_mode(task_params, assembly_clean, fastq, sam)

            sorted_bam = os.path.abspath(sam).split('.sam')[0] + "_sorted.bam"

            # need to get this to work on list
            task_params['sorted_bam'] = sorted_bam

            self.convert_sam_to_sorted_and_indexed_bam(sam, sorted_bam)

        return assembly_clean

    def generate_make_coverage_table_command(self, task_params):
        # create the depth file for this bam
        #
        min_contig_length = task_params['min_contig_length']

        depth_file_path = os.path.join(self.scratch, str('concoct_depth.txt'))
        command = '/kb/module/lib/kb_concoct/bin/jgi_summarize_bam_contig_depths '
        command += '--outputDepth {} --minContigLength {} --minContigDepth 1 '.format(depth_file_path, min_contig_length)
        command += task_params['sorted_bam']

        log('running summarize_bam_contig_depths command: {}'.format(command))
        self.run_command(command)

        return depth_file_path

    # untested, need to check in UI
    def fix_generate_concoct_command_ui_bug(self, task_params):
        # needed to get checkbox for UI to work with string objects, for some reason strings are converted to numerics when running inside KBase UI.
        parameter_no_total_coverage = task_params['no_total_coverage']
        parameter_no_cov_normalization = task_params['no_cov_normalization']

        if task_params['no_total_coverage'] is 1:
            parameter_no_total_coverage = '--no_total_coverage'
        elif task_params['no_total_coverage'] is 0:
            parameter_no_total_coverage = ' '

        if task_params['no_cov_normalization'] is 1:
            parameter_no_cov_normalization = '--no_cov_normalization'
        elif task_params['no_cov_normalization'] is 0:
            parameter_no_cov_normalization = ' '

        return (parameter_no_total_coverage, parameter_no_cov_normalization)



    def generate_concoct_cut_up_fasta_command(self, task_params):
        """
        generate_command: concoct cut_up_fasta
        """
        log("\n\nRunning generate_concoct_cut_up_fasta_command")

        command = 'python {}/scripts/cut_up_fasta.py {} -c {} -o {} --merge_last -b temp.bed > {}/split_contigs.fa'.format(self.CONCOCT_BASE_PATH, task_params['contig_file_path'], task_params['contig_split_size'],task_params['contig_split_overlap'],self.CONCOCT_RESULT_DIRECTORY)
        log('Generated concoct_cut_up_fasta command: {}'.format(command))

        return command



    def generate_concoct_coverage_table_from_bam(self, task_params):
        """
        generate_command: concoct generate coverage table
        """
        log("\n\nRunning generate_concoct_coverage_table_from_bam")
        command = 'python {}/scripts/concoct_coverage_table.py temp.bed {}/*_sorted.bam > {}/coverage_table.tsv'.format(self.CONCOCT_BASE_PATH, self.CONCOCT_RESULT_DIRECTORY, self.CONCOCT_RESULT_DIRECTORY)
        log('Generated concoct generate coverage table from bam command: {}'.format(command))

        return command


    def generate_concoct_command(self, task_params):
        """
        generate_command: concoct
        """
        parameter_no_total_coverage, parameter_no_cov_normalization = self.fix_generate_concoct_command_ui_bug(task_params)

        log("\n\nRunning generate_concoct_command")
        command = 'python {}/bin/concoct --composition_file {}/split_contigs.fa -l {} -b {} --coverage_file {}/coverage_table.tsv -t {} -k {} -c {} -i {} --total_percentage_pca {} {} {}'.format(self.CONCOCT_BASE_PATH, self.CONCOCT_RESULT_DIRECTORY, task_params['min_contig_length'], self.CONCOCT_RESULT_DIRECTORY, self.CONCOCT_RESULT_DIRECTORY,
                                         self.MAPPING_THREADS,task_params['kmer_size'],task_params['max_clusters_for_vgmm'],task_params['max_iterations_for_vgmm'],task_params['total_percentage_pca'],parameter_no_cov_normalization,parameter_no_total_coverage)
        log('Generated concoct command: {}'.format(command))

        return command

    def generate_concoct_post_clustering_merging_command(self, task_params):
        """
        generate_command: concoct post cluster merging
        """

        log("\n\nRunning generate_concoct_post_clustering_merging_command")

        command = 'python {}/scripts/merge_cutup_clustering.py {}/clustering_gt{}.csv > {}/clustering_merged.csv'.format(self.CONCOCT_BASE_PATH, self.CONCOCT_RESULT_DIRECTORY, task_params['min_contig_length'], self.CONCOCT_RESULT_DIRECTORY)
        log('Generated generate_concoct_post_clustering_merging command: {}'.format(command))

        return command

    def generate_concoct_extract_fasta_bins_command(self, task_params):
        """
        generate_command: concoct extract_fasta_bins
        """
        log("\n\nRunning generate_concoct_extract_fasta_bins_command")

        command = 'mkdir {}/{}'.format(self.CONCOCT_RESULT_DIRECTORY, self.CONCOCT_BIN_RESULT_DIR)
        command += ' && '
        command += 'python {}/scripts/extract_fasta_bins.py {} {}/clustering_merged.csv --output_path {}/{}'.format(self.CONCOCT_BASE_PATH, task_params['contig_file_path'],self.CONCOCT_RESULT_DIRECTORY, self.CONCOCT_RESULT_DIRECTORY, self.CONCOCT_BIN_RESULT_DIR)
        log('Generated generate_concoct_extract_fasta_bins_command command: {}'.format(command))

        return command


    def rename_and_standardize_bin_names(self, task_params):
        """
        generate_command: generate renamed bins
        """
        log("\n\nRunning rename_and_standardize_bin_names")
        path_to_concoct_result_bins = os.path.abspath(self.CONCOCT_RESULT_DIRECTORY) + '/' + self.CONCOCT_BIN_RESULT_DIR + '/'
        for dirname, subdirs, files in os.walk(path_to_concoct_result_bins):
            for file in files:
                if file.endswith('.fa'):
                    os.rename(os.path.abspath(path_to_concoct_result_bins) + '/' + file, os.path.abspath(path_to_concoct_result_bins) + '/bin.' + file.split('.fa')[0].zfill(3) + '.fasta')


    def make_binned_contig_summary_file_for_binning_apps(self, task_params):
        """
        generate_command: generate binned contig summary command
        """
        log("\n\nRunning make_binned_contig_summary_file_for_binning_apps")
        path_to_concoct_result = os.path.abspath(self.CONCOCT_RESULT_DIRECTORY)
        path_to_concoct_result_bins = '{}/{}/'.format(path_to_concoct_result, self.CONCOCT_BIN_RESULT_DIR)
        path_to_summary_file = path_to_concoct_result_bins + 'binned_contig.summary'
        with open(path_to_summary_file, 'w+') as f:
            f.write("Bin name\tCompleteness\tGenome size\tGC content\n")
            for dirname, subdirs, files in os.walk(path_to_concoct_result_bins):
                for file in files:
                    if file.endswith('.fasta'):
                         genome_bin_fna_file = os.path.join(self.CONCOCT_BIN_RESULT_DIR, file)
                         bbstats_output_file = os.path.join(self.scratch, self.CONCOCT_RESULT_DIRECTORY, genome_bin_fna_file).split('.fasta')[0] + ".bbstatsout"
                         bbstats_output = self.generate_stats_for_genome_bins(task_params, genome_bin_fna_file, bbstats_output_file)
                         f.write('{}\t0\t{}\t{}\n'.format(genome_bin_fna_file.split("/")[-1], bbstats_output['contig_bp'],bbstats_output['gc_avg']))
        f.close()
        log('Finished make_binned_contig_summary_file_for_binning_apps function')


    def generate_output_file_list(self, result_directory):
        """
        generate_output_file_list: zip result files and generate file_links for report
        """
        log('Start packing result files')
        output_files = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self.mkdir_p(output_directory)
        result_file = os.path.join(output_directory, 'concoct_result.zip')

        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:

            for dirname, subdirs, files in os.walk(result_directory):
                for file in files:
                    if (file.endswith('.sam') or file.endswith('.bam') or file.endswith('.bai') or file.endswith('.summary')):
                            continue
                    if (dirname.endswith(self.CONCOCT_BIN_RESULT_DIR)):
                            continue
                    zip_file.write(os.path.join(dirname, file), file)
                if (dirname.endswith(self.CONCOCT_BIN_RESULT_DIR)):
                    baseDir = os.path.basename(dirname)
                    for file in files:
                        full = os.path.join(dirname, file)
                        zip_file.write(full, os.path.join(baseDir, file))

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'Files generated by CONCOCT App'})

        return output_files

    def generate_html_report(self, result_directory, assembly_ref, binned_contig_obj_ref):
        """
        generate_html_report: generate html summary report
        """

        log('Start generating html report')
        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self.mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'report.html')

        # get summary data from existing assembly object and bins_objects
        Summary_Table_Content = ''
        Overview_Content = ''
        (binned_contig_count, input_contig_count, total_bins_count) = self.generate_overview_info(assembly_ref, binned_contig_obj_ref, result_directory)

        Overview_Content += '<p>Binned contigs: {}</p>'.format(binned_contig_count)
        Overview_Content += '<p>Input contigs: {}</p>'.format(input_contig_count)
        Overview_Content += '<p>Number of bins: {}</p>'.format(total_bins_count)

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          Overview_Content)
                report_template = report_template.replace('Summary_Table_Content',
                                                          Summary_Table_Content)
                result_file.write(report_template)

        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for kb_concoct App'})
        return html_report

    def generate_overview_info(self, assembly_ref, binned_contig_obj_ref, result_directory):
        """
        _generate_overview_info: generate overview information from assembly and binnedcontig
        """

        # get assembly and binned_contig objects that already have some data populated in them
        assembly = self.dfu.get_objects({'object_refs': [assembly_ref]})['data'][0]
        binned_contig = self.dfu.get_objects({'object_refs': [binned_contig_obj_ref]})['data'][0]

        input_contig_count = assembly.get('data').get('num_contigs')
        binned_contig_count = 0
        total_bins_count = 0
        total_bins = binned_contig.get('data').get('bins')
        total_bins_count = len(total_bins)
        for bin in total_bins:
            binned_contig_count += len(bin.get('contigs'))

        return (binned_contig_count, input_contig_count, total_bins_count)

    def generate_report(self, binned_contig_obj_ref, task_params):
        """
        generate_report: generate summary report
        """
        log('Generating report')

        result_directory = os.path.join(self.scratch, "concoct_output_dir")

        task_params['result_directory'] = result_directory

        output_files = self.generate_output_file_list(task_params['result_directory'])

        output_html_files = self.generate_html_report(task_params['result_directory'], task_params['assembly_ref'], binned_contig_obj_ref)

        report_params = {
              'message': '',
              'workspace_name': task_params['workspace_name'],
              'file_links': output_files,
              'html_links': output_html_files,
              'direct_html_link_index': 0,
              'html_window_height': 266,
              'report_object_name': 'kb_concoct_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def create_dict_from_depth_file(self, depth_file_path):
        # keep contig order (required by metabat2)
        depth_file_dict = {}
        with open(depth_file_path, 'r') as f:
            header = f.readline().rstrip().split("\t")
            # print('HEADER1 {}'.format(header))
            # map(str.strip, header)
            for line in f:
                # deal with cases were fastq name has spaces.Assume first
                # non white space word is unique and use this as ID.
                # line = line.rstrip()
                vals = line.rstrip().split("\t")
                if ' ' in vals[0]:
                    ID = vals[0].split()[0]
                else:
                    ID = vals[0]
                depth_file_dict[ID] = vals[1:]
            depth_file_dict['header'] = header
        return depth_file_dict

    def run_concoct(self, task_params):
        """
        run_concoct: concoct app

        required params:
            assembly_ref: Metagenome assembly object reference
            binned_contig_name: BinnedContig object name and output file header
            workspace_name: the name of the workspace it gets saved to.
            reads_list: list of reads object (PairedEndLibrary/SingleEndLibrary)
            upon which CONCOCT will be run

        optional params:
            min_contig_length: minimum contig length; default 1000

            ref: https://github.com/BinPro/CONCOCT/blob/develop/README.md
        """
        log('--->\nrunning ConcoctUtil.run_concoct\n' +
            'task_params:\n{}'.format(json.dumps(task_params, indent=1)))

        self.validate_run_concoct_params(task_params)

        contig_file = self.get_contig_file(task_params['assembly_ref'])
        task_params['contig_file_path'] = contig_file

        reads_list_file = self.stage_reads_list_file(task_params['reads_list'])
        task_params['reads_list_file'] = reads_list_file

        result_directory = os.path.join(self.scratch, self.CONCOCT_RESULT_DIRECTORY)

        self.mkdir_p(result_directory)

        cwd = os.getcwd()
        log('changing working dir to {}'.format(result_directory))
        os.chdir(result_directory)

        #run alignments, and update input contigs to use the clean file

        # this is the clean assembly
        task_params['contig_file_path'] = self.generate_alignment_bams(task_params)

        # not used right now
        depth_file_path = self.generate_make_coverage_table_command(task_params)
        depth_dict = self.create_dict_from_depth_file(depth_file_path)

        # run concoct prep, cut up fasta input
        command = self.generate_concoct_cut_up_fasta_command(task_params)
        self.run_command(command)

        # run concoct make coverage table from bam
        command = self.generate_concoct_coverage_table_from_bam(task_params)
        self.run_command(command)

        # run concoct prep and concoct
        command = self.generate_concoct_command(task_params)
        self.run_command(command)

        # run concoct post cluster merging command
        command = self.generate_concoct_post_clustering_merging_command(task_params)
        self.run_command(command)

        # run extract bins command
        command = self.generate_concoct_extract_fasta_bins_command(task_params)
        self.run_command(command)

        # run fasta renaming
        self.rename_and_standardize_bin_names(task_params)

        self.make_binned_contig_summary_file_for_binning_apps(task_params)


        # file handling and management
        os.chdir(cwd)
        log('changing working dir to {}'.format(cwd))

        log('Saved result files to: {}'.format(result_directory))
        log('Generated files:\n{}'.format('\n'.join(os.listdir(result_directory))))

        # make new BinnedContig object and upload to KBase
        generate_binned_contig_param = {
            'file_directory': os.path.join(result_directory, self.CONCOCT_BIN_RESULT_DIR),
            'assembly_ref': task_params['assembly_ref'],
            'binned_contig_name': task_params['binned_contig_name'],
            'workspace_name': task_params['workspace_name']
        }
        binned_contig_obj_ref = self.mgu.file_to_binned_contigs(generate_binned_contig_param).get('binned_contig_obj_ref')

        # generate report
        reportVal = self.generate_report(binned_contig_obj_ref, task_params)
        returnVal = {
            'result_directory': result_directory,
            'binned_contig_obj_ref': binned_contig_obj_ref
        }
        returnVal.update(reportVal)

        return returnVal
