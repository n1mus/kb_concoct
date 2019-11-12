import errno
import json
import os
import subprocess
import sys
import time
import uuid
import zipfile

from Bio import SeqIO

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from installed_clients.ReadsUtilsClient import ReadsUtils


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class ConcoctUtil:
    CONCOCT_TOOLKIT_PATH = '/kb/deployment/bin/CONCOCT/bin'
    CONCOCT_BIN_DIR = 'final_bins'
    BBMAP_THREADS=16
    BBMAP_MEM='30g'

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.shock_url = config['shock-url']
        self.dfu = DataFileUtil(self.callback_url)
        self.ru = ReadsUtils(self.callback_url)
        self.au = AssemblyUtil(self.callback_url)
        self.mgu = MetagenomeUtils(self.callback_url)


    def validate_run_concoct_params(self, params):
        """
        validate_run_concoct_params:
                validates params passed to run_concoct method

        """
        log('Start validating run_concoct params')

        # check for required parameters
        for p in ['assembly_ref', 'binned_contig_name', 'workspace_name', 'reads_list', 'read_mapping_tool']:
            if p not in params:
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
        output,stderr = pipe.communicate()
        exitCode = pipe.returncode

        if (exitCode == 0):
            log('Executed command:\n{}\n'.format(command) +
                'Exit Code: {}\n'.format(exitCode))
        else:
            error_msg = 'Error running command:\n{}\n'.format(command)
            error_msg += 'Exit Code: {}\nOutput:\n{}\nStderr:\n{}'.format(exitCode, output, stderr)
            raise ValueError(error_msg)
            sys.exit(1)
        return (output,stderr)

    def stage_reads_list_file(self, reads_list):
        """
        stage_reads_list_file: download fastq file associated to reads to scratch area
                          and return result_file_path
        """

        log('Processing reads object list: {}'.format(reads_list))

        result_file_path = []

        # getting from workspace and writing to scratch. The 'reads' dictionary now has file paths to scratch.
        reads = self.ru.download_reads({'read_libraries': reads_list, 'interleaved': None})['files']

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
        get_contig_file: get contif file from GenomeAssembly object
        """

        contig_file = self.au.get_assembly_as_fasta({'ref': assembly_ref}).get('path')

        sys.stdout.flush()
        contig_file = self.dfu.unpack_file({'file_path': contig_file})['file_path']

        return contig_file

    def retrieve_and_clean_assembly(self, task_params):
        if os.path.exists(task_params['contig_file_path']):
            assembly =  task_params['contig_file_path']
            print("FOUND ASSEMBLY ON LOCAL SCRATCH")
        else:
            # we are on njsw so lets copy it over to scratch
            assembly = self.get_contig_file(task_params['assembly_ref'])

        # remove spaces from fasta headers because that breaks bedtools
        assembly_clean = os.path.abspath(assembly).split('.fa')[0] + "_clean.fa"

        command = '/bin/bash reformat.sh in={} out={} addunderscore'.format(assembly,assembly_clean)

        log('running reformat command: {}'.format(command))
        out,err = self.run_command(command)

        return assembly_clean


    def generate_alignment_bams(self, task_params):
        """
            This function runs the bbmap.sh alignments and creates the
            bam files from sam files using samtools.
        """
        # even though read_name is just one file, make read_name a
        # list because that is what stage_reads_list_file expects

        reads_list = task_params['reads_list']

        # we need to copy the fastq file from workspace to local scratch.
        # 'fastq' should be path to one file (if we had read.fwd.fq and read.rev.fq then they
        # should have been interleaved already). However, I should interleave them in case
        # this function gets called other than the narrative.
        # this should return a list of 1 fastq file path on scratch
        read_scratch_path = self.stage_reads_list_file(reads_list)

        task_params['read_scratch_file'] = read_scratch_path

        assembly_clean = self.retrieve_and_clean_assembly(task_params)

        fastq = read_scratch_path[0]
        result_directory = task_params['result_directory']
        sam = os.path.basename(fastq) + '.sam'
        sam = os.path.join(result_directory, sam)
        sorted_bam = sam + '.sorted.bam'

        if task_params['read_mapping_tool'] is 'bbmap':
            command = '/bin/bash bbmap.sh -Xmx{} fast threads={} ref={} in={} out={} mappedonly nodisk overwrite'.format(self.BBMAP_MEM,self.BBMAP_THREADS,assembly_clean,fastq,sam)
        elif task_params['read_mapping_tool'] is 'bwa':
            command = 'bwa index {} && bwa mem -t {} {} {} > {}'.format(assembly_clean,self.BBMAP_THREADS,assembly_clean,fastq,sam)
        elif task_params['read_mapping_tool'] is 'bowtie2':
            bt2index = os.path.basename(assembly_clean) + '.bt2'
            command = 'bowtie2-build -f {} --threads {} {} && bowtie2 -x {} -U {} --threads {} -S {}'.format(assembly_clean,self.BBMAP_THREADS,bt2index,bt2index,fastq,self.BBMAP_THREADS,sam)
        elif task_params['read_mapping_tool'] is 'minimap2':
            command = 'minimap2 -ax sr -t {} {} {} > {}'.format(self.BBMAP_THREADS,assembly_clean,fastq,sam)

        log('running alignment command: {}'.format(command))
        out,err = self.run_command(command)

        task_params['sorted_bam'] = sorted_bam

        # create bam files from sam files
        command = 'samtools view -F 0x04 -uS {} | samtools sort - -o {}'.format(sam,sorted_bam)

        log('running samtools command: {}'.format(command))
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

        log('running samtools command: {}'.format(command))
        self.run_command(command)
        return assembly_clean

    def generate_make_coverage_table_command(self, task_params):
        # create the depth file for this bam
        #

        min_contig_length = task_params['min_contig_length']

        depth_file_path  = os.path.join(self.scratch, str('concoct_depth.txt'))
        command = '/kb/module/lib/kb_concoct/bin/jgi_summarize_bam_contig_depths '
        command += '--outputDepth {} --minContigLength {} --minContigDepth 1 '.format(depth_file_path, min_contig_length)
        command += task_params['sorted_bam']

        log('running summarize_bam_contig_depths command: {}'.format(command))
        self.run_command(command)

        return depth_file_path

    def generate_concoct_command(self, params):
        """
        generate_command: generate concoct params
        """

        # needed to get checkbox for UI to work with string objects, for some reason strings are converted to numerics when running inside KBase UI.
        parameter_no_total_coverage = params.get('no_total_coverage')
        parameter_no_cov_normalization = params.get('no_cov_normalization')

        if params.get('no_total_coverage') is 1:
            parameter_no_total_coverage = '--no_total_coverage'
        elif params.get('no_total_coverage') is 0:
            parameter_no_total_coverage = ' '

        if params.get('no_cov_normalization') is 1:
            parameter_no_cov_normalization = '--no_cov_normalization'
        elif params.get('no_cov_normalization') is 0:
            parameter_no_cov_normalization = ' '


        print("\n\nRunning generate_concoct_command")

        print("no_cov_normalization is: " + str(params.get('no_cov_normalization')))

        command = 'python /kb/deployment/bin/CONCOCT/scripts/cut_up_fasta.py {} -c {} -o {} --merge_last -b temp.bed > concoct_output_dir/split_contigs.fa'.format(params.get('contig_file_path'),params.get('contig_split_size'),params.get('contig_split_overlap'))
        command += ' && '
        command += 'python /kb/deployment/bin/CONCOCT/scripts/concoct_coverage_table.py temp.bed concoct_output_dir/*.sorted.bam > concoct_output_dir/coverage_table.tsv'
        command += ' && '
        command += 'python /kb/deployment/bin/CONCOCT/bin/concoct --composition_file concoct_output_dir/split_contigs.fa -l {} -b concoct_output_dir --coverage_file concoct_output_dir/coverage_table.tsv -t {} -k {} -c {} -i {} --total_percentage_pca {} {} {}'.format(params.get('min_contig_length'),
                                         self.BBMAP_THREADS,params.get('kmer_size'),params.get('max_clusters_for_vgmm'),params.get('max_iterations_for_vgmm'),params.get('total_percentage_pca'),parameter_no_cov_normalization,parameter_no_total_coverage)
        command += ' && '
        command += 'python /kb/deployment/bin/CONCOCT/scripts/merge_cutup_clustering.py concoct_output_dir/clustering_gt{}.csv > concoct_output_dir/clustering_merged.csv'.format(params.get('min_contig_length'))
        command += ' && '
        command += 'mkdir concoct_output_dir/{}'.format(self.CONCOCT_BIN_DIR)
        command += ' && '
        command += 'python /kb/deployment/bin/CONCOCT/scripts/extract_fasta_bins.py {} concoct_output_dir/clustering_merged.csv --output_path concoct_output_dir/{}'.format(params.get('contig_file_path'),self.CONCOCT_BIN_DIR)


        log('Generated run_concoct command: {}'.format(command))

        return command


    def generate_summary_utils_command(self, params):
        """
        generate_command: generate summary utils command
        """

        print("\n\nRunning generate_summary_utils_command")

        command = "/bin/bash /kb/module/lib/kb_concoct/bin/summary_utils.sh"

        log('Generated summary_utils command: {}'.format(command))

        return command


    def generate_output_file_list(self, result_directory):
        """
        generate_output_file_list: zip result files and generate file_links for report
        """
        log('Start packing result files')
        output_files = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self.mkdir_p(output_directory)
        result_file = os.path.join(output_directory, 'concoct_result.zip')
        report_file = None

        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:

            for dirname, subdirs, files in os.walk(result_directory):
                for file in files:
                    if (file.endswith('.sam') or file.endswith('.bam') or file.endswith('.bai') or file.endswith('.summary')):
                            continue
                    if (dirname.endswith(self.CONCOCT_BIN_DIR)):
                            continue
                    zip_file.write(os.path.join(dirname, file),file)
                if (dirname.endswith(self.CONCOCT_BIN_DIR)):
                    baseDir=os.path.basename(dirname)
                    for file in files:
                        full=os.path.join(dirname, file)
                        zip_file.write(full,os.path.join(baseDir,file))


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
        (binned_contig_count, input_contig_count,
         total_bins_count) = self.generate_overview_info(assembly_ref,
                                                          binned_contig_obj_ref,
                                                          result_directory)

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
        bins_directory = os.path.join(self.scratch, result_directory, self.CONCOCT_BIN_DIR)
        binned_contig_count = 0
        total_bins_count = 0
        total_bins = binned_contig.get('data').get('bins')
        total_bins_count = len(total_bins)
        for bin in total_bins:
            binned_contig_count += len(bin.get('contigs'))

        return (binned_contig_count, input_contig_count, total_bins_count)


    def generate_report(self, binned_contig_obj_ref, params):
        """
        generate_report: generate summary report

        """
        log('Generating report')

        result_directory = os.path.join(self.scratch, "concoct_output_dir")

        params['result_directory'] = result_directory

        output_files = self.generate_output_file_list(params['result_directory'])

        output_html_files = self.generate_html_report(params['result_directory'],
                                                       params['assembly_ref'],
                                                       binned_contig_obj_ref)

        report_params = {
              'message': '',
              'workspace_name': params.get('workspace_name'),
              'file_links': output_files,
              'html_links': output_html_files,
              'direct_html_link_index': 0,
              'html_window_height': 266,
              'report_object_name': 'kb_concoct_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output


    def createDictFromDepthfile(self, depth_file_path):
        # keep contig order (required by metabat2)
        depth_file_dict = {}
        with open(depth_file_path,'r') as f:
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


    def run_concoct(self, params):
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
            'params:\n{}'.format(json.dumps(params, indent=1)))

        self.validate_run_concoct_params(params)

        #params['split_contig_length'] = params.get('split_contig_length')

        contig_file = self.get_contig_file(params.get('assembly_ref'))
        params['contig_file_path'] = contig_file

        reads_list_file = self.stage_reads_list_file(params.get('reads_list'))
        params['reads_list_file'] = reads_list_file

        result_directory = os.path.join(self.scratch, "concoct_output_dir")

        params['result_directory'] = result_directory

        self.mkdir_p(result_directory)

        cwd = os.getcwd()
        log('changing working dir to {}'.format(result_directory))
        os.chdir(result_directory)
        #run alignments, and update input contigs to use the clean file

        params['contig_file_path'] = self.generate_alignment_bams(params)

        depth_file_path = self.generate_make_coverage_table_command(params)

        depth_dict = self.createDictFromDepthfile(depth_file_path)


        #run concoct
        command = self.generate_concoct_command(params)
        self.run_command(command)

        #run summary utilities
        command = self.generate_summary_utils_command(params)
        self.run_command(command)

        #file handling and management
        os.chdir(cwd)
        log('changing working dir to {}'.format(cwd))

        log('Saved result files to: {}'.format(result_directory))
        log('Generated files:\n{}'.format('\n'.join(os.listdir(result_directory))))

        generate_binned_contig_param = {
            'file_directory': os.path.join(result_directory,self.CONCOCT_BIN_DIR),
            'assembly_ref': params['assembly_ref'],
            'binned_contig_name': params['binned_contig_name'],
            'workspace_name': params['workspace_name']
        }

        binned_contig_obj_ref = self.mgu.file_to_binned_contigs(generate_binned_contig_param).get('binned_contig_obj_ref')

        reportVal = self.generate_report(binned_contig_obj_ref, params)

        returnVal = {
            'result_directory': result_directory,
            'binned_contig_obj_ref': binned_contig_obj_ref
        }

        returnVal.update(reportVal)

        return returnVal
