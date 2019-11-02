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
        for p in ['assembly_ref', 'binned_contig_name', 'workspace_name', 'reads_list']:
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


    def generate_alignment_bams(self, task_params):
        """
            This function runs the bbmap.sh alignments and creates the
            bam files from sam files using samtools.
        """
        # even though read_name is just one file, make read_name a
        # list because that is what stage_reads_list_file expects
        print("task_params******")
        print(str(task_params))
        reads_list = task_params['reads_list']

        # we need to copy the fastq file from workspace to local scratch.
        # 'fastq' should be path to one file (if we had read.fwd.fq and read.rev.fq then they
        # should have been interleaved already). However, I should interleave them in case
        # this function gets called other than the narrative.
        # this should return a list of 1 fastq file path on scratch
        read_scratch_path = self.stage_reads_list_file(reads_list)

        task_params['read_scratch_file'] = read_scratch_path

        #
        # grab the assembly. It should already exist on 'local' node but
        # needs to be copied over to scratch if on njsw node (assuming we're running
        # jobs in parallel)
        #
        if os.path.exists(task_params['contig_file_path']):
            assembly =  task_params['contig_file_path']
            print("FOUND ASSEMBLY ON LOCAL SCRATCH")
        else:
            # we are on njsw so lets copy it over to scratch
            assembly = self.get_contig_file(task_params['assembly_ref'])

        # remove spaces from fasta headers because that breaks bedtools
        assembly_clean = os.path.basename(assembly).split('.fa')[0] + "_clean.fa"

        command = '/bin/bash reformat.sh in={} out={} addunderscore'.format(assembly,assembly_clean)

        log('running reformat command: {}'.format(command))
        out,err = self.run_command(command)

        #min_contig_length = task_params['min_contig_length']

        fastq = read_scratch_path[0]
        result_directory = task_params['result_directory']
        sam = os.path.basename(fastq) + '.sam'
        sam = os.path.join(result_directory, sam)
        sorted_bam = sam + '.sorted.bam'

        command = '/bin/bash bbmap.sh -Xmx{} fast threads={} ref={} in={} out={} mappedonly nodisk overwrite'.format(self.BBMAP_MEM,self.BBMAP_THREADS,assembly_clean,fastq,sam)
        print("Test PRINT statement")

        log('running alignment command: {}'.format(command))
        out,err = self.run_command(command)

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

    def generate_concoct_command(self, params):
        """
        generate_command: generate concoct params
        """
        print("\n\nRunning generate_concoct_command")

        command = 'python /kb/deployment/bin/CONCOCT/scripts/cut_up_fasta.py {} -c {} -o 0 --merge_last -b temp.bed > concoct_output_dir/split_contigs.fa'.format(params.get('contig_file_path'),params.get('contig_split_size'))
        command += ' && '
        command += 'python /kb/deployment/bin/CONCOCT/scripts/concoct_coverage_table.py temp.bed concoct_output_dir/*.sorted.bam > concoct_output_dir/coverage_table.tsv'
        command += ' && '
        command += '/kb/deployment/bin/CONCOCT/bin/concoct --composition_file concoct_output_dir/split_contigs.fa -l {} -b concoct_output_dir --coverage_file concoct_output_dir/coverage_table.tsv -t {}'.format(params.get('min_contig_length'),
                                         self.BBMAP_THREADS)

        command += ' && '
        command += 'python /kb/deployment/bin/CONCOCT/scripts/merge_cutup_clustering.py concoct_output_dir/clustering_gt{}.csv > concoct_output_dir/clustering_merged.csv'.format(params.get('min_contig_length'))
        command += ' && '
        command += 'mkdir concoct_output_dir/final_bins'
        command += ' && '
        command += 'python /kb/deployment/bin/CONCOCT/scripts/extract_fasta_bins.py {} concoct_output_dir/clustering_merged.csv --output_path concoct_output_dir/final_bins'.format(params.get('contig_file_path'))


        log('Generated run_concoct command: {}'.format(command))

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
            for root, dirs, files in os.walk(result_directory):
                for file in files:
                    if not (file.endswith('.fasta') or file.endswith('.DS_Store')):
                        zip_file.write(os.path.join(root, file), file)
                    if file.endswith('.marker.pdf'):
                        report_file = os.path.join(root, file)

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'File(s) generated by CONCOCT App'})

        if report_file:
            output_files.append({'path': report_file,
                                 'name': os.path.basename(report_file),
                                 'label': os.path.basename(report_file),
                                 'description': 'Visualization of the marker by CONCOCT App'})

        return output_files

    def generate_html_report(self, result_directory, assembly_ref, binned_contig_obj_ref, header):
        """
        generate_html_report: generate html summary report
        """

        log('Start generating html report')
        html_report = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self.mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'report.html')

        Overview_Content = ''

        (binned_contig_count, input_contig_count,
         too_short_count, total_bins_count, total_binned_contig_len,
         totoal_contig_length, total_short_contig_len) = self.generate_overview_info(
                                                            assembly_ref,
                                                            binned_contig_obj_ref,
                                                            result_directory)
        Overview_Content += '<p>Bins: {}</p>'.format(total_bins_count)
        Overview_Content += '<p>Input Contigs: {}</p>'.format(input_contig_count)
        Overview_Content += '<p>Binned Contigs: {} ({:.1%})</p>'.format(
                                            binned_contig_count,
                                            binned_contig_count/float(input_contig_count))
        Overview_Content += '<p>Unbinned Contigs: {} ({:.1%})</p>'.format(
                                            input_contig_count - binned_contig_count,
                                            1 - binned_contig_count/float(input_contig_count))
        Overview_Content += '<p>Contigs Too Short: {} ({:.1%})</p>'.format(
                                            too_short_count,
                                            too_short_count/float(input_contig_count))
        Overview_Content += '<p>Summed Length of Binned Contigs: {} ({:.1%})</p>'.format(
                                            total_binned_contig_len,
                                            total_binned_contig_len/float(totoal_contig_length))
        Overview_Content += '<p>Summed Length of Unbinned Contigs: {} ({:.1%})</p>'.format(
                                            totoal_contig_length - total_binned_contig_len,
                                            1 - total_binned_contig_len/float(totoal_contig_length))
        Overview_Content += '<p>Summed Length of Short Contigs: {} ({:.1%})</p>'.format(
                                            total_short_contig_len,
                                            total_short_contig_len/float(totoal_contig_length))

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          Overview_Content)
                result_file.write(report_template)

        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for CONCOCT App'})
        return html_report

    def generate_overview_info(self, assembly_ref, binned_contig_obj_ref, result_directory):
        """
        generate_overview_info: generate overview information from assembly and binnedcontig
        """

        assembly = self.dfu.get_objects({'object_refs': [assembly_ref]})['data'][0]
        binned_contig = self.dfu.get_objects({'object_refs': [binned_contig_obj_ref]})['data'][0]

        input_contig_count = assembly.get('data').get('num_contigs')

        totoal_contig_length = 0
        for contig_id, contig in assembly.get('data').get('contigs').items():
            totoal_contig_length += int(contig.get('length'))

        binned_contig_count = 0
        total_bins = binned_contig.get('data').get('bins')
        total_binned_contig_len = binned_contig.get('data').get('total_contig_len')
        total_bins_count = len(total_bins)
        for bin in total_bins:
            binned_contig_count += len(bin.get('contigs'))

        too_short_count = 0
        total_short_contig_len = 0
        result_files = os.listdir(result_directory)
        for file_name in result_files:
            if file_name.endswith('.tooshort'):
                for record in SeqIO.parse(os.path.join(result_directory, file_name), "fasta"):
                    total_short_contig_len += len(str(record.seq))
                    too_short_count += 1

        return (binned_contig_count, input_contig_count,
                too_short_count, total_bins_count, total_binned_contig_len,
                totoal_contig_length, total_short_contig_len)

    def generate_report(self, binned_contig_obj_ref, result_directory, params):
        """
        generate_report: generate summary report

        """
        log('Generating report')

        output_files = self.generate_output_file_list(result_directory)

        output_html_files = self.generate_html_report(result_directory,
                                                       params.get('assembly_ref'),
                                                       binned_contig_obj_ref,
                                                       params.get('out_header'))

        created_objects = []
        created_objects.append({"ref": binned_contig_obj_ref,
                                "description": "BinnedContigs from CONCOCT"})

        report_params = {
              'message': '',
              'workspace_name': params.get('workspace_name'),
              'objects_created': created_objects,
              'file_links': output_files,
              'html_links': output_html_files,
              'direct_html_link_index': 0,
              'html_window_height': 266,
              'report_object_name': 'kb_concoct_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

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
        params['out_header'] = 'Bin'

        contig_file = self.get_contig_file(params.get('assembly_ref'))
        params['contig_file_path'] = contig_file

        reads_list_file = self.stage_reads_list_file(params.get('reads_list'))
        params['reads_list_file'] = reads_list_file

        result_directory = os.path.join(self.scratch, str("concoct_output_dir"))

        params['result_directory'] = result_directory

        self.mkdir_p(result_directory)

        cwd = os.getcwd()
        log('changing working dir to {}'.format(result_directory))
        os.chdir(result_directory)
        #run alignments
        self.generate_alignment_bams(params)

        #run concoct
        command = self.generate_concoct_command(params)

        self.run_command(command)

        os.chdir(cwd)
        log('changing working dir to {}'.format(cwd))

        log('Saved result files to: {}'.format(result_directory))
        log('Generated files:\n{}'.format('\n'.join(os.listdir(result_directory))))

        generate_binned_contig_param = {
            'file_directory': result_directory,
            'assembly_ref': params.get('assembly_ref'),
            'binned_contig_name': params.get('binned_contig_name'),
            'workspace_name': params.get('workspace_name')
        }

        binned_contig_obj_ref = self.mgu.file_to_binned_contigs(
                                    generate_binned_contig_param).get('binned_contig_obj_ref')

        reportVal = self.generate_report(binned_contig_obj_ref, result_directory, params)

        returnVal = {
            'result_directory': result_directory,
            'binned_contig_obj_ref': binned_contig_obj_ref
        }

        returnVal.update(reportVal)

        return returnVal
