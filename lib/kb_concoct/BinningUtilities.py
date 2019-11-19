import os
import sys
import time


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class BinningUtil:
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
        print('bbstats_output {}'.format(bbstats_output))
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
        gc_avg = bbstats_output.split('\t')[17] * 100
        gc_std = bbstats_output.split('\t')[18] * 100

        log('Generated generate_stats_for_genome_bins command: {}'.format(command))

        return {'n_scaffolds': n_scaffolds, 'n_contigs': n_contigs, 'scaf_bp': scaf_bp, 'contig_bp': contig_bp, 'gap_pct': gap_pct, 'scaf_N50': scaf_N50, 'scaf_L50': scaf_L50, 'ctg_N50': ctg_N50, 'scaf_N90': scaf_N90, 'ctg_N90': ctg_N90, 'ctg_L90': ctg_L90, 'scaf_max': scaf_max, 'ctg_max': ctg_max, 'scaf_n_gt50K': scaf_n_gt50K, 'gc_avg': gc_avg, 'gc_std':gc_std}



    def run_read_mapping_unpaired_mode(self, task_params, assembly_clean, fastq, sam):
        if task_params['read_mapping_tool'] == 'bbmap':
            command = '/bin/bash bbmap.sh -Xmx{} fast threads={} ref={} in={} out={} mappedonly nodisk overwrite'.format(self.BBMAP_MEM, self.MAPPING_THREADS, assembly_clean, fastq, sam)
        elif task_params['read_mapping_tool'] == 'bwa':
            command = 'bwa index {} && bwa mem -t {} {} {} > {}'.format(assembly_clean, self.MAPPING_THREADS, assembly_clean, fastq, sam)
        elif task_params['read_mapping_tool'] == 'bowtie2':
            bt2index = os.path.basename(assembly_clean) + '.bt2'
            command = 'bowtie2-build -f {} --threads {} {} && bowtie2 -x {} -U {} --threads {} -S {}'.format(assembly_clean, self.MAPPING_THREADS, bt2index, bt2index, fastq, self.MAPPING_THREADS, sam)
        elif task_params['read_mapping_tool'] == 'minimap2':
            command = 'minimap2 -ax sr -t {} {} {} > {}'.format(self.MAPPING_THREADS, assembly_clean, fastq, sam)
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
        path_to_concoct_result_bins = '{}/{}'.format(path_to_concoct_result, self.CONCOCT_BIN_RESULT_DIR)
        path_to_concoct_result_bins = path_to_concoct_result_bins + '/' # need to fix
        path_to_summary_file = '{}notreal.summary'.format(path_to_concoct_result_bins) # need to fix
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
