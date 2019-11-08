# -*- coding: utf-8 -*-
import os
import time
import shutil
import unittest
from configparser import ConfigParser

from kb_concoct.kb_concoctImpl import kb_concoct
from kb_concoct.kb_concoctServer import MethodContext
from kb_concoct.authclient import KBaseAuth as _KBaseAuth

from kb_concoct.Utils.ConcoctUtil import ConcoctUtil


from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.WorkspaceClient import Workspace


class kb_concoctTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        #token = os.environ.get('KB_AUTH_TOKEN', None)
        cls.token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_concoct'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(cls.token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_concoct',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_concoct(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_kb_concoct_" + str(suffix)
        # ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

        cls.ws_info = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa
        cls.dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'], token=cls.token)
        cls.ru = ReadsUtils(os.environ['SDK_CALLBACK_URL'], token=cls.token)
        cls.au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'], token=cls.token)
        cls.concoct_runner = ConcoctUtil(cls.cfg)
        cls.prepare_data()


    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa

    @classmethod
    def prepare_data(cls):
        """
        Lets put everything on workspace
        """
        #
        # READS 1
        # building Interleaved library
        pe1_reads_filename = 'lib1.oldstyle.fastq'
        pe1_reads_path = os.path.join(cls.scratch, pe1_reads_filename)

        # gets put on scratch. "work/tmp" is scratch
        shutil.copy(os.path.join("data", pe1_reads_filename), pe1_reads_path)

        int1_reads_params = {
            'fwd_file': pe1_reads_path,
            'sequencing_tech': 'Unknown',
            'wsname': cls.ws_info[1],
            'name': 'MyInterleavedLibrary1',
            'interleaved': 'true'
        }

        #from scratch upload to workspace

        cls.int1_oldstyle_reads_ref = cls.ru.upload_reads(int1_reads_params)['obj_ref']
        print("cls.int1_oldstyle_reads_ref***")
        print(str(cls.int1_oldstyle_reads_ref))
        # READS 2
        # building Interleaved library
        pe2_reads_filename = 'lib2.oldstyle.fastq'
        pe2_reads_path = os.path.join(cls.scratch, pe2_reads_filename)

        # gets put on scratch. "work/tmp" is scratch
        shutil.copy(os.path.join("data", pe2_reads_filename), pe2_reads_path)

        int2_reads_params = {
            'fwd_file': pe2_reads_path,
            'sequencing_tech': 'Unknown',
            'wsname': cls.ws_info[1],
            'name': 'MyInterleavedLibrary2',
            'interleaved': 'true'
        }

        #from scratch upload to workspace
        cls.int2_oldstyle_reads_ref = cls.ru.upload_reads(int2_reads_params)['obj_ref']
        #
        # building Assembly
        #
        assembly_filename = 'small_arctic_assembly.fa'
        cls.assembly_filename_path = os.path.join(cls.scratch, assembly_filename)
        shutil.copy(os.path.join("data", assembly_filename), cls.assembly_filename_path)

        # from scratch upload to workspace
        assembly_params = {
            'file': {'path': cls.assembly_filename_path},
            'workspace_name': cls.ws_info[1],
            'assembly_name': 'MyAssembly'
        }

        # puts assembly object onto shock
        cls.assembly_ref = cls.au.save_assembly_from_fasta(assembly_params)

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        return self.ws_info[1]

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx


    def test_bad_run_concoct_params(self):
        method_name = 'test_bad_run_concoct_params'
        print ("\n=================================================================")
        print ("RUNNING "+method_name+"()")
        print ("=================================================================\n")

        invalidate_input_params = {
          'missing_assembly_ref': 'assembly_ref',
          'binned_contig_name': 'binned_contig_name',
          'workspace_name': 'workspace_name',
          'reads_list': 'reads_list'
        }
        with self.assertRaisesRegex(
                    ValueError, '"assembly_ref" parameter is required, but missing'):
            self.getImpl().run_kb_concoct(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'assembly_ref': 'assembly_ref',
          'missing_binned_contig_name': 'binned_contig_name',
          'workspace_name': 'workspace_name',
          'reads_list': 'reads_list'
        }
        with self.assertRaisesRegex(
                    ValueError, '"binned_contig_name" parameter is required, but missing'):
            self.getImpl().run_kb_concoct(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'assembly_ref': 'assembly_ref',
          'binned_contig_name': 'binned_contig_name',
          'missing_workspace_name': 'workspace_name',
          'reads_list': 'reads_list'
        }
        with self.assertRaisesRegex(
                    ValueError, '"workspace_name" parameter is required, but missing'):
            self.getImpl().run_kb_concoct(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
          'assembly_ref': 'assembly_ref',
          'binned_contig_name': 'binned_contig_name',
          'workspace_name': 'workspace_name',
          'missing_reads_list': 'reads_list'
        }
        with self.assertRaisesRegex(
                    ValueError, '"reads_list" parameter is required, but missing'):
            self.getImpl().run_kb_concoct(self.getContext(), invalidate_input_params)


    def test_ConcoctUtil_generate_command(self):
        method_name = 'test_ConcoctUtil_generate_command'
        print ("\n=================================================================")
        print ("RUNNING "+method_name+"()")
        print ("=================================================================\n")

        input_params = {
            'contig_file_path': 'mycontig',
            'min_contig_length': '3000',
            'contig_split_size': '10000',
            'contig_split_overlap': '0',
            'kmer_size': '4',
            'max_clusters_for_vgmm': '400',
            'max_iterations_for_vgmm': '500',
            'total_percentage_pca': '90',
            'no_cov_normalization': '--no_cov_normalization',
            'no_total_coverage': '--no_total_coverage'
        }

        expect_command = 'python /kb/deployment/bin/CONCOCT/scripts/cut_up_fasta.py mycontig -c 10000 -o 0 --merge_last '
        expect_command += '-b temp.bed > concoct_output_dir/split_contigs.fa && python /kb/deployment/bin/CONCOCT/scripts/concoct_coverage_table.py '
        expect_command += 'temp.bed concoct_output_dir/*.sorted.bam > concoct_output_dir/coverage_table.tsv && python /kb/deployment/bin/CONCOCT/bin/concoct '
        expect_command += '--composition_file concoct_output_dir/split_contigs.fa -l 3000 -b concoct_output_dir '
        expect_command += '--coverage_file concoct_output_dir/coverage_table.tsv -t 16 -k 4 -c 400 -i 500 --total_percentage_pca 90 --no_cov_normalization --no_total_coverage && python '
        expect_command += '/kb/deployment/bin/CONCOCT/scripts/merge_cutup_clustering.py concoct_output_dir/clustering_gt3000.csv '
        expect_command += '> concoct_output_dir/clustering_merged.csv && mkdir '
        expect_command += 'concoct_output_dir/final_bins && python /kb/deployment/bin/CONCOCT/scripts/extract_fasta_bins.py '
        expect_command += 'mycontig concoct_output_dir/clustering_merged.csv '
        expect_command += '--output_path concoct_output_dir/final_bins'
        command = self.concoct_runner.generate_concoct_command(input_params)
        self.assertEqual(command, expect_command)

    def test_run_concoct(self):
        method_name = 'test_run_concoct'
        print ("\n=================================================================")
        print ("RUNNING "+method_name+"()")
        print ("=================================================================\n")

        # concoct should run to completion here
        ret = self.getImpl().run_kb_concoct(self.getContext(),
                                            {'workspace_name': self.getWsName(),
                                             'assembly_ref': self.assembly_ref,
                                             'min_contig_length': 3000,
                                             'contig_split_size': 10000,
                                             'contig_split_overlap': 0,
                                             'kmer_size': 4,
                                             'max_clusters_for_vgmm': 400,
                                             'max_iterations_for_vgmm': 500,
                                             'total_percentage_pca': 90,
                                             'no_cov_normalization': '--no_cov_normalization',
                                             'no_total_coverage': '--no_total_coverage',
                                             'binned_contig_name': 'concoct_bin_obj',
                                             'reads_list': [self.int1_oldstyle_reads_ref, self.int2_oldstyle_reads_ref] })







    # def test_ConcoctUtil_generate_output_file_list(self):
    #     method_name = 'test_ConcoctUtil_generate_output_file_list'
    #     print ("\n=================================================================")
    #     print ("RUNNING "+method_name+"()")
    #     print ("=================================================================\n")

    #     result_file_file_directory = 'Concoct_result'
    #     result_file_path = os.path.join(self.scratch, result_file_file_directory)
    #     if not os.path.exists(result_file_path):
    #         os.makedirs(result_file_path)

    #     for item in os.listdir(os.path.join("data", "Concoct_Result_Sample")):
    #         shutil.copy(os.path.join("data", "Concoct_Result_Sample", item),
    #                     os.path.join(result_file_path, item))

    #     output_file = self.concoct_runner._generate_output_file_list(result_file_path)[0]

    #     self.assertTrue('path' in output_file)
    #     output_file_path = output_file.get('path')
    #     expect_file_set = {'out_header.abund1', 'out_header.abund2', 'out_header.abund3',
    #                        'out_header.marker', 'out_header.marker_of_each_bin.tar.gz',
    #                        'out_header.summary', 'out_header.noclass', 'out_header.log',
    #                        'out_header.abundance', 'out_header.tooshort'}

    #     with zipfile.ZipFile(output_file_path) as z:
    #         self.assertEqual(set(z.namelist()), expect_file_set)

    #     self.assertEqual(output_file.get('description'), 'File(s) generated by Concoct App')
    #     self.assertEqual(output_file.get('name'), 'concoct_result.zip')
    #     self.assertEqual(output_file.get('label'), 'concoct_result.zip')

    # def test_ConcoctUtil_stage_reads_list_file(self):
    #     method_name = 'test_ConcoctUtil_stage_reads_list_file'
    #     print ("\n=================================================================")
    #     print ("RUNNING "+method_name+"()")
    #     print ("=================================================================\n")

    #     # test SingleEndLibrary
    #     reads_list = [self.se_reads_ref, self.se_reads_ref]

    #     reads_list_file = self.concoct_runner._stage_reads_list_file(reads_list)

    #     with open(reads_list_file) as file:
    #         result_file_list = file.readlines()

    #     self.assertEqual(len(result_file_list), len(reads_list))
    #     for item in result_file_list:
    #         self.assertRegex(item, r'.*\.single\.fastq.*')

    #     # test KBaseAssembly SingleEndLibrary
    #     reads_list = [self.KBA_se_reads_ref, self.KBA_se_reads_ref]

    #     reads_list_file = self.concoct_runner._stage_reads_list_file(reads_list)

    #     with open(reads_list_file) as file:
    #         result_file_list = file.readlines()

    #     self.assertEqual(len(result_file_list), len(reads_list))
    #     for item in result_file_list:
    #         self.assertRegex(item, r'.*\.single\.fastq.*')

    #     # test PairedEndLibrary
    #     reads_list = [self.pe_reads_ref, self.pe_reads_ref]

    #     reads_list_file = self.concoct_runner._stage_reads_list_file(reads_list)

    #     with open(reads_list_file) as file:
    #         result_file_list = file.readlines()

    #     self.assertEqual(len(result_file_list), len(reads_list))
    #     for item in result_file_list:
    #         self.assertRegex(item, r'.*\.inter\.fastq.*')

    #     # test KBaseAssembly PairedEndLibrary
    #     reads_list = [self.KBA_pe_reads_ref, self.KBA_pe_reads_ref]

    #     reads_list_file = self.concoct_runner._stage_reads_list_file(reads_list)

    #     with open(reads_list_file) as file:
    #         result_file_list = file.readlines()

    #     self.assertEqual(len(result_file_list), len(reads_list))
    #     for item in result_file_list:
    #         self.assertRegex(item, r'.*\.inter\.fastq.*')

    # def test_ConcoctUtil_get_contig_file(self):
    #     method_name = 'test_ConcoctUtil_get_contig_file'
    #     print ("\n=================================================================")
    #     print ("RUNNING "+method_name+"()")
    #     print ("=================================================================\n")

    #     contig_file = self.concoct_runner._get_contig_file(self.assembly_ref)

    #     with open(contig_file, 'r') as file:
    #         contig_file_content = file.readlines()

    #     with open(self.assembly_fasta_file_path, 'r') as file:
    #         expect_contig_file_content = file.readlines()

    #     self.assertCountEqual(contig_file_content, expect_contig_file_content)

    # def test_run_concoct_single_reads(self):
    #     method_name = 'test_run_concoct_single_reads'
    #     print ("\n=================================================================")
    #     print ("RUNNING "+method_name+"()")
    #     print ("=================================================================\n")

    #     input_params = {
    #         'assembly_ref': self.assembly_ref,
    #         'binned_contig_name': 'out_header'+'single',
    #         'workspace_name': self.getWsName(),
    #         'reads_list': [self.pe_reads_ref],
    #         'thread': 1,
    #         'prob_threshold': 0.5,
    #         'markerset': 170,
    #         'min_contig_length': 2000,
    #         'plotmarker': 1
    #     }

    #     result = self.getImpl().run_kb_concoct(self.getContext(), input_params)[0]

    #     self.assertTrue('result_directory' in result)
    #     self.assertTrue('binned_contig_obj_ref' in result)
    #     self.assertTrue('report_name' in result)
    #     self.assertTrue('report_ref' in result)

    # def test_run_concoct_multi_reads(self):
    #     method_name = 'test_run_concoct_multi_reads'
    #     print ("\n=================================================================")
    #     print ("RUNNING "+method_name+"()")
    #     print ("=================================================================\n")

    #     input_params = {
    #         'assembly_ref': self.assembly_ref,
    #         'binned_contig_name': 'out_header'+'multi',
    #         'workspace_name': self.getWsName(),
    #         'reads_list': [self.pe_reads_ref, self.pe_reads_ref, self.se_reads_ref],
    #         'thread': 4,
    #         'prob_threshold': 0.7,
    #         'markerset': 40,
    #         'min_contig_length': 1500,
    #         'plotmarker': 1
    #     }

    #     result = self.getImpl().run_kb_concoct(self.getContext(), input_params)[0]

    #     self.assertTrue('result_directory' in result)
    #     self.assertTrue('binned_contig_obj_ref' in result)
    #     self.assertTrue('report_name' in result)
    #     self.assertTrue('report_ref' in result)

    # def test_run_concoct_single_reads_KBaseAssembly_reads(self):
    #     method_name = 'test_run_concoct_single_reads_KBaseAssembly_reads'
    #     print ("\n=================================================================")
    #     print ("RUNNING "+method_name+"()")
    #     print ("=================================================================\n")

    #     input_params = {
    #         'assembly_ref': self.assembly_ref,
    #         'binned_contig_name': 'out_header'+'single'+'KBA',
    #         'workspace_name': self.getWsName(),
    #         'reads_list': [self.KBA_pe_reads_ref],
    #         'thread': 1,
    #         'prob_threshold': 0.5,
    #         'markerset': 170,
    #         'min_contig_length': 2000,
    #         'plotmarker': 1
    #     }

    #     result = self.getImpl().run_kb_concoct(self.getContext(), input_params)[0]

    #     self.assertTrue('result_directory' in result)
    #     self.assertTrue('binned_contig_obj_ref' in result)
    #     self.assertTrue('report_name' in result)
    #     self.assertTrue('report_ref' in result)

    # def test_run_concoct_multi_reads_KBaseAssembly_reads(self):
    #     method_name = 'test_run_concoct_multi_reads_KBaseAssembly_reads'
    #     print ("\n=================================================================")
    #     print ("RUNNING "+method_name+"()")
    #     print ("=================================================================\n")

    #     input_params = {
    #         'assembly_ref': self.assembly_ref,
    #         'binned_contig_name': 'out_header'+'multi'+'KBA',
    #         'workspace_name': self.getWsName(),
    #         'reads_list': [self.KBA_pe_reads_ref, self.KBA_pe_reads_ref, self.KBA_se_reads_ref],
    #         'thread': 4,
    #         'prob_threshold': 0.7,
    #         'markerset': 40,
    #         'min_contig_length': 1500,
    #         'plotmarker': 1
    #     }

    #     result = self.getImpl().run_kb_concoct(self.getContext(), input_params)[0]

    #     self.assertTrue('result_directory' in result)
    #     self.assertTrue('binned_contig_obj_ref' in result)
    #     self.assertTrue('report_name' in result)
    #     self.assertTrue('report_ref' in result)
