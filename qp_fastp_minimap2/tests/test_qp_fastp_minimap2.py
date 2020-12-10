# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import main
from qiita_client.testing import PluginTestCase
from qiita_client import ArtifactInfo
from os import remove
from os.path import exists, isdir, join
from shutil import rmtree, copyfile
from tempfile import mkdtemp
from json import dumps

from qp_fastp_minimap2 import plugin
from qp_fastp_minimap2.qp_fastp_minimap2 import (
    get_dbs_list, fastp_minimap2, QC_REFERENCE_DB)


class FastpMinimap2Tests(PluginTestCase):
    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')

        out_dir = mkdtemp()
        self.maxDiff = None
        self.out_dir = out_dir
        self.dbs = get_dbs_list()
        self.db_path = QC_REFERENCE_DB
        self.params = {'reference': 'artifacts', 'threads': 2}
        self._clean_up_files = []
        self._clean_up_files.append(out_dir)

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_get_dbs_list(self):
        dbs = get_dbs_list()
        self.assertCountEqual(dbs, ['artifacts.mmi', 'empty.mmi'])

    def test_fastp_minimap2(self):
        # inserting new prep template
        prep_info_dict = {
            'SKB8.640193': {'run_prefix': 'S22205_S104'},
            'SKD8.640184': {'run_prefix': 'S22282_S102'}}
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': 'Metagenomic'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        in_dir = mkdtemp()
        self._clean_up_files.append(in_dir)

        fp1_1 = join(in_dir, 'S22205_S104_L001_R1_001.fastq.gz')
        fp1_2 = join(in_dir, 'S22205_S104_L001_R2_001.fastq.gz')
        fp2_1 = join(in_dir, 'S22282_S102_L001_R1_001.fastq.gz')
        fp2_2 = join(in_dir, 'S22282_S102_L001_R2_001.fastq.gz')
        source_dir = 'qp_fastp_minimap2/support_files/raw_data'
        copyfile(f'{source_dir}/S22205_S104_L001_R1_001.fastq.gz', fp1_1)
        copyfile(f'{source_dir}/S22205_S104_L001_R2_001.fastq.gz', fp1_2)
        copyfile(f'{source_dir}/S22282_S102_L001_R1_001.fastq.gz', fp2_1)
        copyfile(f'{source_dir}/S22282_S102_L001_R2_001.fastq.gz', fp2_2)

        data = {
            'filepaths': dumps([
                (fp1_1, 'raw_forward_seqs'),
                (fp1_2, 'raw_reverse_seqs'),
                (fp2_1, 'raw_forward_seqs'),
                (fp2_2, 'raw_reverse_seqs')]),
            'type': "per_sample_FASTQ",
            'name': "Test artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid

        data = {'user': 'demo@microbio.me',
                'command': dumps(['qp-fastp-minimap2', '2021.01',
                                  'Adapter and host filtering']),
                'status': 'running',
                'parameters': dumps(self.params)}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = fastp_minimap2(
            self.qclient, job_id, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        files = [(f'{out_dir}/S22205_S104_L001_R1_001.fastq.gz',
                  'raw_forward_seqs'),
                 (f'{out_dir}/S22205_S104_L001_R2_001.fastq.gz',
                  'raw_reverse_seqs'),
                 (f'{out_dir}/S22282_S102_L001_R1_001.fastq.gz',
                  'raw_forward_seqs'),
                 (f'{out_dir}/S22282_S102_L001_R2_001.fastq.gz',
                  'raw_reverse_seqs')]
        exp = ArtifactInfo('Filtered files', 'per_sample_FASTQ', files)

        self.assertEqual(ainfo, exp)


if __name__ == '__main__':
    main()
