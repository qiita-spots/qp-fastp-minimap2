# -----------------------------------------------------------------------------
# Copyright (c) 2020--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from unittest import main
from qiita_client.testing import PluginTestCase
from os import remove, environ
from os.path import exists, isdir, join, dirname
from shutil import rmtree, copyfile
from tempfile import mkdtemp
from json import dumps
from itertools import zip_longest
from functools import partial

from qp_fastp_bwamem import plugin
from qp_fastp_bwamem.utils import plugin_details
from qp_fastp_bwamem.qp_fastp_bwamem import (
    get_references_list, _generate_commands, fastp_bwamem_to_array,
    QC_REFERENCE, FASTP_CMD, COMBINED_CMD, FASTP_CMD_SINGLE,
    COMBINED_CMD_SINGLE)


class FastpBwamemTests(PluginTestCase):
    def setUp(self):
        """Connects to the server and sets up the plugin"""
        plugin("https://localhost:21174", 'register', 'ignored')

        out_dir = mkdtemp()
        self.maxDiff = None
        self.out_dir = out_dir
        self.refs = get_references_list()
        self.db_path = QC_REFERENCE
        self.params = {'reference': 'genome', 'threads': 2}
        self._clean_up_files = []
        self._clean_up_files.append(out_dir)

    def tearDown(self):
        """Deletes all the files/directories created"""
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_get_references_list(self):
        """Tests that get_references_list returns just genome.fasta"""
        refs = get_references_list()
        self.assertCountEqual(refs, ['genome.fasta'])

    def test_generate_commands(self):
        """Tests that _generate_commands function returns correct commands"""

        # test parameters
        params = {'reference': 'genome', 'nprocs': 2,
                  'out_dir': '/foo/bar/output'}

        fwd_seqs = ['sz1.fastq.gz', 'sc1.fastq.gz',
                    'sa1.fastq.gz', 'sd1.fastq.gz']
        rev_seqs = ['sz2.fastq.gz', 'sc2.fastq.gz',
                    'sa2.fastq.gz', 'sd2.fastq.gz']

        # runs generate_commands with test params
        # testing for genome reference with rev_seqs
        obs = _generate_commands(fwd_seqs, rev_seqs, params['reference'],
                                 params['nprocs'], params['out_dir'])
        # formats the command to generate "truths"
        cmd = COMBINED_CMD.format(**params)

        # command truths (formats the commands with the fastq sequences)
        ecmds = [cmd % (f, r, f, r, f, r)
                 for f, r in zip_longest(fwd_seqs, rev_seqs)]
        # out_file truths (list of paths to output files in sorted order)
        eof = [(f'{params["out_dir"]}/{f}', 'raw_forward_seqs')
               for f in sorted(fwd_seqs)]
        # sorts rev_seqs and adds them as well to the eof
        for f in sorted(rev_seqs):
            eof.append((f'{params["out_dir"]}/{f}', 'raw_reverse_seqs'))
        # check for equality
        self.assertCountEqual(obs[0], ecmds)
        self.assertCountEqual(obs[1], eof)

        # modifies test params
        # testing for no reference with rev_seqs
        params['reference'] = None
        obs = _generate_commands(fwd_seqs, rev_seqs, params['reference'],
                                 params['nprocs'], params['out_dir'])
        cmd = FASTP_CMD.format(**params)
        ecmds = [cmd % (f, r, f, r)
                 for f, r in zip_longest(fwd_seqs, rev_seqs)]
        self.assertCountEqual(obs[0], ecmds)
        self.assertCountEqual(obs[1], list(eof))

        # tests the genome reference with no rev_seqs
        params['reference'] = 'genome'
        obs = _generate_commands(fwd_seqs, [], params['reference'],
                                 params['nprocs'], params['out_dir'])
        cmd = COMBINED_CMD_SINGLE.format(**params)
        ecmds = [cmd % (f, f, f) for f in fwd_seqs]
        eof = [(f'{params["out_dir"]}/{f}', 'raw_forward_seqs')
               for f in sorted(fwd_seqs)]
        self.assertCountEqual(obs[0], ecmds)
        self.assertCountEqual(obs[1], eof)

        # tests no reference with no rev_seqs
        params['reference'] = None
        obs = _generate_commands(fwd_seqs, [], params['reference'],
                                 params['nprocs'], params['out_dir'])
        cmd = FASTP_CMD_SINGLE.format(**params)
        ecmds = [cmd % (f, f) for f in fwd_seqs]
        self.assertCountEqual(obs[0], ecmds)
        self.assertCountEqual(obs[1], eof)

    def test_fastp_bwamem(self):
        """Tests that fastp_bwamem function returns correct job results"""

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
        source_dir = 'qp_fastp_bwamem/support_files/raw_data'
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
                'command': dumps([plugin_details['name'],
                                  plugin_details['version'],
                                  'Adapter and host filtering']),
                'status': 'running',
                'parameters': dumps(self.params)}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # adding extra parameters
        self.params['environment'] = environ["ENVIRONMENT"]

        # Get the artifact filepath information
        artifact_info = self.qclient.get("/qiita_db/artifacts/%s/" % aid)

        # Get the artifact metadata
        prep_info = self.qclient.get('/qiita_db/prep_template/%s/' % pid)
        prep_file = prep_info['prep-file']

        url = 'this-is-my-url'

        main_qsub_fp, finish_qsub_fp, out_files_fp = fastp_bwamem_to_array(
            artifact_info['files'], out_dir, self.params, prep_file,
            url, job_id)

        od = partial(join, out_dir)
        self.assertEqual(od(f'{job_id}.qsub'), main_qsub_fp)
        self.assertEqual(od(f'{job_id}.finish.qsub'), finish_qsub_fp)
        self.assertEqual(od(f'{job_id}.out_files.tsv'), out_files_fp)

        with open(main_qsub_fp) as f:
            main_qsub = f.readlines()
        with open(finish_qsub_fp) as f:
            finish_qsub = f.readlines()
        with open(out_files_fp) as f:
            out_files = f.readlines()
        with open(f'{out_dir}/fastp_bwamem.array-details') as f:
            commands = f.readlines()

        exp_main_qsub = [
            '#!/bin/bash\n',
            '#PBS -M qiita.help@gmail.com\n',
            f'#PBS -N {job_id}\n',
            '#PBS -l nodes=1:ppn=2\n',
            '#PBS -l walltime=30:00:00\n',
            '#PBS -l mem=16g\n',
            f'#PBS -o {out_dir}/{job_id}_${{PBS_ARRAYID}}.log\n',
            f'#PBS -e {out_dir}/{job_id}_${{PBS_ARRAYID}}.err\n',
            '#PBS -t 1-2%8\n',
            '#PBS -l epilogue=/home/qiita/qiita-epilogue.sh\n',
            'set -e\n',
            f'cd {out_dir}\n',
            'source ~/.bash_profile; source activate qp-fastp-bwamem; '
            f'export QC_REFERENCE={QC_REFERENCE}\n',
            'date\n',
            'hostname\n',
            'echo ${PBS_JOBID} ${PBS_ARRAYID}\n',
            'offset=${PBS_ARRAYID}\n', 'step=$(( $offset - 0 ))\n',
            f'cmd=$(head -n $step {out_dir}/fastp_bwamem.array-details | '
            'tail -n 1)\n',
            'eval $cmd\n',
            'set +e\n',
            'date\n']
        self.assertEqual(main_qsub, exp_main_qsub)

        exp_finish_qsub = [
            '#!/bin/bash\n',
            '#PBS -M qiita.help@gmail.com\n',
            f'#PBS -N finish-{job_id}\n',
            '#PBS -l nodes=1:ppn=1\n',
            '#PBS -l walltime=10:00:00\n',
            '#PBS -l mem=10g\n',
            f'#PBS -o {out_dir}/finish-{job_id}.log\n',
            f'#PBS -e {out_dir}/finish-{job_id}.err\n',
            '#PBS -l epilogue=/home/qiita/qiita-epilogue.sh\n',
            'set -e\n',
            f'cd {out_dir}\n',
            'source ~/.bash_profile; source activate qp-fastp-bwamem; '
            f'export QC_REFERENCE={QC_REFERENCE}\n',
            'date\n',
            'hostname\n',
            'echo $PBS_JOBID\n',
            f'finish_qp_fastp_bwamem this-is-my-url {job_id} {out_dir}\n',
            'date\n']
        self.assertEqual(finish_qsub, exp_finish_qsub)

        exp_out_files = [
            f'{out_dir}/S22205_S104_L001_R1_001.fastq.gz\traw_forward_seqs\n',
            f'{out_dir}/S22205_S104_L001_R2_001.fastq.gz\traw_reverse_seqs\n',
            f'{out_dir}/S22282_S102_L001_R1_001.fastq.gz\traw_forward_seqs\n',
            f'{out_dir}/S22282_S102_L001_R2_001.fastq.gz\traw_reverse_seqs']
        self.assertEqual(out_files, exp_out_files)

        # the easiest to figure out the location of the artifact input files
        # is to check the first file of the raw forward reads
        apath = dirname(artifact_info['files']['raw_forward_seqs'][0])
        exp_commands = [
            f'fastp -l 100 -i {apath}/S22205_S104_L001_R1_001.fastq.gz -w 2  '
            f'-I {apath}/S22205_S104_L001_R2_001.fastq.gz --stdout | '
            f'bwa mem -t 2 {QC_REFERENCE}genome.fasta  '
            f'{apath}/S22205_S104_L001_R1_001.fastq.gz '
            f'{apath}/S22205_S104_L001_R2_001.fastq.gz | '
            'samtools fastq -@ 2 -f  12 -F 256 -1 '
            f'{out_dir}/S22205_S104_L001_R1_001.fastq.gz -2 '
            f'{out_dir}/S22205_S104_L001_R2_001.fastq.gz\n',
            f'fastp -l 100 -i {apath}/S22282_S102_L001_R1_001.fastq.gz -w 2  '
            f'-I {apath}/S22282_S102_L001_R2_001.fastq.gz --stdout | '
            f'bwa mem -t 2 {QC_REFERENCE}genome.fasta  '
            f'{apath}/S22282_S102_L001_R1_001.fastq.gz '
            f'{apath}/S22282_S102_L001_R2_001.fastq.gz | '
            'samtools fastq -@ 2 -f  12 -F 256 -1 '
            f'{out_dir}/S22282_S102_L001_R1_001.fastq.gz -2 '
            f'{out_dir}/S22282_S102_L001_R2_001.fastq.gz']
        self.assertEqual(commands, exp_commands)

    def test_fastp_bwamem_just_fwd(self):
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
        fp2_1 = join(in_dir, 'S22282_S102_L001_R1_001.fastq.gz')
        source_dir = 'qp_fastp_bwamem/support_files/raw_data'
        copyfile(f'{source_dir}/S22205_S104_L001_R1_001.fastq.gz', fp1_1)
        copyfile(f'{source_dir}/S22282_S102_L001_R1_001.fastq.gz', fp2_1)

        data = {
            'filepaths': dumps([
                (fp1_1, 'raw_forward_seqs'),
                (fp2_1, 'raw_forward_seqs')]),
            'type': "per_sample_FASTQ",
            'name': "Test artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['input'] = aid

        data = {'user': 'demo@microbio.me',
                'command': dumps([plugin_details['name'],
                                  plugin_details['version'],
                                  'Adapter and host filtering']),
                'status': 'running',
                'parameters': dumps(self.params)}
        job_id = self.qclient.post(
            '/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # adding extra parameters
        self.params['environment'] = environ["ENVIRONMENT"]

        # Get the artifact filepath information
        artifact_info = self.qclient.get("/qiita_db/artifacts/%s/" % aid)

        # Get the artifact metadata
        prep_info = self.qclient.get('/qiita_db/prep_template/%s/' % pid)
        prep_file = prep_info['prep-file']

        url = 'this-is-my-url'

        main_qsub_fp, finish_qsub_fp, out_files_fp = fastp_bwamem_to_array(
            artifact_info['files'], out_dir, self.params, prep_file,
            url, job_id)

        od = partial(join, out_dir)
        self.assertEqual(od(f'{job_id}.qsub'), main_qsub_fp)
        self.assertEqual(od(f'{job_id}.finish.qsub'), finish_qsub_fp)
        self.assertEqual(od(f'{job_id}.out_files.tsv'), out_files_fp)

        with open(main_qsub_fp) as f:
            main_qsub = f.readlines()
        with open(finish_qsub_fp) as f:
            finish_qsub = f.readlines()
        with open(out_files_fp) as f:
            out_files = f.readlines()
        with open(f'{out_dir}/fastp_bwamem.array-details') as f:
            commands = f.readlines()

        exp_main_qsub = [
            '#!/bin/bash\n',
            '#PBS -M qiita.help@gmail.com\n',
            f'#PBS -N {job_id}\n',
            '#PBS -l nodes=1:ppn=2\n',
            '#PBS -l walltime=30:00:00\n',
            '#PBS -l mem=16g\n',
            f'#PBS -o {out_dir}/{job_id}_${{PBS_ARRAYID}}.log\n',
            f'#PBS -e {out_dir}/{job_id}_${{PBS_ARRAYID}}.err\n',
            '#PBS -t 1-2%8\n',
            '#PBS -l epilogue=/home/qiita/qiita-epilogue.sh\n',
            'set -e\n',
            f'cd {out_dir}\n',
            'source ~/.bash_profile; source activate qp-fastp-bwamem; '
            f'export QC_REFERENCE={QC_REFERENCE}\n',
            'date\n',
            'hostname\n',
            'echo ${PBS_JOBID} ${PBS_ARRAYID}\n',
            'offset=${PBS_ARRAYID}\n', 'step=$(( $offset - 0 ))\n',
            f'cmd=$(head -n $step {out_dir}/fastp_bwamem.array-details | '
            'tail -n 1)\n',
            'eval $cmd\n',
            'set +e\n',
            'date\n']
        self.assertEqual(main_qsub, exp_main_qsub)

        exp_finish_qsub = [
            '#!/bin/bash\n',
            '#PBS -M qiita.help@gmail.com\n',
            f'#PBS -N finish-{job_id}\n',
            '#PBS -l nodes=1:ppn=1\n',
            '#PBS -l walltime=10:00:00\n',
            '#PBS -l mem=10g\n',
            f'#PBS -o {out_dir}/finish-{job_id}.log\n',
            f'#PBS -e {out_dir}/finish-{job_id}.err\n',
            '#PBS -l epilogue=/home/qiita/qiita-epilogue.sh\n',
            'set -e\n',
            f'cd {out_dir}\n',
            'source ~/.bash_profile; source activate qp-fastp-bwamem; '
            f'export QC_REFERENCE={QC_REFERENCE}\n',
            'date\n',
            'hostname\n',
            'echo $PBS_JOBID\n',
            f'finish_qp_fastp_bwamem this-is-my-url {job_id} {out_dir}\n',
            'date\n']
        self.assertEqual(finish_qsub, exp_finish_qsub)

        exp_out_files = [
            f'{out_dir}/S22205_S104_L001_R1_001.fastq.gz\traw_forward_seqs\n',
            f'{out_dir}/S22282_S102_L001_R1_001.fastq.gz\traw_forward_seqs']
        self.assertEqual(out_files, exp_out_files)

        # the easiest to figure out the location of the artifact input files
        # is to check the first file of the raw forward reads
        apath = dirname(artifact_info['files']['raw_forward_seqs'][0])
        exp_commands = [
            f'fastp -l 100 -i {apath}/S22205_S104_L001_R1_001.fastq.gz -w 2  '
            f'--stdout | bwa mem -t 2 {QC_REFERENCE}genome.fasta  '
            f'{apath}/S22205_S104_L001_R1_001.fastq.gz | '
            'samtools fastq -@ 2 -f  4 -0 '
            f'{out_dir}/S22205_S104_L001_R1_001.fastq.gz\n',
            f'fastp -l 100 -i {apath}/S22282_S102_L001_R1_001.fastq.gz -w 2  '
            f'--stdout | bwa mem -t 2 {QC_REFERENCE}genome.fasta  '
            f'{apath}/S22282_S102_L001_R1_001.fastq.gz | '
            'samtools fastq -@ 2 -f  4 -0 '
            f'{out_dir}/S22282_S102_L001_R1_001.fastq.gz']
        self.assertEqual(commands, exp_commands)


if __name__ == '__main__':
    main()
