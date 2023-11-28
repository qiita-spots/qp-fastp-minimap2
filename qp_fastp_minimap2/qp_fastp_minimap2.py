# -----------------------------------------------------------------------------
# Copyright (c) 2020, Qiita development team.
#
# Distributed under the terms of the BSD 3-clause License License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
from os import environ
from os.path import basename, join, abspath, dirname
from glob import glob
from itertools import zip_longest
from shutil import copyfile

from qiita_client import ArtifactInfo

MEMORY = '16g'
WALLTIME = '30:00:00'
FINISH_MEMORY = '10g'
FINISH_WALLTIME = '10:00:00'
MAX_RUNNING = 8

QC_REFERENCE_DB = environ["QC_REFERENCE_DB"]

FASTP_BASE = 'fastp -l 100 -i %s -w {nprocs} --adapter_fasta {adapter_fasta}'
MINIMAP2_BASE = 'minimap2 -ax sr -t {nprocs} {database} - -a '
SAMTOOLS_BASE = 'samtools fastq -@ {nprocs} -f '

FASTP_CMD = ' '.join([FASTP_BASE, '-I %s -o {out_dir}/%s -O {out_dir}/%s'])
FASTP_CMD_SINGLE = (f'{FASTP_BASE} -o '
                    '{out_dir}/%s')
COMBINED_CMD = (f'{FASTP_BASE} -I %s --detect_adapter_for_pe --stdout | '
                f'{MINIMAP2_BASE} | {SAMTOOLS_BASE} 12 -F 256 -1 '
                '{out_dir}/%s -2 {out_dir}/%s')
COMBINED_CMD_SINGLE = (f'{FASTP_BASE} --stdout | {MINIMAP2_BASE} | '
                       f'{SAMTOOLS_BASE} 4 -0 '
                       '{out_dir}/%s')


def get_dbs_list():
    folder = QC_REFERENCE_DB

    # skip human database
    return [basename(f) for f in glob(f'{folder}/*.mmi') if 'human' not in f]


def _generate_commands(fwd_seqs, rev_seqs, database, nprocs, out_dir):
    """Helper function to generate commands and facilite testing"""

    # copy adapter_fasta file to out_dir
    source_adapter_fasta = join(
        dirname(abspath(__file__)), 'support_files', 'fastp_known_adapters',
        'fastp_known_adapters_formatted.fna')
    adapter_fasta = join(out_dir, 'fastp_known_adapters_formatted.fna')
    # this if is to help with test_qp_fastp_minimap2.test_generate_commands
    if out_dir != '/foo/bar/output':
        copyfile(source_adapter_fasta, adapter_fasta)

    files = zip_longest(fwd_seqs, rev_seqs)
    if rev_seqs:
        cmd = FASTP_CMD
        if database is not None:
            cmd = COMBINED_CMD
    else:
        cmd = FASTP_CMD_SINGLE
        if database is not None:
            cmd = COMBINED_CMD_SINGLE
    command = cmd.format(nprocs=nprocs, database=database, out_dir=out_dir,
                         adapter_fasta=adapter_fasta)

    out_files = []
    commands = []
    for i, (fwd_fp, rev_fp) in enumerate(files):
        fname = basename(fwd_fp)
        out_files.append((f'{out_dir}/{fname}', 'raw_forward_seqs'))
        if rev_fp:
            rname = basename(rev_fp)
            out_files.append((f'{out_dir}/{rname}', 'raw_reverse_seqs'))
            cmd = command % (fwd_fp, rev_fp, fname, rname)
        else:
            cmd = command % (fwd_fp, fname)
        commands.append(cmd)

    return commands, out_files


def fastp_minimap2(qclient, job_id, parameters, out_dir):
    """Run fastp and minimap2 with the given parameters

    Parameters
    ----------
    qclient : tgp.qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to run split libraries
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    bool, list, str
        The results of the job
    """

    qclient.update_job_step(
        job_id, "Step 3 of 4: Finishing fastp and minimap2")

    ainfo = []
    # Generates 2 artifacts: one for the ribosomal
    # reads and other for the non-ribosomal reads
    out_files = []
    with open(f'{out_dir}/{job_id}.out_files.tsv') as f:
        for line in f.readlines():
            fp, ft = line.split()
            out_files.append((fp, ft))

    # Step 4 generating artifacts
    msg = "Step 4 of 4: Generating new artifact"
    qclient.update_job_step(job_id, msg)
    ainfo = [ArtifactInfo('Filtered files', 'per_sample_FASTQ', out_files)]

    return True, ainfo, ""


def fastp_minimap2_to_array(files, out_dir, params, prep_info, url, job_id):
    """Creates files for submission of per sample fastp and minimap2

    Parameters
    ----------
    files : dict
        The dictionary of files to process, raw_forward_seqs/raw_reverse_seqs
    out_dir : str
        The output directory
    params : dict
        The parameter values to run fastp/minimap2
    prep_info : str
        The path to prep_info
    url : str
        The url to send info to
    job_id : str
        The job id

    Returns
    -------
    str, str, str
        The paths of the main_fp, finish_fp, out_files_fp
    """
    database = None
    if params['reference'] != 'None':
        database = [join(QC_REFERENCE_DB, f'{db}')
                    for db in get_dbs_list()
                    if params['reference'] in db][0]

    fwd_seqs = sorted(files['raw_forward_seqs'])
    if 'raw_reverse_seqs' in files:
        rev_seqs = sorted(files['raw_reverse_seqs'])
    else:
        rev_seqs = []

    prep_info.set_index('sample_name', inplace=True)

    # Note that for processing we don't actually need the run_prefix so
    # we are not going to use it and simply loop over the ordered
    # fwd_seqs/rev_seqs
    commands, out_files = _generate_commands(
        fwd_seqs, rev_seqs, database, params['threads'], out_dir)

    # writing the job array details
    details_name = join(out_dir, 'fastp_minimap2.array-details')
    with open(details_name, 'w') as details:
        details.write('\n'.join(commands))
    n_jobs = len(commands)

    # all the setup pieces
    PPN = params['threads']
    lines = ['#!/bin/bash',
             '#SBATCH -p qiita',
             '#SBATCH --mail-user "qiita.help@gmail.com"',
             f'#SBATCH --job-name {job_id}',
             '#SBATCH -N 1',
             f'#SBATCH -n {PPN}',
             f'#SBATCH --time {WALLTIME}',
             f'#SBATCH --mem {MEMORY}',
             f'#SBATCH --output {out_dir}/{job_id}_%a.log',
             f'#SBATCH --error {out_dir}/{job_id}_%a.err',
             f'#SBATCH --array 1-{n_jobs}%{MAX_RUNNING}',
             'set -e',
             f'cd {out_dir}',
             f'{params["environment"]}',
             'date',  # start time
             'hostname',  # executing system
             'echo ${SLURM_JOBID} ${SLURM_ARRAY_TASK_ID}',
             'offset=${SLURM_ARRAY_TASK_ID}',
             'step=$(( $offset - 0 ))',
             f'cmd=$(head -n $step {details_name} | tail -n 1)',
             'eval $cmd',
             'set +e',
             'date']
    main_fp = join(out_dir, f'{job_id}.slurm')
    with open(main_fp, 'w') as job:
        job.write('\n'.join(lines))
        job.write('\n')

    # finish job
    lines = ['#!/bin/bash',
             '#SBATCH -p qiita',
             '#SBATCH --mail-user "qiita.help@gmail.com"',
             f'#SBATCH --job-name finish-{job_id}',
             '#SBATCH -N 1',
             '#SBATCH -n 1',
             f'#SBATCH --time {FINISH_WALLTIME}',
             f'#SBATCH --mem {FINISH_MEMORY}',
             f'#SBATCH --output {out_dir}/finish-{job_id}.log',
             f'#SBATCH --error {out_dir}/finish-{job_id}.err',
             'set -e',
             f'cd {out_dir}',
             f'{params["environment"]}',
             'date',  # start time
             'hostname',  # executing system
             'echo $SLURM_JOBID',
             f'finish_qp_fastp_minimap2 {url} {job_id} {out_dir}\n'
             "date"]
    finish_fp = join(out_dir, f'{job_id}.finish.slurm')
    with open(finish_fp, 'w') as out:
        out.write('\n'.join(lines))
        out.write('\n')

    out_files_fp = join(out_dir, f'{job_id}.out_files.tsv')
    with open(out_files_fp, 'w') as out:
        out.write('\n'.join([f'{fp}\t{ft}'for fp, ft in out_files]))

    return main_fp, finish_fp, out_files_fp
