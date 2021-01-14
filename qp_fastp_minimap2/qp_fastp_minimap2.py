# -----------------------------------------------------------------------------
# Copyright (c) 2020, Qiita development team.
#
# Distributed under the terms of the BSD 3-clause License License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------
import pandas as pd

from os import environ
from os.path import basename, join
from glob import glob
from itertools import zip_longest

from qiita_client import ArtifactInfo
from qiita_client.util import system_call


QC_REFERENCE_DB = environ["QC_REFERENCE_DB"]


FASTP_BASE = 'fastp -l 100 -i %s -w {nprocs} '
MINIMAP2_BASE = 'minimap2 -ax sr -t {nprocs} {database} - -a '
SAMTOOLS_BASE = 'samtools fastq -@ {nprocs} -f 12 -F 256'

FASTP_CMD = ' '.join([FASTP_BASE, '-I %s -o {out_dir}/%s -O {out_dir}/%s'])
FASTP_CMD_SINGLE = (f'{FASTP_BASE} -o '
                    '{out_dir}/%s')
COMBINED_CMD = (f'{FASTP_BASE} -I %s --stdout | {MINIMAP2_BASE} | '
                f'{SAMTOOLS_BASE} -1 '
                '{out_dir}/%s -2 {out_dir}/%s')
COMBINED_CMD_SINGLE = (f'{FASTP_BASE} --stdout | {MINIMAP2_BASE} | '
                       f'{SAMTOOLS_BASE} -1 '
                       '{out_dir}/%s')


def get_dbs_list():
    folder = QC_REFERENCE_DB

    # skip human database
    return [basename(f) for f in glob(f'{folder}/*.mmi') if 'human' not in f]


def _generate_commands(fwd_seqs, rev_seqs, database, nprocs, out_dir):
    """Helper function to generate commands and facilite testing"""
    files = zip_longest(fwd_seqs, rev_seqs)
    if rev_seqs:
        cmd = FASTP_CMD
        if database is not None:
            cmd = COMBINED_CMD
    else:
        cmd = FASTP_CMD_SINGLE
        if database is not None:
            cmd = COMBINED_CMD_SINGLE
    command = cmd.format(nprocs=nprocs, database=database, out_dir=out_dir)

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
    # Step 1 get the rest of the information need to run Bowtie2
    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = parameters['input']
    del parameters['input']

    # Get the artifact filepath information
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    fwd_seqs = sorted(artifact_info['files']['raw_forward_seqs'])
    if 'raw_reverse_seqs' in artifact_info['files']:
        rev_seqs = sorted(artifact_info['files']['raw_reverse_seqs'])
    else:
        rev_seqs = []

    # Get the artifact metadata
    prep_info = qclient.get('/qiita_db/prep_template/%s/'
                            % artifact_info['prep_information'][0])
    df = pd.read_csv(prep_info['prep-file'], sep='\t', dtype='str',
                     na_values=[], keep_default_na=True)
    df.set_index('sample_name', inplace=True)
    if 'run_prefix' not in df.columns:
        raise ValueError('Missing run_prefix column in your preparation')

    database = None
    if parameters['reference'] != 'None':
        database = [join(QC_REFERENCE_DB, f'{db}')
                    for db in get_dbs_list()
                    if parameters['reference'] in db][0]

    # Step 2 generating command
    qclient.update_job_step(job_id, "Step 2 of 4: Generating commands")

    # Note that for processing we don't actually need the run_prefix so
    # we are not going to use it and simply loop over the ordered
    # fwd_seqs/rev_seqs
    commands, out_files = _generate_commands(
        fwd_seqs, rev_seqs, database, parameters['threads'], out_dir)
    len_commands = len(commands) + 1

    msg = f"Step 3 of 4: Executing QC_Filter job (%d/{len_commands})"
    for i, cmd in enumerate(commands):
        qclient.update_job_step(job_id, msg % (i+1))
        std_out, std_err, return_value = system_call(cmd)
        if return_value != 0:
            error_msg = (f"Error running: {cmd}\n"
                         f"Std out: {std_out}\nStd err: {std_err}\n")
            return False, None, error_msg

    # Step 4 generating artifacts
    msg = "Step 4 of 4: Generating new artifact"
    qclient.update_job_step(job_id, msg)
    ainfo = ArtifactInfo('Filtered files', 'per_sample_FASTQ', out_files)

    return True, ainfo, ""
