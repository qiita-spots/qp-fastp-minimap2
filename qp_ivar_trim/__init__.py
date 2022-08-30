# -----------------------------------------------------------------------------
# Copyright (c) 2020, Qiita development team.
#
# Distributed under the terms of the BSD 3-clause License License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand
from .qp_ivar_trim import get_primer, ivar_trim
from .utils import plugin_details
from os.path import splitext


PRIMER_POS_OFFSET = 5


# Initialize the plugin
plugin = QiitaPlugin(**plugin_details)

# Define the command
primer = get_primer()
primer_without_extension = splitext(primer)[0]
req_params = {'input': ('artifact', ['BAM'])}
opt_params = {
    'primer': [
        'string', primer_without_extension],
    'primer_offset': ['integer', f'{PRIMER_POS_OFFSET}']}

outputs = {'Trimmed files': 'BAM'}
# default_params = {
#     'auto-detect adapters only filtering [not recommended]': {
#         'primer': "None", 'threads': THREADS}}
# for db in dbs_without_extension:
#     name = f'auto-detect adapters and {db} + phix filtering'
#     default_params[name] = {'primer': db, 'threads': THREADS}
default_params = {
    'default params': {
        'primer': primer_without_extension,
        'primer_offset': 5}}

ivar_trim_cmd = QiitaCommand(
    'ivar trim', "Trimming with ivar",
    ivar_trim, req_params, opt_params, outputs, default_params)

plugin.register_command(ivar_trim_cmd)
