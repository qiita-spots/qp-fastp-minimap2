# -----------------------------------------------------------------------------
# Copyright (c) 2020, Qiita development team.
#
# Distributed under the terms of the BSD 3-clause License License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand
from .qp_minimap2 import get_ref_list, minimap2
from .utils import plugin_details
from os.path import splitext


THREADS = 15

# comment
# Initialize the plugin
plugin = QiitaPlugin(**plugin_details)

# Define the command
dbs = get_ref_list()
dbs_without_extension = [splitext(db)[0] for db in dbs]
dbs_defaults = ', '.join([f'"{x}"' for x in dbs_without_extension])
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    'reference': [
        f'choice:["None", {dbs_defaults}]', dbs_without_extension[0]],
    'threads': ['integer', f'{THREADS}']}

outputs = {'Filtered files': 'per_sample_FASTQ'}
default_params = {
    'auto-detect adapters only filtering [not recommended]': {
        'reference': "None", 'threads': THREADS}}
for db in dbs_without_extension:
    name = f'auto-detect adapters and {db} + phix filtering'
    default_params[name] = {'reference': db, 'threads': THREADS}

minimap2_cmd = QiitaCommand(
    'Adapter and host filtering', "Sequence adapter and host filtering",
    minimap2, req_params, opt_params, outputs, default_params)

plugin.register_command(minimap2_cmd)
