# -----------------------------------------------------------------------------
# Copyright (c) 2020, Qiita development team.
#
# Distributed under the terms of the BSD 3-clause License License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand
from .qp_ivar_trim import get_dbs_list, ivar_trim
from .utils import plugin_details
from os.path import splitext


THREADS = 15


# Initialize the plugin
plugin = QiitaPlugin(**plugin_details)

# Define the command
dbs = get_dbs_list()
dbs_without_extension = [splitext(db)[0] for db in dbs]
dbs_defaults = ', '.join([f'"{x}"' for x in dbs_without_extension])
req_params = {'input': ('artifact', ['BAM'])}
opt_params = {
    'primer': [
        f'choice:["None", {dbs_defaults}]', dbs_without_extension[0]],
    'threads': ['integer', f'{THREADS}']}

outputs = {'trimmed files': 'BAM'}
default_params = {
    'auto-detect adapters only filtering [not recommended]': {
        'primer': "None", 'threads': THREADS}}
for db in dbs_without_extension:
    name = f'auto-detect adapters and {db} + phix filtering'
    default_params[name] = {'primer': db, 'threads': THREADS}

ivar_trim_cmd = QiitaCommand(
    'Adapter and host filtering', "Sequence adapter and host filtering",
    ivar_trim, req_params, opt_params, outputs, default_params)

plugin.register_command(ivar_trim_cmd)
