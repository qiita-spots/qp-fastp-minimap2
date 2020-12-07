# -----------------------------------------------------------------------------
# Copyright (c) 2020, Qiita development team.
#
# Distributed under the terms of the BSD 3-clause License License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand
from .qp_fastp_minimap2 import get_dbs_list, fastp_minimap2

# TODO: include relevant imports here
# from .mycommand import my_command_function

# Initialize the plugin
plugin = QiitaPlugin('qp-fastp-minimap2', '2021.01',
                     'fastp + minimap2 pipeline')

# Define the command
default_db_list = get_dbs_list()
req_params = {'input': ('artifact', ['per_sample_FASTQ'])}
opt_params = {
    'adapter_1': ["choice: ['AAAA', 'BBBB', 'CCCC']", 'AAAA'],
    'adapter_2': ["choice: ['ZZZZ', 'YYYY', 'XXXX']", 'ZZZZ'],
    'reference': ["choice: [%s]" % default_db_list, default_db_list[0]],
    'threads': 15}
outputs = {'Filtered files': 'per_sample_FASTQ'}
default_params = {'adapter_1': 'AAAA', 'adapter_2': 'ZZZZ',
                  'reference': default_db_list[0], 'threads': 15}

fastp_minimap2_cmd = QiitaCommand(
    'Adapter and host filtering', "Sequence adapter and host filtering",
    fastp_minimap2, req_params, opt_params, outputs, default_params)

plugin.register_command(fastp_minimap2_cmd)
