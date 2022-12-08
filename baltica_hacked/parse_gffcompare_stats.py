#!/usr/bin/env python3
import sys
import pkgutil
from ttp import ttp

try:
    if 'snakemake' in globals():
        input_ = snakemake.input[0]
        output = snakemake.output[0]
    else:
        input_ = sys.argv[1]
        output = sys.argv[2]

    with open(input_) as fin:
        data = fin.read()
except (TypeError, IndexError):
    print(
        f'Please use the script as: \n    {sys.argv[0]} input_ output \n or use it within snakemake')
    sys.exit(1)

template = pkgutil.get_data(__name__, "resources/gff_compare_template.txt")

parser = ttp(data=data, template=template)
parser.parse()
res = parser.result(format="json")[0]

with open(output, 'w') as out_file:
    out_file.writelines(res)
