#!/usr/bin/env python3
import sys
import json
from ttp import ttp

if 'snakemake' in globals():
    input = snakemake.input
    output = snakemake.output
else:
    input = sys.argv[1]
    output = sys.argv[2]


template = """
<group name="{{ file_name }}">

<group name="header">
#= Summary for dataset: {{ file_name }} 
#     Query mRNAs :      {{ query_mrna }} in      {{ query_loci }} loci  ({{ query_multiexon }} multi-exon transcripts)
#            ({{ query_muti_transcript_loci }} multi-transcript loci, {{ transcripts_per_locus }} transcripts per locus)
# Reference mRNAs :      {{ ref_mrna }} in      {{ ref_loci }} loci  ({{ ref_mutiexon }} multi-exon)
# Super-loci w/ reference transcripts:       {{ super_loci }}

</group>
<group name="table">
{{ ignore("\s+") }}{{ level | _start_ }} level:    {{ sensitivity }}     |    {{ precision }}    |{{ _end_ }}
</group>

<group name="matching">
{{ ignore("\s+") }}Matching {{level | ORPHRASE}}:      {{number}}\n
</group>

<group name="missed">
{{ ignore("\s+") }}Missed {{level}}:{{ ignore("\s+") }}{{ num }}/{{ den }}\t({{ ignore("\s+") }}{{ perc }}%)
</group>

<group name="novel">
{{ ignore("\s+") }}Novel {{level}}:{{ ignore("\s+") }}{{ num }}/{{ den }}\t({{ ignore("\s+") }}{{ perc }}%)
</group>


</group>
"""

try:
    with open(input) as fin:
        data = fin.read()
except TypeError:
    print(
        f'Please the script as: \n    {sys.argv[0]} input output \n or use within snakemake')
    sys.exit(1)


parser = ttp(data=data, template=template)
parser.parse()
res = parser.result(structure="json")

with open(output, 'w') as out_file:
    json.dump(res, out_file)
