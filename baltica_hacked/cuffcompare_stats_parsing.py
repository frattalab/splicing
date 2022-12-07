from ttp import ttp
import snakemake

template = """
<group name="header">
#= Summary for dataset: {{file_name}} :
# {{ ignore("\s+") }}Query mRNAs :  {{query_mrna}} in   {{query_loci}} loci  ({{query_multiexon}} multi-exon transcripts)
# {{ ignore("\s+") }}({{query_muti_transcript_loci}} multi-transcript loci, ~{{transcripts_per_locus}} transcripts per locus)
# Reference mRNAs :  {{ref_mrna}} in   {{ref_loci}} loci  ({{ref_mutiexon}} multi-exon)
# Super-loci w/ reference transcripts:    {{super_loci}}
</group>
<group name="table">
{{ ignore("\s+") }}{{level}} level:{{ ignore("\s+") }}{{sn}}{{ ignore("\s+") }}{{sp}}{{ ignore("\s+") }}{{fsn}}{{ ignore("\s+") }}{{fsp}}
</group>
{{ ignore("\s+") }}Missed exons:{{ ignore("\s+") }}{{ a }}/ {{ b}}
<group name="summary">
{{ ignore("\s+") }}{{ summary | PHRASE | _line_ | contains("intron", "exon", "loci")  }}
</group>
"""

with open(snakemake.input) as fin:
    data = fin.read()
parser = ttp(data=data, template=template)
parser.parse()
res = parser.result(structure="json")
with open(snakemake.output) as fout:
    fout.write(parser.result(format="json"))
