comp = config["contrasts"].keys()
workdir: config["path"]
config_path = config["config_path"]


subworkflow qc:
    workdir:
        config["path"]
    snakefile:
        "qc.smk"
    configfile: 
        config_path

subworkflow junctionseq:
    workdir:
        config["path"]
    snakefile:
        "junctionseq.smk"
    configfile: 
        config_path

subworkflow leafcutter:
    workdir:
        config["path"]
    snakefile:
        "leafcutter.smk"
    configfile: 
        config_path

subworkflow majiq:
    workdir:
        config["path"]
    snakefile:
        "majiq.smk"
    configfile: 
        config_path

subworkflow stringtie:
    workdir:
        config["path"]
    snakefile:
        "stringtie.smk"
    configfile: 
        config_path

subworkflow rmats:
    workdir:
        config["path"]
    snakefile:
        "rmats.smk"
    configfile: 
        config_path

subworkflow analysis:
    workdir:
        config["path"]
    snakefile:
        "analysis.smk"
    configfile: 
        config_path

rule final:
    input:
        qc("qc/multiqc/multiqc_report.html"),
        majiq(
            expand("majiq/voila/{comp}_voila.tsv",comp=comp)),
        rmats(
            expand("rmats/{comp}/",comp=comp)),
        leafcutter(
            expand("leafcutter/{comp}/{comp}_cluster_significance.txt",comp=comp)),
        junctionseq(
            expand("junctionseq/analysis/{comp}_sigGenes.results.txt.gz", comp=comp)),
        stringtie("stringtie/merged/merged.combined.gtf"),
        analysis(
            expand(
            "results/baltica_report{project_title}.html", 
            project_title="_" + config.get("project_title", "").replace(' ', '_'))),
        analysis(
            expand(
            "results/baltica_table{project_title}.xlsx", 
            project_title="_" + config.get("project_title", "").replace(' ', '_')))
