import pandas as pd

configfile: "config.yaml"
metadata = pd.read_csv(config["samplemetadata"])

#seperate out normal bam file from tumor cell bams
metadata_normal = metadata[metadata["normal"] == "T"]
metadata = metadata[metadata["normal"] == "F"]
samples = metadata["sample"].unique()

rule all:
    input:
        expand("results/{sample}/hmmcopy_results/metrics.csv.gz", sample = samples),
        expand("results/{sample}/hmmcopy_results/reads.csv.gz", sample = samples),
        expand("results/{sample}/counthaps/allele_counts_perblock.csv.gz", sample = samples),
        expand("results/{sample}/schnapps/schnapps.Rdata", sample = samples),
        expand("results/{sample}/counthaps/allele_counts_all.csv.gz", sample = samples),

include: "rules/hmmcopy.smk"
include: "rules/haplotyping.smk"
include: "rules/ascn.smk"