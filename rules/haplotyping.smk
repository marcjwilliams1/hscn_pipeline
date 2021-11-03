def getnormalbam(wildcards):
    x = metadata_normal[metadata_normal["sample"] == wildcards.sample]
    return x["bamfiles"]

rule create_yaml_normal_inferhaps:
    input: getnormalbam
    output: temp("results/{sample}/infer_haps.yaml")
    shell:
        """
        printf "normal:\n  bam: {input[0]}" > {output[0]}
        """

rule inferhaps:
    input:
        yaml = "results/{sample}/infer_haps.yaml",
    output:
        inferhapsdir = directory("results/{sample}/inferhaps/"),
        haplotypes = "results/{sample}/inferhaps/haplotypes.csv.gz",
        pipelinedir = temp(directory("results/{sample}/inferhaps/pipeline/")),
    params:
        config = "config.yaml",
    threads: 40
    resources:
        mem_mb=1024 * 4
    singularity: "docker://quay.io/singlecellpipeline/single_cell_pipeline_haplotypes:v0.8.9"
    shell:
        """
        single_cell infer_haps \
            --input_yaml {input.yaml} \
            --maxjobs {threads} \
            --nocleanup \
            --config_file {params.config} \
            --sentinel_only \
            --submit local \
            --loglevel DEBUG \
            --pipelinedir {output.pipelinedir} \
            --out_dir {output.inferhapsdir}
        """

rule createvcf:
    input: haplotypes = "results/{sample}/inferhaps/haplotypes.csv.gz",
    output: haplotypesvcf = "results/{sample}/inferhaps/haplotypes.vcf",
    threads: 1
    resources: mem_mb=1024 * 50
    singularity: "docker://marcjwilliams1/schnapps"
    script: "../scripts/create_vcf_file.R"

rule genotypecells:
    input:
        haplotypesvcf = "results/{sample}/inferhaps/haplotypes.vcf",
        bam = getbamfile
    params:
        results_dir = directory("results/{sample}/counthaps/{cell_id}/"),
    output:
        results_DP_mtx = temp("results/{sample}/counthaps/{cell_id}/cellSNP.tag.DP.mtx"),
        results_AD_mtx = temp("results/{sample}/counthaps/{cell_id}/cellSNP.tag.AD.mtx"),
        results_OTH_mtx = temp("results/{sample}/counthaps/{cell_id}/cellSNP.tag.OTH.mtx"),
        results_sample = temp("results/{sample}/counthaps/{cell_id}/cellSNP.samples.tsv"),
        vcf = temp("results/{sample}/counthaps/{cell_id}/cellSNP.base.vcf.gz"),
    conda: "../envs/cellsnp.yml"
    threads: 10
    resources: mem_mb=1024 * 1
    shell:
        """
        cellsnp-lite -s {input.bam} \
            -I {wildcards.cell_id} \
            -O {params.results_dir} \
            -R {input.haplotypesvcf} \
            -p {threads} \
            --cellTAG None \
            --UMItag None \
            --gzip \
            --minCOUNT 0 \
            --exclFLAG UNMAP,SECONDARY,QCFAIL,DUP
        """

rule format_per_cell_counts:
    input: 
        vcf = "results/{sample}/counthaps/{cell_id}/cellSNP.base.vcf.gz",
        results_DP_mtx = "results/{sample}/counthaps/{cell_id}/cellSNP.tag.DP.mtx",
        results_AD_mtx = "results/{sample}/counthaps/{cell_id}/cellSNP.tag.AD.mtx",
        results_OTH_mtx = "results/{sample}/counthaps/{cell_id}/cellSNP.tag.OTH.mtx",
        results_sample = "results/{sample}/counthaps/{cell_id}/cellSNP.samples.tsv",
        haplotypes = "results/{sample}/inferhaps/haplotypes.csv.gz"
    output: 
        alldata = temp("results/{sample}/counthaps/{cell_id}/allele_counts_all.csv.gz"),
        perblock = temp("results/{sample}/counthaps/{cell_id}/allele_counts_perblock.csv.gz")
    threads: 5
    resources: mem_mb=1024 * 5
    singularity: "docker://marcjwilliams1/schnapps"
    script: "../scripts/format_vcf_haps.R"

def allsnpsfiles(wildcards):
    sample = metadata[metadata['sample'] == wildcards.sample]
    x = expand("results/{{sample}}/counthaps/{cell_id}/allele_counts_all.csv.gz", cell_id = sample["cell_id"])
    return x

rule mergesnps:
    input: 
        allsnpsfiles
    output:
        merged = "results/{sample}/counthaps/allele_counts_all.csv.gz"
    threads: 4
    resources: mem_mb=1024 * 25
    run:
        dist_list = []
        for f in input:
            print(f)
            dist_temp = pd.read_csv(f)
            dist_list.append(dist_temp)
        dist = pd.concat(dist_list, ignore_index=True)
        dist.to_csv(output.merged, sep=',')

def allblockfiles(wildcards):
    sample = metadata[metadata['sample'] == wildcards.sample]
    x = expand("results/{{sample}}/counthaps/{cell_id}/allele_counts_perblock.csv.gz", cell_id = sample["cell_id"])
    return x

rule mergeblocks:
    input: 
        allblockfiles
    output:
        merged = "results/{sample}/counthaps/allele_counts_perblock.csv.gz"
    threads: 4
    resources: mem_mb=1024 * 24
    run:
        dist_list = []
        for f in input:
            print(f)
            dist_temp = pd.read_csv(f)
            dist_list.append(dist_temp)
        dist = pd.concat(dist_list, ignore_index=True)
        dist.to_csv(output.merged, sep=',')