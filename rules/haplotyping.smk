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
        #export LSF_SERVERDIR=/admin/lsfjuno/lsf/10.1/linux3.10-glibc2.17-x86_64/etc 
        #export PATH=/common/juno/OS7/10.1/linux3.10-glibc2.17-x86_64/bin:$PATH 
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

def get_bam_list(wildcards):
    x = metadata[metadata["sample"] == wildcards.sample]
    bamlist = ','.join(x["bamfiles"])
    return bamlist

def get_cellid_list(wildcards):
    x = metadata[metadata["sample"] == wildcards.sample]
    cellidlist = ','.join(x["cell_id"])
    return cellidlist

rule genotypecells:
    input:
        haplotypesvcf = "results/{sample}/inferhaps/haplotypes.vcf",
    params:
        bams = get_bam_list,
        cell_ids = get_cellid_list,
        results_dir = directory("results/{sample}/counthaps/"),
    output:
        results_DP_mtx = "results/{sample}/counthaps/cellSNP.tag.DP.mtx",
        results_AD_mtx = "results/{sample}/counthaps/cellSNP.tag.AD.mtx",
        results_OTH_mtx = "results/{sample}/counthaps/cellSNP.tag.OTH.mtx",
        results_sample = "results/{sample}/counthaps/cellSNP.samples.tsv",
        vcf = "results/{sample}/counthaps/cellSNP.base.vcf.gz",
    conda: "../envs/cellsnp.yml"
    threads: 50
    resources: mem_mb=1024 * 1
    shell:
        """
        cellsnp-lite -s {params.bams} \
            -I {params.cell_ids} \
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
        vcf = "results/{sample}/counthaps/cellSNP.base.vcf.gz",
        results_DP_mtx = "results/{sample}/counthaps/cellSNP.tag.DP.mtx",
        results_AD_mtx = "results/{sample}/counthaps/cellSNP.tag.AD.mtx",
        results_OTH_mtx = "results/{sample}/counthaps/cellSNP.tag.OTH.mtx",
        results_sample = "results/{sample}/counthaps/cellSNP.samples.tsv",
        haplotypes = "results/{sample}/inferhaps/haplotypes.csv.gz"
    output: 
        alldata = "results/{sample}/counthaps/allele_counts_all.csv.gz",
        perblock = "results/{sample}/counthaps/allele_counts_perblock.csv.gz"
    threads: 40
    resources: mem_mb=1024 * 1
    singularity: "docker://marcjwilliams1/schnapps"
    script: "../scripts/format_vcf_haps.R"
