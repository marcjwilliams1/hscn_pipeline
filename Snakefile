import pandas as pd

configfile: "config.yaml"
CELL_IDS, = glob_wildcards(config["bamdir"] + "{id}.bam")
metadata = pd.read_csv(config["samplemetadata"])

#seperate out normal bam file from tumor cells
metadata_normal = metadata[metadata["normal"] == "T"]
metadata = metadata[metadata["normal"] == "F"]
samples = metadata["sample"].unique()


rule all:
    input:
        expand("results/{sample}/hmmcopy_results/metrics.csv.gz", sample = samples),
        expand("results/{sample}/hmmcopy_results/reads.csv.gz", sample = samples),
        expand("results/{sample}/counthaps/allele_counts_perblock.csv.gz", sample = samples)
        #"results/inferhaps/haplotypes.csv.gz"

def getbamfile(wildcards):
    x = metadata[metadata["cell_id"] == wildcards.cell_id]
    x = x[x["sample"] == wildcards.sample]
    return x["bamfiles"]

rule read_counter:
    input: getbamfile
    output: 
        mydir = directory("results/{sample}/read_counts/{cell_id}/"),
        file1 = "results/{sample}/read_counts/{cell_id}/test.txt",
        file2 = "results/{sample}/read_counts/{cell_id}/wig.txt",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 30
    conda: "envs/scp.yml"
    params:
        exclude_list=config["hmmcopy"]["exclude_list"],
        bin_size=config["hmmcopy"]["bin_size"],
        mapq=config["hmmcopy"]["min_mqual"]
    shell:
        """
        mkdir -p {output[0]}
        python scripts/read_counter.py {input[0]} {output[2]} -w {params.bin_size} --exclude_list {params.exclude_list} --mapping_quality_threshold {params.mapq}
        touch {output[1]}
        """

rule correct_read_counts:
    input: "results/{sample}/read_counts/{cell_id}/wig.txt"
    output: "results/{sample}/corrected_read_counts/{cell_id}.csv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 30
    conda: "envs/scp.yml"
    params:
        gc_wig_file=config["hmmcopy"]["gc_wig_file"],
        map_wig_file=config["hmmcopy"]["map_wig_file"]
    shell:
        """
        python scripts/corrected_read_counts.py {params.gc_wig_file} {params.map_wig_file} {input[0]} {output[0]}
        """

rule runhmmcopy:
    input: "results/{sample}/corrected_read_counts/{cell_id}.csv"
    output: 
        celldir = temp(directory("results/{sample}/hmmcopy_results_temp/{cell_id}/")),
        reads = temp("results/{sample}/hmmcopy_results_temp/{cell_id}/0/reads.csv"),
        metrics = temp("results/{sample}/hmmcopy_results_temp/{cell_id}/0/metrics.csv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 30
    conda: "envs/scp.yml"
    params: 
        strength=config["hmmcopy"]["strength"],
        e=config["hmmcopy"]["e"],
        mu=config["hmmcopy"]["mu"],
        lambdap=config["hmmcopy"]["lambda"],
        nu=config["hmmcopy"]["nu"],
        kappa=config["hmmcopy"]["kappa"],
        m=config["hmmcopy"]["m"],
        eta=config["hmmcopy"]["eta"],
        g=config["hmmcopy"]["g"],
        s=config["hmmcopy"]["s"],
        multiplier=config["hmmcopy"]["multipliers"],
    shell:
        """
        mkdir {output[0]}
        Rscript scripts/run_hmmcopy.R --corrected_data={input[0]} \
                                      --outdir={output[0]} \
                                      --sample_id={wildcards.cell_id} \
                                      --param_str={params.strength} \
                                      --param_e={params.e} \
                                      --param_mu={params.mu} \
                                      --param_l={params.lambdap} \
                                      --param_nu={params.nu} \
                                      --param_k={params.kappa} \
                                      --param_m={params.m} \
                                      --param_eta={params.eta} \
                                      --param_g={params.g} \
                                      --param_s={params.s} \
                                      --param_multiplier={params.multiplier}
        """

rule mvhmmcopy:
    input:
        reads = "results/{sample}/hmmcopy_results_temp/{cell_id}/0/reads.csv",
        celldir = "results/{sample}/hmmcopy_results_temp/{cell_id}/",
    output:
        reads = "results/{sample}/hmmcopy_results/percell/{cell_id}_reads.csv",
    shell:
        """
        mv {input.reads} {output.reads}
        """

rule add_metrics:
    input: 
        metrics = "results/{sample}/hmmcopy_results_temp/{cell_id}/0/metrics.csv",
        insrt = "results/{sample}/stats/{cell_id}.isize.txt",
        flgstat = "results/{sample}/stats/{cell_id}.bam.flagstat",
        celldir = "results/{sample}/hmmcopy_results_temp/{cell_id}/"
    output:
        metrics = "results/{sample}/hmmcopy_results/percell/{cell_id}_metrics.csv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 30
    conda: "envs/scp.yml"
    script: "scripts/extend_qc.py"


rule insert_size:
    input:
        getbamfile
    output:
        txt="results/{sample}/stats/{cell_id}.isize.txt",
        pdf="results/{sample}/stats/{cell_id}.isize.pdf"
    params:
        # optional parameters (e.g. relax checks as below)
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024
    wrapper:
        "0.78.0/bio/picard/collectinsertsizemetrics"

rule alignment_summary:
    input:        
        ref=config["ref_genome"],
        bam=getbamfile
    output:
        "results/{sample}/stats/{cell_id}.summary.txt"
    params:
        # optional parameters (e.g. relax checks as below)
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024
    wrapper:
        "0.78.0/bio/picard/collectalignmentsummarymetrics"

rule samtools_flagstat:
    input:
        getbamfile
    output:
        "results/{sample}/stats/{cell_id}.bam.flagstat"
    wrapper:
        "0.78.0/bio/samtools/flagstat"


def allreadsfile(wildcards):
    sample = metadata[metadata['sample'] == wildcards.sample]
    x = expand("results/{{sample}}/hmmcopy_results/percell/{cell_id}_reads.csv", cell_id = sample["cell_id"])
    return x

rule mergereads:
    input: 
        allreadsfile
    output:
        merged = "results/{sample}/hmmcopy_results/reads.csv.gz"
    run:
        dist_list = []
        for f in input:
            print(f)
            dist_temp = pd.read_csv(f)
            dist_list.append(dist_temp)
        dist = pd.concat(dist_list, ignore_index=True)
        dist.to_csv(output.merged, sep=',')

def allmetricsfile(wildcards):
    sample = metadata[metadata['sample'] == wildcards.sample]
    x = expand("results/{{sample}}/hmmcopy_results/percell/{cell_id}_metrics.csv", cell_id = sample["cell_id"])
    return x


rule mergemetrics:
    input: 
        allmetricsfile
    output:
        merged = "results/{sample}/hmmcopy_results/metrics.csv.gz"
    run:
        dist_list = []
        for f in input:
            print(f)
            dist_temp = pd.read_csv(f)
            dist_list.append(dist_temp)
        dist = pd.concat(dist_list, ignore_index=True)
        dist.to_csv(output.merged, sep=',')

def getnormalbam(wildcards):
    x = metadata[metadata_normal["sample"] == wildcards.sample]
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
        singularity = config["infer_haps_singularity"]
    threads: 40
    resources:
        mem_mb=1024 * 4
    shell:
        """
        module load singularity
        export LSF_SERVERDIR=/admin/lsfjuno/lsf/10.1/linux3.10-glibc2.17-x86_64/etc 
        export PATH=/common/juno/OS7/10.1/linux3.10-glibc2.17-x86_64/bin:$PATH 
        singularity run --bind /work/shah/users/william1/projects/schnapps_input_pipeline --bind /juno --bind /admin --bind /common \
            {params.singularity}  \
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
    script: "scripts/create_vcf_file.R"

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
    conda: "envs/cellsnp.yml"
    threads: 40
    resources: mem_mb=1024 * 4
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
    resources: mem_mb=1024 * 4
    singularity: "docker://marcjwilliams1/schnapps"
    script: "scripts/format_vcf_haps.R"

        # cellsnp-lite -s testdata/SA1090-A96213A-R20-C06.bam,testdata/SA1090-A96213A-R20-C08.bam \
        #     -I SA1090-A96213A-R20-C06,SA1090-A96213A-R20-C08 \
        #     -O cellsnptest \
        #     -R haps.vcf \
        #     -p 1 \
        #     --cellTAG None \
        #     --UMItag None \
        #     --gzip \
        #     --minCOUNT 0