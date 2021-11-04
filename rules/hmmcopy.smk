def getbamfile(wildcards):
    x = metadata[metadata["cell_id"] == wildcards.cell_id]
    x = x[x["sample"] == wildcards.sample]
    return x["bamfiles"]

rule read_counter:
    input: getbamfile
    output: 
        mydir = directory("results/{sample}/read_counts/{cell_id}/"),
        file2 = "results/{sample}/read_counts/{cell_id}/wig.txt",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 30
    params:
        exclude_list=config["hmmcopy"]["exclude_list"],
        bin_size=config["hmmcopy"]["bin_size"],
        mapq=config["hmmcopy"]["min_mqual"]
    singularity: "docker://quay.io/singlecellpipeline/single_cell_pipeline_hmmcopy:v0.8.9"
    shell:
        """
        mkdir -p {output[0]}
        python scripts/read_counter.py {input[0]} {output[1]} \
            -w {params.bin_size} \
            --exclude_list {params.exclude_list} \
            --mapping_quality_threshold {params.mapq}
        """

rule correct_read_counts:
    input: "results/{sample}/read_counts/{cell_id}/wig.txt"
    output: "results/{sample}/corrected_read_counts/{cell_id}.csv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 30
    singularity: "docker://quay.io/singlecellpipeline/single_cell_pipeline_hmmcopy:v0.8.9"
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
    singularity: "docker://quay.io/singlecellpipeline/single_cell_pipeline_hmmcopy:v0.8.9"
    shell:
        """
        mkdir -p {output[0]}
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
    singularity: "docker://quay.io/singlecellpipeline/single_cell_pipeline_hmmcopy:v0.8.9"
    script: "../scripts/extend_qc.py"


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
    resources: mem_mb=lambda wildcards, attempt: attempt * 1024 * 50
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
    resources: mem_mb=lambda wildcards, attempt: attempt * 1024 * 50
    run:
        dist_list = []
        for f in input:
            print(f)
            dist_temp = pd.read_csv(f)
            dist_list.append(dist_temp)
        dist = pd.concat(dist_list, ignore_index=True)
        dist.to_csv(output.merged, sep=',')