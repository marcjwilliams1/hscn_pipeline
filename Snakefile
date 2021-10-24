import pandas as pd

configfile: "config.yaml"
CELL_IDS, = glob_wildcards(config["bamdir"] + "{id}.bam")
metadata = pd.read_csv(config["samplemetadata"])
samples = metadata["samples"].to_list().unique()


rule all:
    input:
        "results/hmmcopy_results/metrics.csv.gz",
        "results/hmmcopy_results/reads.csv.gz",
        "results/inferhaps/haplotypes.csv.gz"

def getbamfile(wildcards):
    x = samples[samples["cell_id"] == wildcards.cell_id]
    x = x[x["sample"] == wildcards.sample]
    return x["bamfiles"][0]

rule read_counter:
    input: "testdata/{cell_id}.bam"
    output: 
        mydir = directory("results/read_counts/{cell_id}/"),
        file1 = "results/read_counts/{cell_id}/test.txt",
        file2 = "results/read_counts/{cell_id}/wig.txt",
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
    input: "results/read_counts/{cell_id}/wig.txt"
    output: "results/corrected_read_counts/{cell_id}.csv"
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
    input: "results/corrected_read_counts/{cell_id}.csv"
    output: 
        celldir = temp(directory("results/hmmcopy_results_temp/{cell_id}/")),
        reads = temp("results/hmmcopy_results_temp/{cell_id}/0/reads.csv"),
        metrics = temp("results/hmmcopy_results_temp/{cell_id}/0/metrics.csv"),
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
        reads = "results/hmmcopy_results_temp/{cell_id}/0/reads.csv",
        celldir = "results/hmmcopy_results_temp/{cell_id}/",
    output:
        reads = "results/hmmcopy_results/percell/{cell_id}_reads.csv",
    shell:
        """
        mv {input.reads} {output.reads}
        """

rule add_metrics:
    input: 
        metrics = "results/hmmcopy_results_temp/{cell_id}/0/metrics.csv",
        insrt = "results/stats/{cell_id}.isize.txt",
        flgstat = "results/stats/{cell_id}.bam.flagstat",
        celldir = "results/hmmcopy_results_temp/{cell_id}/"
    output:
        metrics = "results/hmmcopy_results/percell/{cell_id}_metrics.csv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 30
    conda: "envs/scp.yml"
    script: "scripts/extend_qc.py"


rule insert_size:
    input:
        "testdata/{cell_id}.bam"
    output:
        txt="results/stats/{cell_id}.isize.txt",
        pdf="results/stats/{cell_id}.isize.pdf"
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
        bam="testdata/{cell_id}.bam"
    output:
        "results/stats/{cell_id}.summary.txt"
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
        "testdata/{cell_id}.bam"
    output:
        "results/stats/{cell_id}.bam.flagstat"
    wrapper:
        "0.78.0/bio/samtools/flagstat"


rule mergereads:
    input: 
        expand("results/hmmcopy_results/percell/{cell_id}_reads.csv", cell_id = CELL_IDS)
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

rule mergemetrics:
    input: 
        expand("results/hmmcopy_results/percell/{cell_id}_metrics.csv", cell_id = CELL_IDS)
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

rule create_yaml_normal_inferhaps:
    input: "testdata/OV2295_normal.bam"
    output: temp("results/infer_haps.yaml")
    shell:
        """
        printf "normal:\n  bam: {input[0]}" > {output[0]}
        """

rule inferhaps:
    input:
        yaml = "results/infer_haps.yaml",
    output:
        inferhapsdir = directory("results/inferhaps/"),
        haplotypes = "results/inferhaps/haplotypes.csv.gz",
        pipelinedir = temp(directory("results/inferhaps/pipeline/")),
    params:
        config = "config.yaml",
        singularity = config["infer_haps_singularity"]
    threads: 40
    resources:
        mem_mb=1024 * 4
    shell:
        """
        echo {output}
        echo {input}
        echo {params}
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

#cellsnp-lite -s testdata/SA1090-A96213A-R20-C06.bam,testdata/SA1090-A96213A-R20-C08.bam -I SA1090-A96213A-R20-C06,SA1090-A96213A-R20-C08 -O cellsnptest -R haps.vcf -p 1 --cellTAG None --UMItag None --gzip --minCOUNT 0

