rule signals:
    input:
        haplotypes = "results/{sample}/counthaps/allele_counts_perblock.csv.gz",
        qc = "results/{sample}/hmmcopy_results/metrics.csv.gz",
        hmmcopy = "results/{sample}/hmmcopy_results/reads.csv.gz"
    output:
        qc = "results/{sample}/signals/qc.csv.gz",
        qcplot = report("results/{sample}/signals/qc.png", category = "QC plot"),
        csv = "results/{sample}/signals/signals.csv.gz",
        heatmap = report("results/{sample}/signals/heatmap.png", category = "Heatmaps"),
        heatmapraw = report("results/{sample}/signals/heatmapraw.png", category = "Heatmaps Raw"),
        rdata = "results/{sample}/signals/signals.Rdata"
    threads: 10
    resources: mem_mb=lambda wildcards, attempt: attempt * 1024 * 5
    params:
        mincells=config["signals"]["mincells"],
        qualfilter=config["signals"]["qualfilter"]
    singularity: "docker://marcjwilliams1/signals:latest"
    shell:
        """
        Rscript scripts/run_signals.R  \
            --hmmcopyqc {input.qc} \
            --hmmcopyreads {input.hmmcopy}  \
            --allelecounts {input.haplotypes}  \
            --ncores {threads} \
            --qcplot {output.qcplot} \
            --heatmap {output.heatmap} \
            --heatmapraw {output.heatmapraw} \
            --csvfile {output.csv} \
            --qccsvfile {output.qc} \
            --Rdatafile {output.rdata} \
            --sphasefilter \
            --qualfilter {params.qualfilter} \
            --mincells {params.mincells}
        """
