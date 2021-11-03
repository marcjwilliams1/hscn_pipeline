rule schnapps:
    input:
        haplotypes = "results/{sample}/counthaps/allele_counts_perblock.csv.gz",
        qc = "results/{sample}/hmmcopy_results/metrics.csv.gz",
        hmmcopy = "results/{sample}/hmmcopy_results/reads.csv.gz"
    output:
        qc = "results/{sample}/schnapps/qc.csv.gz",
        qcplot = report("results/{sample}/schnapps/qc.png", category = "QC plot"),
        csv = "results/{sample}/schnapps/schnapps.csv.gz",
        heatmap = report("results/{sample}/schnapps/heatmap.png", category = "Heatmaps"),
        heatmapraw = report("results/{sample}/schnapps/heatmapraw.png", category = "Heatmaps Raw"),
        rdata = "results/{sample}/schnapps/schnapps.Rdata"
    threads: 10
    resources: mem_mb=1024 * 4
    params:
        mincells=config["schnapps"]["mincells"],
        qualfilter=config["schnapps"]["qualfilter"]
    singularity: "docker://marcjwilliams1/schnapps:latest"
    shell:
        """
        Rscript scripts/run_schnapps.R  \
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