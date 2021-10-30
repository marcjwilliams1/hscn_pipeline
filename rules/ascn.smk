rule schnapps:
    input:
        haplotypes = "results/{sample}/counthaps/allele_counts_perblock.csv.gz",
        qc = "results/{sample}/hmmcopy_results/metrics.csv.gz",
        hmmcopy = "results/{sample}/hmmcopy_results/reads.csv.gz"
    output:
        qc = "results/{sample}/schnapps/qc.csv.gz",
        qcplot = "results/{sample}/schnapps/qc.png",
        csv = "results/{sample}/schnapps/schnapps.csv.gz",
        heatmap = "results/{sample}/schnapps/heatmap.png",
        heatmapraw = "results/{sample}/schnapps/heatmapraw.png",
        rdata = "results/{sample}/schnapps/schnapps.Rdata"
    threads: 10
    resources: mem_mb=1024 * 4
    singularity: "docker://marcjwilliams1/schnapps"
    shell:
        """
        pwd
        #module load R/R-3.6.1
        Rscript ../scripts/run_schnapps.R  \
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
            --sphasefilter
        """