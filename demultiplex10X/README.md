## Demultiplex 10X CNV bam

This snakemake pipeline will demultiplex a 10X bam into single cell bams. To run you'll need to change some of the files in the Snakefile. First, set the metrics file (top of the Snakefile) to the *cell_metrics.csv file outputted by the 10X pipeline, secondly change the location of the bam file in the samtools_view rule.