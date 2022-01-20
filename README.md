# Haplotype specific copy number in single cell whole genome sequencing

This is a snakemake pipeline to produce estimates of haplotype specific copy number in single cells using [signals](https://github.com/shahcompbio/signals). signals takes as input total copy number estimates and per cell haplotype block counts as outputted from the Shah lab's [single cell pipeline](https://github.com/shahcompbio/single_cell_pipeline). If possible, it is advised to use the official realease and the steps outlined [here](https://github.com/shahcompbio/single_cell_pipeline/blob/master/docs/source/install.md). The official release will likely be more up to date and continue to have bug fixes. This pipeline is provided as a lightweight alternative to go along with the our paper (link), and for those with other technologies wanting to explore the code.

## Input files

The required input files are single cell tumour bam files and a bulk normal bam file. These should be listed in a csv file along with sample ID's. A normal bam file should also be included here for the haplotype calling, this should have the value `T` in the `normal` column. An example csv is provided in the `metadata` folder. You can have multiple patients per file, just ensure that they have different `sample` names. The location and name of this metadata file should be accurate in the top of the config.yaml file. 

If you have data produced from the 10X CNV assay you can demultiplex the bam into single cell bams. See the `demultiplex10X` folder for a seperate snakemake pipeline to do that.

## Reference files

A number of reference files are required. These can be pulled as follows:

```
wget https://singlecelltestsets.s3.amazonaws.com/refdata_full_genome.tar.gz
tar -xvf refdata_full_genome.tar.gz
```

Make sure this folder is in the same folder as the Snakefile and is called `refdata` (or change the paths in the config.yaml).

## Config options

The default bin size is 0.5Mb, you can change this by modifying the `config.yaml` file. You'll need to change the `bin_size` parameter and the reference wig files. You can find these files at different bin sizes [here](https://zenodo.org/record/5549581).
```
wget https://zenodo.org/record/5549581/files/hmmcopy_wigfiles.tar.gz
tar -xvf hmmcopy_wigfiles.tar.gz
mv forzip/ hmmcopy_wigfiles
```

Reference files are provided for hg19. If you want to use hg38, you'll need to supply a *fa file, change the wig files (see the zenodo link for hg38 files) and change the reference files used by shapeit.

## What does the pipeline do?

The pipeline performs the following steps:

1. Counts the number of reads in each bin in each cell
2. Performs GC correction on the read counts
3. Calls total copy number using HMMcopy based on GC corrected read counts
4. Calculates quality scores and metrics for all the copy number profiles
5. Infer haplotype blocks from a matched normal genome using [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html).
6. Genotypes SNPs in the single cells using [cellsnp-lite](https://cellsnp-lite.readthedocs.io/en/latest/)
7. Merges the per cell SNP counts into haplotype block counts
8. Estimates haplotype specific copy number using signals.

## Things the pipeline does not do...

Please note that currently this pipeline does not classify cell cycle phase, attempt to detect contamination from bacteria or other organisms or call somatic SNV's and breakpoints. See the officially supported pipelines for all these additional features. This is only provided as a lightweight setup to run the necessary steps to produce the inputs for signals.

## Outputs

The final output will appear in `results/{sample}/signals/`. This includes an Rdata file of the  signals object which can be opened in R and then explored interactively using  signals. Also included is a csv file of allele specific calls per cell, qc plots and heatmaps of all the data.

## Running recommendations

Depending on the number of cells you have, the pipeline can take quite a long time so running on a HPC is recommended. Using snakemake profiles is an easy way to achieve this. Once you have set up a profile you can submit a job to your cluster and snakemake will take care of the rest and submit all the steps as seperate jobs. For example the following code snippet will work for a cluster running LSF.

```
#!/bin/bash
#BSUB -J  signals
#BSUB -n 1
#BSUB -R rusage[mem=4]
#BSUB -W 240:00
#BSUB -eo logs/cluster/%J.stderr

source ~/miniconda3/etc/profile.d/conda.sh
conda activate snakemake
module load singularity/3.6.2

snakemake \
  --profile lsf \
  --keep-going \
  --restart-times 0 \
  --singularity-args "--bind /directories/to/bind" \
  --rerun-incomplete

snakemake --report report.html
```

Ideally you'll have singularity installed, if so snakemake will pull all the required software. If this is not possible you can try the conda environments in `envs/`. You'll need to modify the snakemake rules to use conda rather than singularity.
