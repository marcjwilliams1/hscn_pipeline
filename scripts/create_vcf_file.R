library(data.table)
library(tidyverse)

haps <- fread(snakemake@input$haplotypes)
#genomes_1000 <- fread(snakemake@config$thousand_genomes_snps)
#names(genomes_1000) <- c("chromosome", "position", "ref", "alt")

#haps <- merge(haps, genomes_1000, on = c("chromosome", "position"), all.x = T)

haps <- select(haps, chromosome, position, ref, alt) %>% 
  filter(alt != "-")

#rename columns and add columns required by VCF
haps$`#CHROM` <- haps$chromosome
haps$POS <- haps$position
haps$REF <- haps$ref
haps$ALT <- haps$alt
haps$QUAL <- "."
haps$ID <- "."
haps$FILTER <- "PASS"
haps$INFO <- "."

haps <- select(haps, `#CHROM`, POS, ID, REF, ALT,QUAL,FILTER,INFO) %>% 
  distinct()

write("##fileformat=VCFv4.2", snakemake@output$haplotypesvcf)
fwrite(haps, file = snakemake@output$haplotypesvcf, sep = "\t", append = T, col.names = T)

# fread("../../../schnapps_input_pipeline/metadata/samples.csv") %>% 
#   .[, cell_id := str_remove(cell_id, "/work/shah/users/william1/projects//schnapps_input_pipeline/testdata//")] %>% 
#   .[, bamfiles := str_replace(bamfiles, "\\//", "/")] %>% 
#   fwrite("../../../schnapps_input_pipeline/metadata/samples.csv")
