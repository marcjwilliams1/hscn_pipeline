library(data.table)
library(tidyverse)

save.image("mytest.rds")

getdf <- function(x, ...) {
  
  s <- Matrix::summary(x)
  
  row <- s$i
  if (!is.null(rownames(x))) {
    row <- rownames(x)[row]
  }
  col <- s$j
  if (!is.null(colnames(x))) {
    col <- colnames(x)[col]
  }
  
  ret <- data.table(
    snpid = row, 
    cellid = col, 
    value = s$x,
    stringsAsFactors = FALSE
  )
  ret
}

message("Read in haplotypes file")
haps <- fread(snakemake@input$haplotypes) %>%
  pivot_wider(names_from = "allele", values_from = "allele_id", names_prefix = "allele") %>%
  rename(ref_phase = allele0, alt_phase = allele1, CHROM = chromosome, POS = position) %>% 
  select(-ref, -alt) %>% 
  as.data.table()

message("Read in list of barcodes")
samples <- fread(snakemake@input$results_sample, header = FALSE, col.names = "barcode") %>%
  .[, cellid := 1:.N]

message("Read in VCF")
vcf <- fread(snakemake@input$vcf)

message("Format VCF files")
vcf <- vcf %>%
  .[, snpid := 1:.N] %>%
  rename(CHROM = `#CHROM`) %>% 
  select(-ID, -QUAL, -FILTER, -INFO)

message("Read in matrix files")
DP <- Matrix::readMM(snakemake@input$results_DP_mtx)
DPdf <- getdf(DP) %>%
  merge(samples, on = "cell_id", all.x = TRUE) %>%
  rename(DP = value)
AD <- Matrix::readMM(snakemake@input$results_AD_mtx)
ADdf <- getdf(AD) %>%
  merge(samples, on = "cell_id", all.x = TRUE) %>%
  rename(AD = value)

message("Matrix counts and vcf")
newdf <- ADdf[DPdf, nomatch = NA, on = c("cellid", "snpid", "barcode")] %>%
  setnafill(., type = "const", cols = c("AD"), fill = 0) %>% 
  merge(., vcf, all.x = TRUE, on = "snpid") 

message("Merge haplotypes and phase alleles")
newdf <- haps[newdf, on = c("CHROM", "POS")] %>% 
  .[, allele0 := fifelse(ref_phase == 0, DP - AD, AD)] %>% 
  .[, allele1 := fifelse(ref_phase == 0, AD, DP - AD)] %>% 
  dplyr::select(-ref_phase, -alt_phase, -AD, -snpid, -cellid) %>%
  dplyr::rename(cell_id = barcode, position = POS, ref = REF, alt = ALT, chr = CHROM) %>%
  dplyr::select(cell_id, chr, position, ref, alt, hap_label, allele0, allele1) %>%
  .[, totalcounts := allele0 + allele1] %>% 
  .[order(cell_id, chr, position)]

fwrite(x = newdf, file = snakemake@output$alldata)

binsize <- as.numeric(snakemake@config$hmmcopy$bin_size)

summed_data <- newdf %>% 
  .[, binid := floor(position / binsize) * binsize + 1] %>%
  .[, list(allele0 = sum(allele0),
           allele1 = sum(allele1),
           start = min(position),
           end = max(position)), by = c("cell_id", "chr", "hap_label", "binid")] %>% 
  .[, start := floor(start / binsize) * binsize + 1] %>%
  .[, end := start + binsize - 1] %>%
  .[, binid := NULL] %>% 
  .[, totalcounts := allele0 + allele1] %>% 
  .[order(cell_id, chr, hap_label)]
fwrite(summed_data, file = snakemake@output$perblock)
