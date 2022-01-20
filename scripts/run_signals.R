

read_copynumber_dlp <- function(cnpaths,
                                metricspaths,
                                sample_ids = NULL,
                                patient = NULL,
                                filtercells = TRUE,
                                cols_to_keep = c("cell_id", "chr", "start", "end", "map", "copy", "state"),
                                mappability = 0.99,
                                s_phase_filter = FALSE,
                                quality_filter = 0.75,
                                filter_reads = 0.1e6){
  
  metricslist <- list()
  i <- 1
  for (x in metricspaths){
    if (!is.null(sample_ids)) {
      metricslist[[x]] <- data.table::fread(x)[, sample_id := sample_ids[i]]
      i <- i + 1
    } else {
      metricslist[[x]] <- data.table::fread(x)
      metricslist[[x]]$file <- x
    }
  }
  
  #merge metrics files
  metricsdata <- data.table::rbindlist(metricslist, use.names = T, fill = T)
  
  if (filtercells == TRUE){
    
    message(paste0("Total number of cells pre-filtering: ", dim(metricsdata)[1]))
    
    #find cells that pass filtering criteria
    cells_to_keep <- metricsdata %>%
      .[quality > quality_filter] #%>%
      #.[is_contaminated == FALSE] %>%
      #.[total_mapped_reads > filter_reads] %>%
      #.[!stringr::str_detect(experimental_condition, "NTC|NCC|gDNA|GM|CONTROL|control")]
    
    if (s_phase_filter){
      message("Removing s-phase cells")
      cells_to_keep <- cells_to_keep[is_s_phase == FALSE]
    } else{
      message("Total number of s-phase cells: ", sum(cells_to_keep$is_s_phase))
    }
    
    metricsdata <- metricsdata[cell_id %in% cells_to_keep$cell_id]
    message(paste0("Total number of cells post-filtering: ", dim(metricsdata)[1]))
  }
  
  # read in copy number data, only keeping filtered cells
  cnlist <- list()
  for (x in cnpaths){
    cn_ <- data.table::fread(x,
                             colClasses = list(character = c("chr", "cell_id"),
                                               integer = c("start", "end", "reads", "state"),
                                               numeric = c("map", "gc", "copy")))
    cn_ <- cn_[cell_id %in% metricsdata$cell_id]
    cnlist[[x]] <- cn_
  }
  
  #bind copy number data
  cndata <- data.table::rbindlist(cnlist, use.names = T, fill = T)
  cndata <- cndata[, ..cols_to_keep]
  cndata <- cndata[map > mappability]
  
  rm(cnlist)
  
  return(list(cn = cndata %>% as.data.frame(), metrics = metricsdata %>% as.data.frame()))
}

read_copynumber_dlpqc <- function(metricspaths,
                                sample_ids = NULL,
                                patient = NULL,
                                filtercells = TRUE,
                                mappability = 0.99,
                                s_phase_filter = FALSE,
                                quality_filter = 0.75,
                                filtercond = "NTC|NCC|gDNA|GM|CONTROL|control",
                                filter_reads = 0.1e6){
  
  metricslist <- list()
  i <- 1
  for (x in metricspaths){
    if (!is.null(sample_ids)) {
      #read in metrics data
      metricslist[[x]] <- data.table::fread(x)[, sample_id := sample_ids[i]]
      i <- i + 1
    } else {
      metricslist[[x]] <- data.table::fread(x)
      metricslist[[x]]$file <- x
    }
  }
  
  metricsdata <- data.table::rbindlist(metricslist, use.names = T, fill = T)
  
  if (filtercells == TRUE){
    
    message(paste0("Total number of cells pre-filtering: ", dim(metricsdata)[1]))
    
    #find cells that pass filtering criteria
    cells_to_keep <- metricsdata %>%
      .[quality > quality_filter] #%>%
      #.[is_contaminated == FALSE] %>%
      #.[is_s_phase == FALSE] %>%
      #.[total_mapped_reads > filter_reads] %>%
      #.[!stringr::str_detect(experimental_condition, filtercond)]
    
    if (s_phase_filter){
      cells_to_keep <- cells_to_keep[is_s_phase == FALSE]
    }
    
    metricsdata <- metricsdata[cell_id %in% cells_to_keep$cell_id]
    message(paste0("Total number of cells post-filtering: ", dim(metricsdata)[1]))
  }

  
  return(metrics = metricsdata %>% as.data.frame())
}

callhscn <- function(cn,
                     qc,
                     haps,
                     diploidcutoff = 0.05,
                     mincells = 8,
                     frachaps = 0.8,
                     ncores = 1){

  cn <- as.data.table(cn)
  qc <- as.data.table(qc)
  haps <- as.data.table(haps)

  message(paste0("Number of cells (pre filtering for normal cells): ", length(unique(cn$cell_id))))
  #find the fraction of each cell that is diploid and fraction that == cell ploidy
  cn_ploidy <- cn[, list(state_mean = mean(state), 
                         state_var = var(state), 
                         frac_nondiploid = sum(state != 2) / .N, 
                         pct_cell_ploidy = sum(state == schnapps:::Mode(state)) / .N, 
                         mode_state = schnapps:::Mode(state)), by = "cell_id"] %>% 
    .[order(pct_cell_ploidy, decreasing = TRUE)]
  
  qc <- qc %>% 
    left_join(cn_ploidy %>% dplyr::select(cell_id, frac_nondiploid, pct_cell_ploidy, state_mean), by = "cell_id")
  
  #find cells that are probably diploid or that are diploid + misscalled ploidy 
  cn_ploidy <- cn_ploidy[pct_cell_ploidy < (1 - diploidcutoff)]
  cutoff <- mean(cn_ploidy$pct_cell_ploidy) + 3 * sd(cn_ploidy$pct_cell_ploidy) #compute outlier cutoff
  message(paste0("Fraction ploidy cutoff: ", round(cutoff, 3)))

  if (diploidcutoff == 0.0){
    message("Removing misscalled diploid cells")
    cells_to_remove <- cn_ploidy[pct_cell_ploidy > 0.95 & mode_state > 2]
    non_diploidcells <- cn_ploidy[!(cell_id %in% cells_to_remove$cell_id)]
  } else{
    message("Finding tumour cells")
    non_diploidcells <- cn_ploidy[frac_nondiploid > diploidcutoff & pct_cell_ploidy < cutoff]
  }
  
  #filter out diploid cells
  cn <- cn[cell_id %in% non_diploidcells$cell_id]
  qc <- qc %>% 
    dplyr::mutate(pct_cell_ploidy_cutoff = cutoff, is_tumour_cell = cell_id %in% non_diploidcells$cell_id)
  haps <- haps[cell_id %in% unique(cn$cell_id)]
  
  message(paste0("Total number of cells after filtering diploid or misscalled ploidy: ", length(unique(cn$cell_id))))
  
  if (length(unique(cn$cell_id)) < 10){
    stop("Number of cells is < 10, unable to run schnapps")
    hscn <- list(data = NULL, message = "Number of cells is less than 10, unable to run schnapps")
    return(hscn)
  }

  if (length(unique(cn$cell_id)) == 0){
    stop("No cells returning NULL")
    return(NULL)
  }

  message("Infer HSCN")
  hscn <- callHaplotypeSpecificCN(cn,
                                  haps,
                                  likelihood = "auto",
                                  mincells = mincells,
                                  cluster_per_chr = TRUE, 
                                  ncores = ncores)
  
  #add library and sample IDs to qc
  qc$sample_id <- unlist(schnapps:::get_library_labels(qc$cell_id, idx=1))
  qc$library_id <- unlist(schnapps:::get_library_labels(qc$cell_id, idx=2))

  #add qc to schnapps object
  hscn$qc_per_cell <- dplyr::left_join(hscn$qc_per_cell, qc, by = "cell_id")
  
  return(hscn)
}

callascn <- function(hscn,
                     qc,
                     haps,
                     diploidcutoff = 0.05,
                     ncores = 1){
  
  message("Infer ASCN")
  ascn <- callAlleleSpecificCNfromHSCN(hscn, ncores = ncores)
  ascn <- filtercn(ascn)
  ascn$qc_per_cell <- dplyr::left_join(ascn$qc_per_cell, hscn$qc_per_cell)
  
  return(ascn)
}

library(argparse)
library(tidyverse)
library(data.table)
library(schnapps)
library(ggplot2)
library(cowplot)

message(paste0("schnapps version: ", packageVersion("schnapps")))

parser <- ArgumentParser()

parser$add_argument("--hmmcopyreads", default="character", nargs = "+",
                    help = "hmmcopy reads files")
parser$add_argument("--hmmcopyqc", default=NULL, type="character", nargs = "+",
                    help="hmmcopy QC files")
parser$add_argument("--allelecounts", default=NULL, type="character", nargs = "+",
                    help="hmmcopy QC files")
parser$add_argument("--ncores", default=NULL, type="integer",
                    help="Number of cores in schnapps inference")
parser$add_argument("--csvfile", default=NULL, type="character",
                    help="output csvfile")
parser$add_argument("--qccsvfile", default=NULL, type="character",
                    help="output csvfile")
parser$add_argument("--Rdatafile", default=NULL, type="character",
                    help="output Rdata file")
parser$add_argument("--heatmap", default=NULL, type="character",
                    help="heatmap plot")
parser$add_argument("--heatmapraw", default=NULL, type="character",
                    help="heatmap plot of raw data")
parser$add_argument("--qcplot", default=NULL, type="character",
                    help="QC plot")
parser$add_argument("--diploidcutoff", default=0.05, type="double",
                    help="Cutoff for filtering diploid cells, this is the fraction of the genome that is non diploid")
parser$add_argument("--maxcellsplotting", default=2500, type="integer",
                    help="Max number of cells to plot in the heatmap")
parser$add_argument("--mincells", default=8, type="integer",
                    help="Max number of cells to plot in the heatmap")
parser$add_argument("--qualfilter", default=0.75, type="double",
                    help="Quality score cutoff")
parser$add_argument("--sphasefilter", action = "store_false",
                    help="Filter out sphase cells or not")

args <- parser$parse_args()
print(args)

cndata <- read_copynumber_dlp(cnpaths = args$hmmcopyreads,
                              metricspaths = args$hmmcopyqc, 
                              s_phase_filter = args$sphasefilter,
                              quality_filter = args$qualfilter)
cell_ids <- unique(cndata$metrics$cell_id)
message(paste0("Number of cells: ", length(cell_ids)))

bin_ids <- as.data.table(cndata$cn)[, c("chr", "start", "end")] %>%
  unique(., by = c("chr", "start", "end")) %>%
  .[, binid := paste(chr, start, end, sep = "_")] %>%
  .$binid

message("Read in haplotypes data")
haplotypes <- fread(args$allelecounts, colClasses = list(character = c("chr")))

message("Call haplotype specific copy number")
res <- callhscn(cndata$cn,
                cndata$metrics,
                haplotypes,
                frachaps = 0.8,
                mincells = args$mincells,
                diploidcutoff = args$diploidcutoff,
                ncores = args$ncores)
print(res)

message(paste0("Number of bins: ", dim(res$data)[1]))

message("Call allele specific copy number")
res2 <- callascn(res, ncores = args$ncores)
res2$likelihood <- res$likelihood
print(res2)

message("Make QC plots")
g1 <- plotBAFperstate(res$data %>% filter(totalcounts > 9))+
  theme(
    panel.background = element_rect(fill = "white"), # bg of the panel
    plot.background = element_rect(fill = "white"), # bg of the plot
    legend.background = element_rect(fill = "white"), # get rid of legend bg
    legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
  )

g2 <- plot_variance_state(res$data %>% filter(totalcounts > 9))+
  theme(
    panel.background = element_rect(fill = "white"), # bg of the panel
    plot.background = element_rect(fill = "white"), # bg of the plot
    legend.background = element_rect(fill = "white"), # get rid of legend bg
    legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
  )

g3 <- res$qc_per_cell %>% 
  mutate(Ploidy = paste0(ploidy)) %>% 
  ggplot(aes(x = average_distance, fill = Ploidy)) + 
  geom_histogram(alpha = 0.7, bins = 50) +
  xlab("Average distance raw BAF\nto expected BAF per cell") +
  ylab("Counts") +
  theme(legend.position = "none") +
  theme(
    panel.background = element_rect(fill = "white"), # bg of the panel
    plot.background = element_rect(fill = "white"), # bg of the plot
    legend.background = element_rect(fill = "white"), # get rid of legend bg
    legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
  )

g4 <- res$qc_per_cell %>% 
  mutate(Ploidy = paste0(ploidy)) %>% 
  ggplot(aes(x = totalhapcounts / 1e6, fill = Ploidy)) + 
  geom_histogram(alpha = 0.7, bins = 50) +
  xlab("Total number of hap counts per cell (Millions)") +
  ylab("Counts") +
  theme(legend.position = "none") +
  theme(
    panel.background = element_rect(fill = "white"), # bg of the panel
    plot.background = element_rect(fill = "white"), # bg of the plot
    legend.background = element_rect(fill = "white"), # get rid of legend bg
    legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
  )

g5 <- res$qc_per_cell %>% 
  mutate(Ploidy = paste0(ploidy)) %>% 
  ggplot(aes(x = pct_cell_ploidy, fill = Ploidy)) + 
  geom_histogram( alpha = 0.7, bins = 50) +
  xlab("Fraction of genome = ploidy per cell") +
  ylab("Counts") +
  theme(legend.position = c(0.1, 0.9)) +
  theme(
    panel.background = element_rect(fill = "white"), # bg of the panel
    plot.background = element_rect(fill = "white"), # bg of the plot
    legend.background = element_rect(fill = "white"), # get rid of legend bg
    legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
  )

g6 <- plot_clusters_used_for_phasing(res) +
  theme(
    panel.background = element_rect(fill = "white"), # bg of the panel
    plot.background = element_rect(fill = "white"), # bg of the plot
    legend.background = element_rect(fill = "white"), # get rid of legend bg
    legend.box.background = element_rect(fill = "white") # get rid of legend panel bg
  )

gall <- cowplot::plot_grid(cowplot::plot_grid(g1, g2, align = "h", axis = "tb"), cowplot::plot_grid(g3, g4, g5, ncol = 3), g6, ncol = 1, rel_heights = c(0.25, 0.25, 0.5))

cowplot::save_plot(args$qcplot,
                   gall,
                   base_width = 15, base_height = 20)

#filter out cells where BAF and inferred state indicate possible errors
cell_ids <- unique(res$data$cell_id)
message(paste0("Number of cells (pre final QC filtering): ", length(cell_ids)))
res <- filtercn(res)
cell_ids <- unique(res$data$cell_id)
message(paste0("Number of cells (after final QC filtering): ", length(cell_ids)))

message("Cluster data")
cl <- umap_clustering(res$data, 
                      field = "copy", 
                      hscn = TRUE,
                      minPts = round(0.03 * length(cell_ids)))

message("Write csv file")
fwrite(x = res$data, file = args$csvfile)
fwrite(x = res$qc_per_cell, file = args$qccsvfile)

message("Write Rdata file")
saveRDS(file = args$Rdatafile, object = list(hscn = res, ascn = res2, cl = cl))

message("Make heatmap")

if (dim(cl$clustering)[1] > args$maxcellsplotting){
  message(paste0(dim(cl$clustering)[1], " total cells, downsampling to ", args$maxcellsplotting, " for plotting"))
  cells <- dplyr::sample_n(cl$clustering, args$maxcellsplotting) %>% pull(cell_id)
} else{
  cells <- cl$clustering %>% pull(cell_id)
}

if (length(unique(cl$clustering$clone_id)) > 20){
  show_clone_label <- FALSE
} else{
  show_clone_label <- TRUE
}



h1 <- plotHeatmap(res$data %>% dplyr::filter(cell_id %in% cells),
                  tree = ape::keep.tip(cl$tree, cells),
                  reorderclusters = T,
                  show_clone_label = show_clone_label,
                  clusters = cl$clustering %>% dplyr::filter(cell_id %in% cells),
                  plottree = FALSE,
                  plotcol = "state",
                  str_to_remove = "SPECTRUM-OV-")
h2 <- plotHeatmap(res$data %>% dplyr::filter(cell_id %in% cells),
                  tree = ape::keep.tip(cl$tree, cells),
                  reorderclusters = T,
                  clusters = cl$clustering %>% dplyr::filter(cell_id %in% cells),
                  plottree = FALSE,
                  plotcol = "state_phase", 
                  show_clone_label = F, 
                  show_library_label = F,
                  str_to_remove = "SPECTRUM-OV-")

h3 <- plotHeatmap(res$data %>% dplyr::filter(cell_id %in% cells),
                  tree = ape::keep.tip(cl$tree, cells),
                  reorderclusters = T,
                  show_clone_label = show_clone_label,
                  clusters = cl$clustering %>% dplyr::filter(cell_id %in% cells),
                  plottree = FALSE,
                  plotcol = "copy",
                  str_to_remove = "SPECTRUM-OV-")

h4 <- plotHeatmap(res$data %>% dplyr::filter(cell_id %in% cells),
                  tree = ape::keep.tip(cl$tree, cells),
                  reorderclusters = T,
                  clusters = cl$clustering %>% dplyr::filter(cell_id %in% cells),
                  plottree = FALSE,
                  plotcol = "BAF",
                  show_clone_label = F,
                  show_library_label = F,
                  str_to_remove = "SPECTRUM-OV-")


png(args$heatmap, width = 25, height = 12, units = "in", res = 300)
print(ComplexHeatmap::draw(h1 + h2,
                            ht_gap = unit(0.6, "cm"),
                            column_title = paste0("Total number of cells: ", length(unique(res$data$cell_id)), "\nNumber of cells plotted: ", length(cells)),
                            column_title_gp = grid::gpar(fontsize = 20),
                            heatmap_legend_side = "bottom",
                            annotation_legend_side = "bottom",
                            show_heatmap_legend = TRUE))
dev.off()

png(args$heatmapraw, width = 25, height = 12, units = "in", res = 300)
print(ComplexHeatmap::draw(h3 + h4,
                            ht_gap = unit(0.6, "cm"),
                            column_title = paste0("Total number of cells: ", length(unique(res$data$cell_id)), "\nNumber of cells plotted: ", length(cells)),
                            column_title_gp = grid::gpar(fontsize = 20),
                            heatmap_legend_side = "bottom",
                            annotation_legend_side = "bottom",
                            show_heatmap_legend = TRUE))
dev.off()

message("Finished")

print(sessionInfo())
