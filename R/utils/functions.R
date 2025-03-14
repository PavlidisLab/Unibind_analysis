## Functions used throughout analysis
## -----------------------------------------------------------------------------

library(GenomicRanges)

## TODO: range table functions -> package

# source("/home/amorin/Projects/TR_aggregation/R/utils/range_table_functions.R")


# Read unibind experiment (peak table) as a data frame

load_unibind <- function(path) {
  
  peak <- read.delim(
    file = path,
    header = FALSE,
    col.names = c("Chromosome", "Start", "End", "Motif", "Score", "Strand"),
    sep = "\t"
  )
  
  return(peak)
}




# remove 'chr' prefix of chromosome identifiers, coerce mitochondrial to 
# 'MT' and only keep standard autosomal and sex chromosomes

get_standard_chr <- function(range_table) {
  
  stopifnot("Chromosome" %in% names(range_table))
  
  range_table$Chromosome <- str_replace(range_table$Chromosome, "^chr", "")
  range_table$Chromosome <- str_replace(range_table$Chromosome, "^M$", "MT")
  range_table <- filter(range_table, Chromosome %in% c(1:22, "MT", "X", "Y"))
  return(range_table)
}




# Convert a peak table to a GR object. If format is TRUE, the table will
# first be processed to be consistent with standard used for overlapping with
# a gene annotation table. 

peak_to_gr <- function(peak_table) {
  
  gr <- makeGRangesFromDataFrame(peak_table,
                                 keep.extra.columns = TRUE, 
                                 ignore.strand = TRUE)
  return(gr)
}



# Convert a unibind peak df into a 1bp Genomic Ranges object and filter
# black listed regions

unibind_to_gr <- function(peak_table, bl_gr) {
  
  peak_gr <- peak_table %>%
    get_standard_chr() %>%
    peak_to_gr() %>%
    filter_blacklist(bl_gr)
  
  # Hacky means of fixing range to single bp "summit"
  mid <- floor(summary(width(peak_gr))["Median"] / 2)
  start(peak_gr) <- end(peak_gr) <- start(peak_gr) + mid
  
  return(peak_gr)
  
}



# change strand information from 1/-1 to +/-

strand_to_plusminus <- function(range_table) {
  
  stopifnot("Strand" %in% names(range_table))
  
  range_table$Strand <- str_replace(as.character(range_table$Strand), "^1$", "+")
  range_table$Strand <- str_replace(as.character(range_table$Strand), "^-1$", "-")
  return(range_table)
}




# Coerce start and end of gene annotation table to just the TSS - used to fix
# TSS in protein coding tables as a single point when performing overlaps

startend_to_tss <- function(range_table) {
  
  stopifnot(c("Start", "End", "Transcription_start_site") %in% names(range_table))
  
  range_table$Start <- range_table$End <- range_table$Transcription_start_site
  
  return(range_table)
}



# Convert a gene annotation table to a GR object. If TSS is TRUE, will fix
# gene start and end coordinates to be the 1bp TSS

pc_to_gr <- function(range_table, TSS = TRUE) {
  
  range_table <- strand_to_plusminus(range_table)
  
  if (TSS) {
    range_table <- startend_to_tss(range_table)
  }

  gr <- makeGRangesFromDataFrame(range_table,
                                 keep.extra.columns = TRUE,
                                 ignore.strand = FALSE)
  
  return(gr)
}



# Convert a blacklisted regions table to a GR object

bl_to_gr <- function(range_table) {
  
  stopifnot(c("Chromosome", "Start", "End") %in% names(range_table))
  
  gr <- makeGRangesFromDataFrame(range_table, 
                                 keep.extra.columns = TRUE, 
                                 ignore.strand = TRUE)
  return(gr)
  
}



# Given GR objects corresponding to a peak table and a blacklist table,
# return peak_gr with any ranges overlapping bl_gr removed

filter_blacklist <- function(peak_gr, bl_gr) {
  
  hits <- suppressWarnings(
    findOverlaps(
      query = peak_gr,
      subject = bl_gr,
      ignore.strand = TRUE,
      type = "any",
      select = "all"
    )
  )
  if (length(hits) > 0) {
    peak_gr <- peak_gr[-hits@from]
  }
  
  return(peak_gr)
}



# Generate a gene binding score of peak-gene distance using an exponential decay
# function proposed in Ouyang et al., 2009 https://www.pnas.org/content/106/51/21521
# The original formulation scaled the score by the MACS2 score, which is 
# excluded here following contemporary practices
# distance: A vector of integers of the basepairs between a gene TSS and peak summits.
# decay_constant: An integer controlling how steeply the score decreases
# returns: An integer

ouyang <- function(distance, decay_constant = 5e3) {
  
  stopifnot(is.numeric(distance), length(distance) > 0)
  
  scores <- lapply(distance, function(x) {
    exp(-(abs(x) / decay_constant))
  })
  
  return(sum(unlist(scores)))
}



# Generate binding scores for every gene from a provided peak GR object.
# pc_gr: A GR object of protein coding gene with range fixed to the 1bp TSS
# peak_gr: A GR object of the peaks with range fixed to the 1bp summit
# max_dist: An integer specifying the cutoff (in bps) of peaks to consider
# ncore: An integer of how many cores parallel will use
# returns a numeric vector the length of pc_gr of the gene binding scores

binding_scores <- function(pc_gr, 
                           peak_gr, 
                           max_dist = 1e6, 
                           ncore = 1) {
  
  stopifnot(class(pc_gr) == "GRanges", class(peak_gr) == "GRanges")

  score_l <- mclapply(1:length(pc_gr), function(x) {
    
    dist <- GenomicRanges::distance(pc_gr[x], peak_gr, select = "all")
    dist <- dist[!is.na(dist) & dist < max_dist]
    
    if (length(dist) == 0) {
      return(0)
    }
  
    score <- ouyang(dist)
    
  }, mc.cores = ncore)
  
  names(score_l) <- pc_gr$Symbol
  
  return(unlist(score_l))
}



# Given an input df that has the file names of Unibind experiments, 
# iteratively load and gene bind score each experiment then bind into a matrix

load_and_score <- function(dir,  # top level directory of unibind experiments
                           input_df,  
                           bl_gr, # GR object of blacklisted regions
                           pc_gr, # GR object of protein coding genes
                           ncores = 1) {
  
  stopifnot("File" %in% colnames(input_df))
  
  bscore_l <- mclapply(1:nrow(input_df), function(x) {
    
    path <- file.path(dir, input_df$Symbol[x], input_df$File[x])
    
    tryCatch({
      
      message(paste(input_df$File[x], Sys.time()))
      peak_gr <- load_unibind(path) %>% unibind_to_gr(bl_gr)
      binding_scores(pc_gr = pc_gr, peak_gr = peak_gr, ncore = 1)
      
    }, error = function(e) NULL)
    
  }, mc.cores = ncores)
  
  names(bscore_l) <- meta$File
  bscore_mat <- do.call(cbind, bscore_l)
  
  return(bscore_mat)
}
