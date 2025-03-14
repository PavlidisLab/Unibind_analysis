## Download Unibind data
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")

options(timeout = 1000)  # default 60 was insufficient


# Unibind data: Robust versus Permissive collection of peaks are each saved
# in a single directory for all experiments. Permissive collection is all 
# experiments; robust are experiments where the canonical motif is enriched
# near the summits of the given experiment's peaks. Note that experiments
# present in both collections are identical -- robust does not filter away peaks,
# just datasets.
# https://unibind.uio.no/docs/
# ------------------------------------------------------------------------------


rob_link_hg <- "https://unibind.uio.no/static/data/20220914/bulk_Robust/Homo_sapiens/damo_hg38_TFBS_per_TF.tar.gz"
rob_link_mm <- "https://unibind.uio.no/static/data/20220914/bulk_Robust/Mus_musculus/damo_mm10_TFBS_per_TF.tar.gz"
perm_link_hg <- "https://unibind.uio.no/static/data/20220914/bulk_Permissive/Homo_sapiens/damo_hg38_TFBS_per_TF.tar.gz"
perm_link_mm <- "https://unibind.uio.no/static/data/20220914/bulk_Permissive/Mus_musculus/damo_mm10_TFBS_per_TF.tar.gz"




input_df <- data.frame(
  link = c(rob_link_mm, perm_link_mm, rob_link_hg, perm_link_hg),
  path = c(rob_path_mm, perm_path_mm, rob_path_hg, perm_path_hg),
  stringsAsFactors = FALSE
)


input_df$path <- paste0(input_df$path, ".tar.gz")


# NOTE: The untar creates a dir with another dir one level down with identical 
# name, which then contains the individual TF dirs. I ended up just manually 
# moving these dirs a level up. 


for (i in 1:nrow(input_df)) {
  
  if (!file.exists(input_df$path[i])) {
    
    download.file(input_df$link[i], input_df$path[i])
    
    untar(tarfile = input_df$path[i], 
          exdir = str_replace(input_df$path[i], ".tar.gz", ""))
  }
}



# Refseq select protein coding genes
# Using ensembl range formatting: no 'chr' prefix and strand as 1/-1
# https://www.ncbi.nlm.nih.gov/refseq/refseq_select/
# Table schema (mouse) http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1827432270_fzD31oCmH8TcANzQ82Yzz7h0RKlz&hgta_doSchemaDb=mm10&hgta_doSchemaTable=ncbiRefSeqSelect
# --- Raw columns
# bin	
# V2 == name == Refseq_ID	
# V3 == chrom == Chromsome	
# V4 == strand == Strand	
# V5 == txStart == Start	
# V6 == txEnd == End	
# cdsStart	
# cdsEnd	
# exonCount	
# exonStarts	
# exonEnds	
# score	
# V13 == name2 == Symbol	
# cdsStartStat	
# cdsEndStat	
# exonFrames
# ------------------------------------------------------------------------------


download_refseq <- function(outfile, 
                            species) {  # Human|Mouse
  
  if (file.exists(outfile)) return(message(outfile, " already exists!"))
  
  if (species == "Mouse") {
    link <- "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/ncbiRefSeqSelect.txt.gz"
    chr <- c(1:19, "MT", "X", "Y")
  } else if (species == "Human") {
    link <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz"
    chr <- c(1:22, "MT", "X", "Y")
  }
  
  download.file(link, outfile)
  
  refseq <- read.delim(outfile, stringsAsFactors = FALSE, header = FALSE)
  
  refseq <- dplyr::select(refseq, c(V3, V5, V6, V4, V2, V13))
  
  colnames(refseq) <- c("Chromosome",
                        "Start",
                        "End",
                        "Strand",
                        "Refseq_ID",
                        "Symbol")
  
  refseq <- refseq %>% 
    mutate(
      Strand = ifelse(Strand == "+", 1, -1),
      Transcription_start_site = ifelse(Strand == 1, Start, End),
      Chromosome = str_replace(Chromosome, "chr", "")) %>% 
    filter(Chromosome %in% chr) %>% 
    arrange(match(Chromosome, chr), Transcription_start_site) %>% 
    dplyr::relocate(Transcription_start_site, .after = Chromosome)
  
  
  write.table(refseq,
              quote = FALSE,
              row.names = FALSE,
              sep = "\t",
              file = outfile)
  
}



download_refseq(outfile = ref_path_hg, species = "Human")
download_refseq(outfile = ref_path_mm, species = "Mouse")


# ENCODE blacklisted regions
# ------------------------------------------------------------------------------


chr_hg <- c(1:22, "MT", "X", "Y")
chr_mm <- c(1:19, "MT", "X", "Y")

bl_url_hg <- "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
bl_url_mm <- "https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"


download_blacklist <- function(outfile, url, chr) {
  
  if (!file.exists(outfile)) {
    
    download.file(url, destfile = outfile)
    bl <- read.delim(outfile, header = FALSE)
    colnames(bl) <- c("Chromosome", "Start", "End")
    bl$Chromosome <- gsub("chr", "", bl$Chromosome)
    bl <- bl[order(match(bl$Chromosome, chr)),]
    
    write.table(bl,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                file = outfile)
  }
  
}


download_blacklist(bl_path_hg, bl_url_hg, chr_hg)
download_blacklist(bl_path_mm, bl_url_mm, chr_mm)
