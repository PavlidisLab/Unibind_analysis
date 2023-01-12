## Download Unibind data, ...
## TODO: port over pcoding/ccre/etc table calls from TR agg OR make package
## TODO: robust/permissive def
## -----------------------------------------------------------------------------

source("R/00_config.R")

# Unibind data: Robust versus Permissive collection of peaks are each saved
# in a single directory for all experiments.
# Permissive:
# Robust: 


rob_link_mm <- "https://testunibind.uio.no/static/data/bulk_Robust/Mus_musculus/mm10_TFBS_per_TF.tar.gz"
perm_link_mm <- "https://testunibind.uio.no/static/data/bulk_Permissive/Mus_musculus/mm10_TFBS_per_TF.tar.gz"
rob_link_hg <- "https://testunibind.uio.no/static/data/bulk_Robust/Homo_sapiens/hg38_TFBS_per_TF.tar.gz"
perm_link_hg <- "https://testunibind.uio.no/static/data/bulk_Permissive/Homo_sapiens/hg38_TFBS_per_TF.tar.gz"


input_df <- data.frame(
  link = c(rob_mm_link, perm_mm_link, rob_hg_link, perm_hg_link),
  path = c(rob_mm_path, perm_mm_path, rob_hg_path, perm_hg_path),
  stringsAsFactors = FALSE
)





# UK Biobank depletion ranks

outfile <- "~/scratch/UK_biobank_depletion_ranks.gz"
url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9329122/bin/41586_2022_4965_MOESM3_ESM.gz"

if (!file.exists(outfile)) {
  options(timeout = 1000)  # default 60 was insufficient
  download.file(url, outfile)  
}
