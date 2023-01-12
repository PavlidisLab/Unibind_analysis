## Downloading mouse and human ChIP-seq BED files from Unibind
## https://testunibind.uio.no/downloads/

library(tidyverse)

# There are two bulk download sets for each species. Permissive and Robust, 
# where Robust has had the DAMO algorithm applied to find motif enriched peaks


rob_mm_link <- "https://testunibind.uio.no/static/data/bulk_Robust/Mus_musculus/mm10_TFBS_per_TF.tar.gz"
rob_mm_path <- "~/Data/Peak_files/Unibind/Mm10/Robust/mm10_TFBS_per_TF.tar.gz"

perm_mm_link <- "https://testunibind.uio.no/static/data/bulk_Permissive/Mus_musculus/mm10_TFBS_per_TF.tar.gz"
perm_mm_path <- "~/Data/Peak_files/Unibind/Mm10/Permissive/mm10_TFBS_per_TF.tar.gz"

rob_hg_link <- "https://testunibind.uio.no/static/data/bulk_Robust/Homo_sapiens/hg38_TFBS_per_TF.tar.gz"
rob_hg_path <- "~/Data/Peak_files/Unibind/Hg38/Robust/hg38_TFBS_per_TF.tar.gz"

perm_hg_link <- "https://testunibind.uio.no/static/data/bulk_Permissive/Homo_sapiens/hg38_TFBS_per_TF.tar.gz"
perm_hg_path <- "~/Data/Peak_files/Unibind/Hg38/Permissive/hg38_TFBS_per_TF.tar.gz"


input_df <- data.frame(
  link = c(rob_mm_link, perm_mm_link, rob_hg_link, perm_hg_link),
  path = c(rob_mm_path, perm_mm_path, rob_hg_path, perm_hg_path),
  stringsAsFactors = FALSE
)


# Because of the structure of the unzipped files, an error I was getting with 
# untar() when trying to point to a destination directory, and that there are 
# just 4 overall directories, I manually restructured the downloaded dirs into
# Data/Peak_files/Unibind



download.file(input_df$link[3], input_df$path[3])


untar(tarfile = input_df$path[4])


# Because the permissive/robust have the same dir naming, and because untar
# is giving an error when I try and point to a non-default directory, can't
# loop because it will overwrite the untar'd dirs.
# for (i in 1:nrow(input_df)) {
#   download.file(input_df$link[i], input_df$path[i])
#   untar(tarfile = "~/Data/Peak_files/Unibind/mm10_TFBS_per_TF.tar.gz",
#         exdir = "~/Data/Peak_files/Unibind/")
# }

