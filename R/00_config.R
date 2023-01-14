## Establishing pathing
## -----------------------------------------------------------------------------


dat_dir <- "/space/grp/amorin/"


# Path of compressed Unibind download
rob_path_hg <- paste0(dat_dir, "Peak_files/Unibind/Hg38/Robust/damo_hg38_TFBS_per_TF")
rob_path_mm <- paste0(dat_dir, "Peak_files/Unibind/Mm10/Robust/damo_mm10_TFBS_per_TF")
perm_path_hg <- paste0(dat_dir, "Peak_files/Unibind/Hg38/Permissive/damo_hg38_TFBS_per_TF")
perm_path_mm <- paste0(dat_dir, "Peak_files/Unibind/Mm10/Permissive/damo_mm10_TFBS_per_TF")


# TODO: formalize package or download functions
ref_hg <- paste(dat_dir, "Metadata/refseq_select_hg38.tsv")
ref_mm <- paste0(dat_dir, "Metadata/refseq_select_mm10.tsv")

bl_hg <-  paste0(dat_dir, "Chromosome_info/blacklist_hg38.tsv")
bl_mm <- paste0(dat_dir, "Chromosome_info/blacklist_hg38.tsv")


# Path of UK biobank depletion ranks
dr_path <- "/space/scratch/amorin/UK_biobank_depletion_ranks.gz"
