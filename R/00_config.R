## Establishing pathing. NOTE: A lot of these paths were established in
## preceding projects, hence the rather disjointed directory structure...
## -----------------------------------------------------------------------------


# This is the root dir where the raw peak files are stored
peak_dir <- "/space/grp/amorin/Peak_files/Unibind"

# This is the root dir where processed output data was dumped
dat_dir <- "/space/scratch/amorin/R_objects"

# This is the root dir where chromosome data like blacklists are stored
chrom_dir <- "/space/grp/amorin/Chromosome_info"

# This is the root dir where information like protein coding tables are stored
meta_dir <- "/space/grp/amorin/Metadata"

# For use in parallel
cores <- 8

# Minimum peaks to be considered for analysis
min_peaks <- 100

# Path of compressed robust and permissive unibind download
rob_path_hg <- file.path(peak_dir, "Hg38/Robust/damo_hg38_TFBS_per_TF")
rob_path_mm <- file.path(peak_dir, "Mm10/Robust/damo_mm10_TFBS_per_TF")
perm_path_hg <- file.path(peak_dir, "Hg38/Permissive/damo_hg38_TFBS_per_TF")
perm_path_mm <- file.path(peak_dir, "Mm10/Permissive/damo_mm10_TFBS_per_TF")

# Raw bindscore matrices
bmat_path_hg <- file.path(dat_dir, "unibind_all_scores_hg.RDS")
bmat_path_mm <- file.path(dat_dir, "unibind_all_scores_mm.RDS")

# List of processed matrices and metadata
bind_dat_path <- file.path(dat_dir, "processed_unibind_data.RDS")
meta_path <- file.path(dat_dir, "unibind_metadata.RDS")

# Average bind scores and output of binding specificity model
bind_summary_path <- file.path(dat_dir, "unibind_bindscore_summary.RDS")
bind_model_path <- file.path("unibind_bindscore_modelfit.RDS")

# Experiments saved as list of GRange objects
gr_perm_path_hg <- file.path(dat_dir, "unibind_grlist_perm_human.RDS")
gr_rob_path_hg <- file.path(dat_dir, "unibind_grlist_rob_human.RDS")
gr_perm_path_mm <- file.path(dat_dir, "unibind_grlist_perm_mouse.RDS")
gr_rob_path_mm <- file.path(dat_dir, "unibind_grlist_rob_mouse.RDS")

# Refseq select protein coding tables used for all gene-based analysis
ref_path_hg <- file.path(meta_dir, "refseq_select_hg38_jan2024.tsv")
ref_path_mm <- file.path(meta_dir, "refseq_select_mm10_jan2024.tsv")

# ENCODE black lists, used to filter peaks in blacklisted regions
bl_path_hg <-  file.path(chrom_dir, "blacklist_hg38.tsv")
bl_path_mm <- file.path(chrom_dir, "blacklist_mm10.tsv")

