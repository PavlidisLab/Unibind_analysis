## Establishing pathing
## -----------------------------------------------------------------------------


dat_dir <- "/space/grp/amorin/"

# For use in parallel
cores <- 8

# Minimum peaks to be considered for analysis
min_peaks <- 100

# Unibind metadata
# TODO: reconsider pathing
meta_outfile <- "~/scratch/R_objects/Unibind_metadata.RDS"

# Path of compressed Unibind download
rob_path_hg <- paste0(dat_dir, "Peak_files/Unibind/Hg38/Robust/damo_hg38_TFBS_per_TF")
rob_path_mm <- paste0(dat_dir, "Peak_files/Unibind/Mm10/Robust/damo_mm10_TFBS_per_TF")
perm_path_hg <- paste0(dat_dir, "Peak_files/Unibind/Hg38/Permissive/damo_hg38_TFBS_per_TF")
perm_path_mm <- paste0(dat_dir, "Peak_files/Unibind/Mm10/Permissive/damo_mm10_TFBS_per_TF")

# Raw bindscore matrices
bmat_path_hg <- "/space/scratch/amorin/R_objects/unibind_all_scores_hg.RDS"
bmat_path_mm <- "/space/scratch/amorin/R_objects/unibind_all_scores_mm.RDS"

# Processed list of meta and matrices
bind_dat_path <- "/space/scratch/amorin/R_objects/processed_unibind_data.RSD"

# Experiments saved as list of GRange objects
gr_path_perm_hg <- "/space/scratch/amorin/R_objects/unibind_grlist_perm_human.RDS"
gr_path_rob_hg <- "/space/scratch/amorin/R_objects/unibind_grlist_rob_human.RDS"
gr_path_perm_mm <- "/space/scratch/amorin/R_objects/unibind_grlist_perm_mouse.RDS"
gr_path_rob_mm <- "/space/scratch/amorin/R_objects/unibind_grlist_rob_mouse.RDS"

# TODO: formalize package or download functions
ref_hg <- paste0(dat_dir, "Metadata/refseq_select_hg38.tsv")
ref_mm <- paste0(dat_dir, "Metadata/refseq_select_mm10.tsv")
bl_path_hg <-  paste0(dat_dir, "Chromosome_info/blacklist_hg38.tsv")
bl_path_mm <- paste0(dat_dir, "Chromosome_info/blacklist_hg38.tsv")

# Path of UK biobank depletion ranks
dr_path <- "/space/scratch/amorin/UK_biobank_depletion_ranks.gz"
