outfile <- "~/scratch/UK_biobank_depletion_ranks.gz"
url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9329122/bin/41586_2022_4965_MOESM3_ESM.gz"

if (!file.exists(outfile)) {
  options(timeout = 1000)  # default 60 was insufficient
  download.file(url, outfile)  
}
