## Functions used throughout analysis
## -----------------------------------------------------------------------------

library(GenomicRanges)

## TODO: range table functions -> package

source("/home/amorin/Projects/TR_aggregation/R/utils/range_table_functions.R")



load_unibind <- function(path) {
  
  peak <- read.delim(
    file = path,
    header = FALSE,
    col.names = c("Chromosome", "Start", "End", "Motif", "Score", "Strand"),
    sep = "\t"
  )
  
  return(peak)
}



unibind_to_gr <- function(peak, bl_gr) {
  
  peak_gr <- peak %>%
    get_standard_chr() %>%
    peak_to_gr(format = FALSE) %>%
    filter_blacklist(bl_gr)
  
  # Hacky means of fixing range to single bp "summit"
  mid <- floor(summary(width(peak_gr))["Median"] / 2)
  start(peak_gr) <- end(peak_gr) <- start(peak_gr) + mid
  
  return(peak_gr)
  
}



# TODO: internal bind_score() call also uses parallel, need to test if quicker
# to parallel the outer mclapply call or the binding_scores() call.

load_and_score <- function(dir, input_df, bl_gr, pc_gr, ncores = 1) {
  
  scores <- mclapply(1:nrow(input_df), function(x) {
    
    path <- paste0(dir, "/", input_df$Symbol[x], "/", input_df$File[x])
    
    tryCatch({
      
      peak_gr <- load_unibind(path) %>% unibind_to_gr(bl_gr)
      
      
      bscore <- binding_scores(pc_gr = pc_gr, 
                               peak_gr = peak_gr, 
                               method = "Ouyang", 
                               ncore = 1)
      
      message(input_df$File[x], " complete ", Sys.time())
      
      return(bscore)
      
    }, error = function(e) NULL)
    
    
    
  }, mc.cores = ncores)
  
  return(scores)
}



# Assumes m x n mat with named rows and columns. Returns a dataframe of all
# the row-col elements - if symmetric is TRUE, then only keep unique pairs
# https://stackoverflow.com/questions/28035001/

mat_to_df <- function(mat, symmetric = TRUE) {
  
  if (symmetric) {
    df <- data.frame(
      Row = rownames(mat)[row(mat)[lower.tri(mat)]],
      Col = colnames(mat)[col(mat)[lower.tri(mat)]],
      Value = mat[lower.tri(mat)],
      stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      Row = rownames(mat)[row(mat)],
      Col = colnames(mat)[col(mat)],
      Value = c(mat),
      stringsAsFactors = FALSE
    )
  }
  return(df)
}



tri_cormat <- function(cor_mat) {
  
  hc <- hclust(as.dist(1 - cor_mat))
  cor_mat_tri <- cor_mat[hc$order, hc$order]
  cor_mat_tri[lower.tri(cor_mat_tri)] <-  NA
  diag(cor_mat_tri) <- NA
  cor_mat_tri <- cor_mat_tri[1:nrow(cor_mat_tri) - 1, 2:ncol(cor_mat_tri)]
  
  return(cor_mat_tri)
}
