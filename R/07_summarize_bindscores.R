##
## TODO: unify meta between matrices and GR
## -----------------------------------------------------------------------------


        


# Generate average binding scores across experiments: 1) Across all experiments;
# 2) Grouped by TF; 3) Mean of means across TFs, to summarize across all 
# experiments but to reduce effect of imbalanced TF counts.
# ------------------------------------------------------------------------------


# Return a df of the average bind score per gene across all experiments

get_all_mean <- function(mat) {
  data.frame(Symbol = rownames(mat), Mean = rowMeans(mat))
}


# Return a gene by TF matrix of the average binding scores

get_tf_mean <- function(mat, meta) {
  
  tfs <- unique(meta$Symbol)
  
  mean_l <- lapply(unique(tfs), function(x) {
    meta <- filter(meta, Symbol == x)
    rowMeans(mat[, meta$File, drop = FALSE])
  })
  
  mean_mat <- as.matrix(do.call(cbind, mean_l))
  colnames(mean_mat) <- tfs
  
  return(mean_mat)
  
}


# Mean summaries into a list

mean_l <- list(
  Human_all = get_all_mean(bind_l$Human),
  Human_TF = get_tf_mean(bind_l$Human, dat$Permissive_hg$Meta),
  Human_mom = get_all_mean(get_tf_mean(bind_l$Human, dat$Permissive_hg$Meta)),
  Mouse_all = get_all_mean(bind_l$Mouse),
  Mouse_TF = get_tf_mean(bind_l$Mouse, dat$Permissive_mm$Meta),
  Mouse_mom = get_all_mean(get_tf_mean(bind_l$Mouse, dat$Permissive_mm$Meta))
)


# plot(mean_l$Human_all$Mean, mean_l$Human_mom$Mean)
# cor(mean_l$Human_all$Mean, mean_l$Human_mom$Mean)
# plot(mean_l$Mouse_all$Mean, mean_l$Mouse_mom$Mean)
# cor(mean_l$Mouse_all$Mean, mean_l$Mouse_mom$Mean)


# Top bound genes, using mean of means.
# ------------------------------------------------------------------------------


filter_top <- function(mean_df, qtl = 0.99) {
  filter(mean_df, Mean >= quantile(Mean, qtl)) %>% arrange(desc(Mean))
}



# Human
top_bound_hg <- filter_top(mean_l$Human_mom)
top_dist_hg <- dist_hg[top_bound_hg$Symbol, ]
top_summ_hg <- dist_summ_hg[top_bound_hg$Symbol, ]

# Mouse
top_bound_mm <- filter_top(mean_l$Mouse_mom)
top_dist_mm <- dist_mm[top_bound_mm$Symbol, ]
top_summ_mm <- dist_summ_mm[top_bound_mm$Symbol, ]


# Relationship between count of DHSs around genes and mean binding


tt <- left_join(mean_l$Human_mom, pc_dhs_count, by = "Symbol")

ggplot(tt, aes(x = Count, y = Mean)) +
  geom_point(alpha = 0.4, colour = "royalblue") +
  xlab("Count of DHS elements around TSS") +
  ylab("Mean binding score") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.margin = margin(10, 20, 10, 10))

cor(tt$Mean, tt$Count, method = "spearman")


# Binding but no DHS

tt %>% 
  filter(Count == 0) %>% 
  arrange(desc(Mean)) %>% 
  head(20)

pc_dhs_count_top <- mutate(pc_dhs_count, Top = Symbol %in% top_bound_hg$Symbol)
boxplot(pc_dhs_count_top$Count ~ pc_dhs_count_top$Top)
wilcox.test(pc_dhs_count_top$Count ~ pc_dhs_count_top$Top)



# Rarely/never bound genes
# ------------------------------------------------------------------------------


filter_bottom <- function(mean_df, qtl = 0.01) {
  filter(mean_df, Mean <= quantile(Mean, qtl)) %>% arrange(Mean)
}


# TODO: finalize mat for infreq calc
mat_raw <- dat$Permissive_hg$Mat_raw
mat_qnl <- dat$Permissive_hg$Mat_qnl
mean_raw <- get_all_mean(get_tf_mean(mat_raw, dat$Permissive_hg$Meta))
mean_qnl <- get_all_mean(get_tf_mean(mat_qnl, dat$Permissive_hg$Meta))
# plot(mean_raw$Mean, mean_qnl$Mean)
# cor(mean_raw$Mean, mean_qnl$Mean)

btm_bound_hg1 <- filter_bottom(mean_raw)
btm_bound_hg2 <- filter_bottom(mean_qnl)
length(intersect(btm_bound_hg1$Symbol, btm_bound_hg2$Symbol))


dist_summ_hg[btm_bound_hg1$Symbol, ]


# Example of infreq bound gene that still has a proximal peak
dist_summ_hg["CT45A6", ]
summary(mat_raw["CT45A6", ])
head(sort(mat_raw["CT45A6", ], decreasing = TRUE))
summary(mat_qnl["CT45A6", ])
head(sort(mat_qnl["CT45A6", ], decreasing = TRUE))


# Example of infreq bound for DHS

pc_dhs_count2 <- mutate(pc_dhs_count,
                        Raw = Symbol %in% btm_bound_hg1$Symbol,
                        QNL = Symbol %in% btm_bound_hg2$Symbol)

pc_dhs_count2 %>% head
boxplot(pc_dhs_count2$Count ~ pc_dhs_count2$Raw)
boxplot(pc_dhs_count2$Count ~ pc_dhs_count2$QNL)


# For frequent bound, are peaks dispersed?

pc_sub <- pc_hg[pc_hg$Symbol == "MIDN"]

downstream <- dhs_gr[dhs_gr@seqnames == as.numeric(pc_sub@seqnames) & dhs_gr@ranges < pc_sub@ranges]
upstream <- dhs_gr[dhs_gr@seqnames == as.numeric(pc_sub@seqnames) & dhs_gr@ranges > pc_sub@ranges]
downstream$Distance <- distance(pc_sub, downstream)
upstream$Distance <- distance(pc_sub, upstream)





# Plots
# ------------------------------------------------------------------------------


# Histogram of mean binding scores

plot_hist <- function(mean_df, species, qtl_upper = 0.99) {
  
  if (species == "Human") {
    fill <- "royalblue"
  } else if (species == "Mouse") {
    fill <- "goldenrod"
  }
  
  ggplot(mean_df, aes(x = Mean)) +
    geom_histogram(alpha = 0.6, fill = fill, bins = 100) +
    geom_vline(xintercept = quantile(mean_df$Mean, qtl_upper), 
               linewidth = 1.2,
               colour = "firebrick3") +
    theme_classic() +
    ylab("Density") +
    xlab("Mean binding score") +
    ggtitle(species) +
    theme(
      axis.text = element_text(size = 25),
      axis.title = element_text(size = 25),
      plot.title = element_text(hjust = 0.5, size = 30),
      plot.margin = margin(10, 20, 10, 10) 
    )
}


p1a <- plot_hist(mean_l$Human_mom, "Human")
p1b <- plot_hist(mean_l$Mouse_mom, "Mouse")
