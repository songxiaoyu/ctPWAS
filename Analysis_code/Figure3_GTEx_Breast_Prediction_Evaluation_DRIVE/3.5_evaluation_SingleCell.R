
rm(list = ls())

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(SingleCellExperiment)
library(Matrix)

setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/Density_Song/Paper_PWAS")

paper_dir <- getwd()
data_dir <- file.path(paper_dir, "Data")
results_dir <- file.path(paper_dir, "Results")
figure_dir <- file.path(paper_dir, "Figure")

# --------------------------------------------------
# 1. Load snRNA-seq data and prediction models
# --------------------------------------------------

a <- readRDS(file.path( data_dir, "GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.RData"))

b <- a[, !is.na(a$tissue) & a$tissue == "breast"]

load(file.path(  results_dir,"GTEx_Breast_pi2_ThreeModel_Summary.Rdata"))

ensembl38 <- read_csv( file.path(data_dir, "ensembl38.txt"),  show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  dplyr::select(gene_stable_id, gene_name) %>%
  distinct(gene_stable_id, .keep_all = TRUE)

# Use both models
models <- list(S = filtered_weights_s,M = filtered_weights_m_old)

models <- lapply(models, function(model) {
  model %>%
    mutate(gene_stable_id = sub("\\..*$", "", gene)) %>%
    left_join(ensembl38, by = "gene_stable_id")
})

# Restrict to genes present in both models
common_genes <- Reduce(intersect, lapply(models, function(x) {unique(x$gene_name[!is.na(x$gene_name)])}))
common_genes <- intersect(common_genes, rownames(b))
c <- b[rownames(b) %in% common_genes, ]

# --------------------------------------------------
# 2. Calculate observed epithelial snRNA-seq means
# --------------------------------------------------

c$`Participant ID` <- droplevels(c$`Participant ID`)

cellinfo <- colData(c)
geneinfo <- rowData(c)

IDs <- c("GTEX-1R9PN", "GTEX-1CAMS", "GTEX-1MCC2")

# Only epithelial cells from the three individuals
keep_cells <- c$`Participant ID` %in% IDs & c$`Tissue composition` == "Epithelial"
c_epi <- c[, keep_cells]

# Compute mean logcounts by individual
logcount <- assay(c_epi, "logcounts")

snMean <- sapply(IDs, function(id) {
  rowMeans(logcount[, c_epi$`Participant ID` == id, drop = FALSE])
})

rownames(snMean) <- rownames(c_epi)
colnames(snMean) <- paste0(IDs, "_Epithelial")

# Rank-based inverse normal transformation
d_norm <- apply(snMean, 2, function(x) {
  qnorm((rank(x, ties.method = "average") - 0.5) /length(x))
})

rownames(d_norm) <- rownames(snMean)
colnames(d_norm) <- colnames(snMean)

# --------------------------------------------------
# 3. Load genotype data
# --------------------------------------------------

geno1 <- fread(file.path(data_dir,"GTEx_Breast/shapeit_data_for_predictdb_variants-r2"))
geno <- as.data.frame(geno1[, c("ID", IDs), with = FALSE])

# Convert genotype columns once
geno[IDs] <- lapply(geno[IDs], function(x) {
  as.numeric(replace(as.character(x),
                     as.character(x) == ".", NA))
})

geno_mat <- as.matrix(geno[, IDs, drop = FALSE])
rownames(geno_mat) <- geno$ID

# --------------------------------------------------
# 4. Predict expression for each model
# --------------------------------------------------

predict_expression <- function(model) {

  model2 <- model %>% filter(gene_name %in% common_genes, !is.na(gene_name))

  # Keep SNPs available in genotype data
  model2 <- model2[model2$varID %in% rownames(geno_mat), ,drop = FALSE ]

  # Only genes with at least one matched SNP
  GeneID <- unique(model2$gene)

  snp_idx <- match(model2$varID, rownames(geno_mat))
  gene_idx <- match(model2$gene, GeneID)

  G <- geno_mat[snp_idx, , drop = FALSE]

  # Sparse gene-by-SNP weight matrix
  W <- Matrix::sparseMatrix(
    i = gene_idx,
    j = seq_len(nrow(model2)),
    x = model2$Cell1,
    dims = c(length(GeneID), nrow(model2))
  )

  # Genes x individuals
  P <- as.matrix(W %*% G)

  annotation <- model2[match(GeneID, model2$gene),c("gene", "gene_stable_id", "gene_name") ]

  result <- cbind(annotation, as.data.frame(P))
  colnames(result)[4:6] <- IDs

  result
}

# Epithelial prediction for each method
epithelial_s <- predict_expression(models$S)
epithelial_m <- predict_expression(models$M)

# --------------------------------------------------
# 5. Compare predictions against observed snRNA-seq
# --------------------------------------------------

top <- 100

comparisons <- list(
  "White vs AA" = c("GTEX-1R9PN", "GTEX-1MCC2"),
  "White vs Asian" = c("GTEX-1R9PN", "GTEX-1CAMS"),
  "AA vs Asian" = c("GTEX-1MCC2", "GTEX-1CAMS")
)

analyze_comparison <- function(epithelial, comparison,
                               method, top = 100) {

  id1 <- comparison[1]
  id2 <- comparison[2]

  # Predicted differences
  diff <- epithelial[[id1]] - epithelial[[id2]]

  # Keep genes with observed expression in both individuals
  gene <- epithelial$gene_name

  valid <- !is.na(gene) &
    gene %in% rownames(d_norm) &
    is.finite(diff)

  dat_pred <- data.frame(
    gene_name = gene[valid],
    diff = diff[valid]
  )

  # One entry per gene name
  dat_pred <- dat_pred %>%
    arrange(desc(abs(diff))) %>%
    distinct(gene_name, .keep_all = TRUE)

  # Top positive and negative predictions
  pos <- dat_pred %>%
    filter(diff > 0) %>%
    arrange(desc(diff)) %>%
    slice_head(n = top)

  neg <- dat_pred %>%
    filter(diff < 0) %>%
    arrange(diff) %>%
    slice_head(n = top)

  # Observed differences for selected genes
  observed_diff <- d_norm[, paste0(id1, "_Epithelial")] -
    d_norm[, paste0(id2, "_Epithelial")]

  dat1 <- observed_diff[pos$gene_name]
  dat2 <- observed_diff[neg$gene_name]

  dat1 <- dat1[is.finite(dat1)]
  dat2 <- dat2[is.finite(dat2)]

  # One-sided Wilcoxon signed-rank tests
  p1 <- if (length(dat1) > 0) {
    wilcox.test(dat1, alternative = "greater")$p.value
  } else NA_real_

  p2 <- if (length(dat2) > 0) {
    wilcox.test(dat2, alternative = "less")$p.value
  } else NA_real_

  dat <- bind_rows(
    data.frame(Group = "Predicted higher", Value = dat1),
    data.frame(Group = "Predicted lower", Value = dat2)
  )

  dat$Method <- method
  dat$Comparison <- paste(
    c("White vs AA", "White vs Asian", "AA vs Asian")[
      vapply(comparisons, function(x) identical(x, comparison),
             logical(1))
    ]
  )

  list(
    data = dat,
    stats = data.frame(
      Method = method,
      Comparison = dat$Comparison[1],
      P_higher = p1,
      P_lower = p2,
      N_higher = length(dat1),
      N_lower = length(dat2),
      Mean_higher = if (length(dat1)) mean(dat1) else NA_real_,
      Mean_lower = if (length(dat2)) mean(dat2) else NA_real_
    )
  )
}

# Analyze all six combinations
results <- list()

for (method in names(models)) {

  epithelial <- if (method == "S") {
    epithelial_s
  } else {
    epithelial_m
  }

  for (comp in names(comparisons)) {

    key <- paste(method, comp, sep = "_")

    results[[key]] <- analyze_comparison(
      epithelial = epithelial,
      comparison = comparisons[[comp]],
      method = method,
      top = top
    )
  }
}


# Save analysis results
save(results, file = file.path(results_dir, "snRNAseq_comparison_S_vs_M.RData"))
