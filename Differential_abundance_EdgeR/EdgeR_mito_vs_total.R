#EdgeR for differential abundance for tRNAs in the P. citri mitochondrial samples vs. the whole insect samples. 

  library(edgeR)
  library(dplyr)
  library(ggplot2)

setwd("/Users/home/Documents/Mealybug_tRNA/Mitochondrial_tRNAs/Manuscript/Combined_Parsing")

# File: 75_strandcountmatrix_justTCME.txt
dat_raw <- read.delim("75_strandcountmatrix_justTCME.txt", check.names = FALSE)

# Keep the Library names as rownames and drop the column
stopifnot(colnames(dat_raw)[1] == "Library")
rownames(dat_raw) <- dat_raw$Library
dat_raw <- dat_raw[ , -1, drop = FALSE]
# rounded nearest integer (common, safe for small fractions).

count_mat <- round(as.matrix(dat_raw))

# Define groups from rownames: Mito_* vs Total_*
sample_names <- rownames(count_mat)
group <- ifelse(grepl("^Mito_", sample_names, ignore.case = FALSE), "Mito", "Total")
group <- factor(group, levels = c("Total", "Mito"))  # set reference level = Total

# DGEList, filtering, normalization

y <- DGEList(counts = t(count_mat))  # transpose so genes are rows, samples are columns
# Filter low-expression genes (edgeR recommended)
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]

# TMM normalization
y <- calcNormFactors(y, method = "TMM")


# Design, dispersion, fit, test

design <- model.matrix(~ group)
y <- estimateDisp(y, design = design)

fit <- glmQLFit(y, design = design)
qlf <- glmQLFTest(fit, coef = "groupMito")  # Mito vs Total

res <- topTags(qlf, n = Inf)$table
# res has: logFC (Mito vs Total), logCPM (base-2), F, PValue, FDR
res <- res %>% tibble::rownames_to_column(var = "Gene")

# Parse Genome from gene names (prefix before first "-")
res <- res %>%
  mutate(
    Genome = sub("^([^\\-]+)-.*$", "\\1", Gene),
  )

# cpm() returns CPM per sample; take mean across samples for each gene
cpm_all <- cpm(y, normalized.lib.sizes = TRUE)  # genes x samples
avg_cpm <- rowMeans(cpm_all)
res$avgCPM <- avg_cpm[match(res$Gene, rownames(cpm_all))]

# 6) Plot: logFC (y) vs expression (x), colored by Genome

p <- ggplot(res, aes(x = logCPM, y = logFC, color = Genome)) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_point(alpha = 0.9, size = 2) +
  labs(
    title = "Differential Expression: Mito vs Total (edgeR)",
    x = "logCPM (log2 counts per million)",
    y = "logFC (Mito / Total)",
    color = "Genome"
  ) +
  theme_minimal(base_size = 12)
p <- p + scale_color_manual(values = c(
  Mitochondria = "#B22222",
  Moranella    = "#FFD54D",
  Nuclear      = "#cfe8ff",  # very light blue
  Tremblaya    = "#2C7A7B"
))


print(p)
