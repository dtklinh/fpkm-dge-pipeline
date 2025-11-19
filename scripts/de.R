#!/usr/bin/env Rscript

# scripts/de.R

suppressPackageStartupMessages({
  library(limma)
  library(argparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# Argument Parsing
parser <- ArgumentParser(description = "Differential Expression on FPKM using Limma")
parser$add_argument("--input_fpkm", required=TRUE, help="Path to FPKM matrix")
parser$add_argument("--metadata", required=TRUE, help="Path to metadata CSV")
parser$add_argument("--output_table", required=TRUE, help="Path to output CSV")
parser$add_argument("--output_plot", required=TRUE, help="Path to volcano plot")
parser$add_argument("--group_col", default="condition", help="Column in metadata for grouping")
parser$add_argument("--control", default="Control", help="Name of control group")
parser$add_argument("--treatment", default="Treatment", help="Name of treatment group")

args <- parser$parse_args()

# 1. Load Data
# Assumes Gene IDs are in the first column
fpkm_data <- read_csv(args$input_fpkm) %>% 
  tibble::column_to_rownames(var = colnames(.)[1])

meta_data <- read_csv(args$metadata)

# Ensure metadata rows match fpkm columns
if (!all(colnames(fpkm_data) %in% meta_data$sample_id)) {
  stop("Error: Column names in FPKM do not match sample_ids in metadata.")
}

# Reorder metadata to match columns
meta_data <- meta_data[match(colnames(fpkm_data), meta_data$sample_id), ]

# 2. Pre-processing
# Remove lowly expressed genes (sum FPKM > 1 across samples)
keep <- rowSums(fpkm_data) > 1
fpkm_clean <- fpkm_data[keep, ]

# Log2 Transformation (adding pseudo-count of 1 to avoid log(0))
log_fpkm <- log2(fpkm_clean + 1)

# 3. Limma Analysis
group <- factor(meta_data[[args$group_col]])
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit Linear Model
fit <- lmFit(log_fpkm, design)

# Create Contrast
contrast_formula <- paste(args$treatment, args$control, sep="-")
contrasts <- makeContrasts(contrasts = contrast_formula, levels = design)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Get Results
res <- topTable(fit2, number = Inf, adjust.method = "BH")
res <- tibble::rownames_to_column(res, "GeneID")

# 4. Save Results
write_csv(res, args$output_table)

# 5. Volcano Plot
res <- res %>%
  mutate(diffexpressed = case_when(
    logFC > 1 & adj.P.Val < 0.05 ~ "UP",
    logFC < -1 & adj.P.Val < 0.05 ~ "DOWN",
    TRUE ~ "NO"
  ))

p <- ggplot(res, aes(x=logFC, y=-log10(adj.P.Val), color=diffexpressed)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  theme_minimal() +
  labs(title = paste("Volcano Plot:", args$treatment, "vs", args$control))

ggsave(args$output_plot, plot=p, width=6, height=5)

print(paste("Analysis complete. Results saved to", args$output_table))