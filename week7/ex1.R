library(tidyverse)
library(broom)
library(DESeq2)

# set your working directory
setwd("~/qbb2024-answers/week7")

# load the gene expression counts
counts_df <- read_delim("gtex_whole_blood_counts_downsample.txt")
# move the gene_id column to rownames, so that the contents of the tibble is entirely numeric
counts_df <- column_to_rownames(counts_df, var = "GENE_NAME")

# load the metadata
metadata_df <- read_delim("gtex_metadata_downsample.txt")
# move the sample IDs from the first column to rownames
metadata_df <- column_to_rownames(metadata_df, var = "SUBJECT_ID")

# look at the first few rows and columns
counts_df[1:5,]
metadata_df[1:5,]

# check that the columns of the counts are identical and in the same order as the rows of the metadata
colnames(counts_df) == rownames(metadata_df)
table(colnames(counts_df) == rownames(metadata_df))
# all True

# create a DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts_df,
                              colData = metadata_df,
                              design = ~ SEX + AGE + DTHHRDY)
# apply VST normalization
vsd <- vst(dds)

# apply and plot principal components
PCA_SEX = plotPCA(vsd, intgroup = "SEX")
PCA_AGE = plotPCA(vsd, intgroup = "AGE")
PCA_DTHHRDY = plotPCA(vsd, intgroup = c("DTHHRDY"))
ggsave("PCA_SEX.png",PCA_SEX)
ggsave("PCA_AGE.png",PCA_AGE)
ggsave("PCA_DTHHRDY.png",PCA_DTHHRDY)

# What proportion of variance in the gene expression data is explained by each of the first two principal components? 
# 7% and 48%
# Which principal components appear to be associated with which subject-level variables?
# DTHHRDY (cause of death)

# extract the normalized expression data and bind to metadata
vsd_df <- assay(vsd) %>%
  t() %>%
  as_tibble()
vsd_df <- bind_cols(metadata_df, vsd_df)

# test for differential expression of the gene WASH7P
m1 <- lm(formula = WASH7P ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()
# Does WASH7P show significant evidence of sex-differential expression (and if so, in which direction)? 
# No. p value for SEX expression is 2.792437e-01, which is larger than 0.05. This indicates that the result is not statistically significant.

# test for differential expression of the gene SLC25A47
m2 <- lm(formula = SLC25A47 ~ DTHHRDY + AGE + SEX, data = vsd_df) %>%
  summary() %>%
  tidy()
# Does SLC25A47 show significant evidence of sex-differential expression (and if so, in which direction)? 
# Yes.  p value for SEX expression is 2.569926e-02, which is smaller than 0.05. This indicates that the result is statistically significant.
# male tend to have higher expression level

# Perform differential expression analysis
dds <- DESeq(dds)
view(dds)

# Extract the differential expression results for the variable SEX.
SEX_res <- results(dds, name = "SEX_male_vs_female")  %>%
  as_tibble(rownames = "GENE_NAME")

# The rows with a padj < 0.1
SEX_res <- SEX_res %>%
  filter(padj<0.1) %>%
  arrange(padj)
# How many genes exhibit significant differential expression between males and females at a 10% FDR?
# 262

# load the mappings of genes to chromosomes
gene_location = read.delim("gene_locations.txt")
# Merge these data with your differential expression 
merge_data = left_join(gene_location,SEX_res,by="GENE_NAME") %>% 
  arrange(padj) # Order the merged tibble by padj, from smallest to largest.
# Which chromosomes encode the genes that are most strongly upregulated in males versus females, respectively?
# Y and X chromosome (Sex chromosome)
# Are there more male-upregulated genes or female-upregulated genes near the top of the list?
# more male-related genes are upregulated near the top of the list.

# test for differential expression of the gene WASH7P and SLC25A47
WASH7P = merge_data %>%
  filter(GENE_NAME == "WASH7P")
SLC25A47 = merge_data %>%
  filter(GENE_NAME == "SLC25A47")
# consistent

# Perform differential expression analysis on cause of death

# extract the data 
DTHHRDY_res <- results(dds, name = "DTHHRDY_ventilator_case_vs_fast_death_of_natural_causes")  %>%
  as_tibble(rownames = "GENE_NAME")

# The rows with a padj < 0.1
DTHHRDY_res <- DTHHRDY_res %>%
  filter(padj<0.1) %>%
  arrange(padj)
# How many genes exhibit significant differential expression between males and females at a 10% FDR?
# 16069
# it makes sense. PCA plot shows that causes of death is more associated with differential expression levels. Thus more genes are differentially expressed at an FDR of 10% by death classification.

library(ggplot2)

VolcanoPlot = ggplot(data = SEX_res,
      aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = (abs(log2FoldChange) > 1 & -log10(padj)>1))) +
  geom_text(data = SEX_res %>% filter(abs(log2FoldChange) > 2 & -log10(padj) > 10),
            aes(x = log2FoldChange, y = -log10(pvalue) + 5, label = GENE_NAME), size = 3,) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("darkgray", "coral")) +
  labs(y = expression(-log[10]("p-value")),
       x = expression(log[2]("fold change")),
       title="Differential Gene Expression by Sex")
ggsave("VolcanoPlot.png",volcanoPlot)
