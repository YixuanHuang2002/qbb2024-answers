# Load necessary libraries
library(tidyverse)

# Load the data
data <- read_tsv("~/qbb2024-answers/day4/dicts_expr.tsv")

# Combine tissue and geneID into a single variable
data <- data %>%
  mutate(Tissue_Gene = paste0(Tissue, " ", GeneID)) %>%
  mutate(log_expression = log(Expr + 0.01))



# Load ggplot2 if not already loaded
library(ggplot2)

# Create the violin plot
ggplot(data, aes(x = reorder(Tissue_Gene, log_expression), y = log_expression)) +
  geom_violin(aes(fill = Tissue_Gene), alpha = 0.6) +
  geom_boxplot(width = 0.1, alpha = 0.2) +
  labs(
    title = "Expression Variability of Genes Across Tissues",
    x = "Tissue and Gene ID",
    y = "Log-Transformed Expression Level"
  ) +
  theme_minimal() +
  coord_flip()  # Flip coordinates to have categories on y-axis

# NOT surprised
# The plot shows that almost all the pancreas specific genes have high expression level. This suggests that gene expression in these tissue is more stable and tightly regulated.
# Testis, lung, and other tissue also have relatively high specific gene expression, but there is less consistency across the specific gene.
# Why?
# Pancreas is an important organ that have many different roles in endocrine 7 exocrine. The different part of pancreas have different roles like producing different enzymes.
# This may explain why Pancreas has higher expression variablity
# other organ like lung has relatively restricted functions(in breathing). 
# So their gene expression variablity is low due to their specific function.
# the gene expression plot fits with the organ's metabolic needs and physiological role.


