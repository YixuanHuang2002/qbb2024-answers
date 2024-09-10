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
