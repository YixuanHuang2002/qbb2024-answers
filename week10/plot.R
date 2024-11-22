# Load required libraries
library(ggplot2)
library(readr)

# Load the data
data <- read_tsv("~/qbb2024-answers/week10/nucleus_signals.csv")

# Violin plot for Nascent RNA Signal
ggplot(data, aes(x = Gene, y = `Mean Nascent RNA Signal`, fill = Gene)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Distribution of Nascent RNA Signal Across Genes",
       x = "Gene",
       y = "Mean Nascent RNA Signal") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set1")

# Violin plot for PCNA Signal
ggplot(data, aes(x = Gene, y = `Mean PCNA Signal`, fill = Gene)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Distribution of PCNA Signal Across Genes",
       x = "Gene",
       y = "Mean PCNA Signal") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

# Violin plot for Log2 Ratio
ggplot(data, aes(x = Gene, y = `Log2 Ratio`, fill = Gene)) +
  geom_violin(trim = FALSE) +
  theme_minimal() +
  labs(title = "Distribution of Log2 Ratio (Nascent RNA / PCNA) Across Genes",
       x = "Gene",
       y = "Log2 Ratio") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")