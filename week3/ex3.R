### 3.1 ###

library(ggplot2)

# open and read the file
allele_frequencies <- read.delim("~/qbb2024-answers/week3/AF.txt")

# Plot the histogram of allele frequencies using ggplot2
allele_frequency_spectrum <- ggplot(allele_frequencies, aes(x = Allele_Frequency)) +
  geom_histogram(bins = 11,  # Set the number of bins to 11
                 fill = "skyblue", 
                 color = "black",
                 aes(y = (..count..) / sum(..count..) * 100)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(title = "Allele Frequency Spectrum", 
       x = "Allele Frequency", 
       y = "Percent of Variants") +
  theme_minimal()
# save the plot
ggsave("~/qbb2024-answers/week3/allele_frequency_spectrum.png", plot = allele_frequency_spectrum, width = 8, height = 6, dpi = 300)

## Question 3.1
# DESCRIPTION: Few variants has very low or very high allele frequencies, most of them are in the middle. Generally, the percent of variants first increases and then decreases as allele frequency goes up.  
# Expected.
# Lethal or unfavorable alleles will be screened out by selection, so those variants will have low percentage in population
# Most of the mutations are recessive, and do not have destructive effect on the surviving of individuals, so the heterozygous (allele frequency ~0.5) will be the most in the population.
# Very low percentage of variants will have 100% allele frequency.
# normal distribution




### 3.2 ###

# open and read the file
read_depths <- read.delim("~/qbb2024-answers/week3/DP.txt")

# Plot the histogram of read depths using ggplot2
read_depth <- ggplot(read_depths, aes(x = Read_Depth)) +
  geom_histogram(bins = 21,  # Set bins = 21
                 fill = "lightgrey", 
                 color = "black", 
                 aes(y = ..count.. / sum(..count..) * 100)) +  
  # show %
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  xlim(0, 20) +  # Set x-axis limits
  labs(title = "the precentage of variants of different read depth",
       x = "Read Depth",
       y = "Percent of Variants") +
  theme_minimal()

ggsave("~/qbb2024-answers/week3/read_depth.png", plot = read_depth, width = 8, height = 6, dpi = 300)
## Question 3.2
# DESCRIPTION: 3-4 read depths is most common. There is a decline in the number of variants as the read depth increases after 4.
# Expected.
# it is reasonable to see that very low percent of variants has low read depth, most of theme have 3 or more reads covering them. This indicates sufficient coverage across most of the genome.
# The decline at higher depths is typical, as fewer regions are excessively sequenced. 
# poisson distribution