library(tidyverse)
library(ggthemes)

# read the file in step2_2
file <- read.delim("~/qbb2024-answers/9_20/snp_counts.txt") 

#Log transforms enrichment values
snpenrich <- file %>% mutate(logEnrichment = log2(Enrichment)) 

#Filter data by feature type
exons <- snpenrich %>% filter(Feature == "exons")
introns <- snpenrich %>% filter(Feature == "introns")
cCREs <- snpenrich %>% filter(Feature == "cCREs")
other <- snpenrich %>% filter(Feature == "other")

ggplot() +
  geom_line(data=exons,aes(MAF,logEnrichment,color="Exons")) + 
  geom_line(data=introns,aes(MAF,logEnrichment,color="Introns")) + 
  geom_line(data=cCREs,aes(MAF,logEnrichment,color="cCREs")) + 
  geom_line(data=other,aes(MAF,logEnrichment,color="Other")) + 
  labs(color="Legend",x="Minor Allele Frequency",y="SNP Enrichment",title="SNP Enrichment of Each Feature") + 
  ggsave(filename = "~/qbb2024-answers/9_20/snp_enrichments.pdf") 
