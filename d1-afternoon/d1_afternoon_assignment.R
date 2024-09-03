library(tidyverse)

###########################################################
# question 1 : read

df <- read_delim("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

# Question 2 
glimpse(df)
df[1,]

# Question 3
rna_seq <- df %>% 
  filter(SMGEBTCHT == "TruSeq.v1")

# Question 4
ggplot(data = rna_seq,
       mapping = aes(x = SMTSD))+
  geom_bar()+
  ylab("Number of Samples")+
  ggtitle("Number of Samples of Each Tissue")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Question 5
# Refined Histogram
ggplot(rna_seq, aes(x = SMRIN)) +
  geom_histogram(binwidth = 0.5, color = "black", alpha = 0.5) +
  xlab("RNA Integrity Number")+
  ylab("Freqency Number")+
  ggtitle("Distribution of RNA Integrity Numbers") 
# The shape is unimodal

# Question 6
ggplot(rna_seq, aes(x = SMTSD,
                    y = SMRIN)) +
  geom_boxplot() +
  xlab("Tissue Type")+
  ylab("RNA Integrity Number")+
  ggtitle("Relationship between Tissue Type and RNA Integrity Number") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# There are differences.
# Three tissue have specially high RNA integrity number:
# Cell- Cultured fibroblaste, Cells-EBV-transformed lymphocytes, and Cells-CML
# Those cells could have RNAs least degraded, meaning those cells have to produce many proteins, which means they replicate themselves very quickly. 


# Question 7
ggplot(rna_seq, aes( x = SMTSD, y = SMGNSDTC)) +
  geom_boxplot() +
  xlab("Tissue Type")+
  ylab("Number of genes detected ")+
  ggtitle("Relationship between Tissue Type and Number of genes detected") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# There are differences.
# The testis expresses the largest number of genes of any mammalian organ. 
# DNA sequence integrity in the male germline by correcting DNA damage through a mechanism we term transcriptional scanning.
# So Testis need to replicate a lot more genes than other tissue to make sure that they can correct almost all the DNA damage throughout the genome. 


#Q8
ggplot(rna_seq, mapping = aes(x = SMTSISCH,
                    y = SMRIN
                    )) +
  geom_point(size = 0.5 , alpha = 0.5)+
  geom_smooth(method ='lm', color = "blue",size = 0.5)+
  facet_wrap("SMTSD")+
  xlab("ischemic time")+
  ylab("RIN ")+
  ggtitle("Relationship between ischemic time and RIN")
# What relationships do you notice? most of them are negatively correlated in most tissue types
# Does the relationship depend on tissue? No

# Q9

ggplot(data=rna_seq)+
  geom_point(size = 0.5,alpha = 0.5, mapping=aes(x = SMTSISCH,
                         y = SMRIN, color = SMATSSCR
                         )) +
  facet_wrap(vars(SMTSD))+
  geom_smooth(methods='lm', mapping=aes(x = SMTSISCH, y = SMRIN))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("ischemic time")+
  ylab("RIN ")+
  ggtitle("Relationship between ischemic time and RIN")

# What relationships do you notice? a negative correlation between ischemic time and RIN, but a positive correlation between ischemic time and SMATSSCR
# Does the relationship depend on tissue? YES






  