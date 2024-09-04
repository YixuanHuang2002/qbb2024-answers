library("tidyverse")

df <- read_tsv("~/Data/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

df <- df %>%
  mutate( SUBJECT=str_extract( SAMPID, "[^-]+-[^-]+" ), .before=1 )

df_sub <- df %>%
  group_by(SUBJECT) %>%
  summarise(my_counts = n() ) %>%
  arrange(desc(my_counts))



# K-562 has the most samples, second highest number of samples: GTEX-NPJ8
# GTEX-1JMI6, and GTEX-1PAR6 both have the least samples.

df_SMTSSD <- df %>%
  group_by(SMTSD) %>%
  summarise(b = n() ) %>%
  arrange(desc(b))

# Whole Blood has the most samples,second highest number of samples: Muscle - Skeletal 
# Kidney - Medulla has the least samples

#Filter for samples from this subject and save as a new object (e.g. df_npj8)
df_npj8 <- df %>%
  filter(SUBJECT == "GTEX-NPJ8")

df_npj8 %>%
  group_by(SMTSD) %>%
  summarise(b = n() ) %>%
  arrange(desc(b))
# Whole Blood has the most samples

blood <- df_npj8 %>%
  filter(SMTSD == "Whole Blood")

blood[,15:20]
# Difference
# they are sequenced and collected in different method

# Question 7
df_na <- df %>%
  filter( !is.na(SMATSSCR) )

df_0 <- df_na %>%
  group_by(SUBJECT) %>%
  summarise(mean(SMATSSCR)) %>%
  filter(`mean(SMATSSCR)`==0)
# 15 SUBJECTs have a mean SMATSSCR score of 0

df_mean <- df_na %>%
  group_by(SUBJECT) %>%
  summarise(mean(SMATSSCR))

range(df_mean$`mean(SMATSSCR)`) # the range of all the SMATSSCR mean is [0,2.571429]
mean(df_mean$`mean(SMATSSCR)`) # the mean of all the SMATSSCR mean is 0.9050125

# how to present the data?
# draw a dot plot to show the distributions 

