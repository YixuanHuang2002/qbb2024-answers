library(tidyverse)
library(broom)

# Q1

dnm <- read_csv(file = "~/qbb2024-answers/day4/aau1043_dnm.csv")
ages <- read_csv(file = "~/qbb2024-answers/day4/aau1043_parental_age.csv")



dnm_summary <- dnm %>%
  group_by(Proband_id) %>%
  summarize(n_paternal_dnm = sum(Phase_combined == "father", na.rm = TRUE),
            n_maternal_dnm = sum(Phase_combined == "mother", na.rm = TRUE))


dnm_by_parental_age <- left_join(dnm_summary,ages,by = "Proband_id")


# Q2

# the count of maternal de novo mutations vs. maternal age
ggplot(data = dnm_by_parental_age,
       mapping = aes(x = Mother_age, y = n_maternal_dnm))+
  geom_point()+
  stat_smooth(method = 'lm')
# the count of paternal de novo mutations vs. paternal age
ggplot(data = dnm_by_parental_age,
       mapping = aes(x = Father_age, y = n_paternal_dnm))+
  geom_point()+
  stat_smooth(method = 'lm')

lm(data = dnm_by_parental_age,
   formula = n_maternal_dnm ~ 1 + Mother_age ) %>%
  summary()
# the size of the relationship: 0.37757
# When mother is one year older, the DNM will have 0.37757 more mutation.
# Match
# significant, p value is smaller than < 2e-16.This is a strong indicator that the observed differences are not due to random chance.

# 2.3
lm(data = dnm_by_parental_age,
   formula = n_paternal_dnm ~ 1 + Father_age) %>%
  summary()
# the size of the relationship: 1.35384
# When father is one year older, the DNM will have 1.35384 more mutation.
# Match
# significant, p value is smaller than < 2e-16. This is a strong indicator that the observed differences are not due to random chance.



# 2.4 predict the number of paternal DNMs for a proband with a father who was 50.5 years old 
10.32632 + 1.35384 * 50.5
# result: 78.69524


# 2.5  whether the number of paternally inherited DNMs match the number of maternally inherited DNMs. 
ggplot(data = dnm_by_parental_age)+
  geom_histogram(aes(x = n_maternal_dnm, y = ..density..), binwidth = 1, fill = "red", alpha = 0.5, position = "identity") +
  geom_histogram(aes(x = n_paternal_dnm, y = ..density..), binwidth = 1, fill = "blue", alpha = 0.5, position = "identity") +
  labs(title = "Distribution of Maternal and Paternal DNMs per Proband",
       x = "Number of DNMs",
       y = "Density") +
  theme_minimal()
# Overall, male will have more inheritable DNMs than female 
# But DNM happens more frequently in mother than in father

# 2.6 there is a significant difference between the number of maternally vs. paternally inherited DNMs per proband
t.test(dnm_by_parental_age$Father_age,dnm_by_parental_age$Mother_age,paired = TRUE)
# Results:
# t = 15.821, df = 395, p-value < 2.2e-16
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:  3.131405 4.020110
# mean difference: 3.575758 

# What statistical test did you choose? Why?
# t-test.
# Why paired?
# Paired samples occur when each observation in one group is directly related to an observation in another group. 
# In this context, each proband has both maternal and paternal DNMs, making the samples paired.
# This pairing means that the observations are not independent of each other. Thus a paired t-test would fit the model.
# Why two-side?
# A two-sided test assesses whether there is a significant difference in either direction.
# We don't have a hypothesis about which DNMs to be more prevalent.
# We are interested in determining whether there is any significant difference between maternal and paternal DNMs, regardless of direction. A two-sided t-test would be fine.
# In conclusion, a paired, two-sided t-test would fit the model best.

# Yes.According to the results, p-value < 2.2e-16.

