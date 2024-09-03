library(tidyverse)
library(palmerpenguins)
library(ggthemes)

glimpse(penguins)

# pick out only one columes

penguins[2,c("species","island")]
penguins[2,2]



###########################################################
# draw the data

ggplot(data=penguins) + 
  geom_point( mapping = aes(x= flipper_length_mm,
                                  y= body_mass_g, 
                                  colour = species,
                                  shape = species)) +
  geom_smooth(mapping = aes(x= flipper_length_mm,
                            y= body_mass_g),
              method ='lm')+
  scale_color_colorblind()+
  xlab("Flipper Lenghth (mm)")+
  ylab("Body Mass (g)")+
  ggtitle("Relationship between body mass\nand flipper length")

# save the graph
ggsave(filename = "~/qbb2024-answers/d1-afternoon/penguin-plot.pdf")



###########################################################
# Does bill length or depth depend on sex of the penguins?

male_penguin <- penguins %>% 
  filter(sex == "male")

ggplot(data = penguins,
       mapping = aes(x = bill_length_mm, fill = sex))+
  geom_histogram(position = "dodge")


ggplot(data = penguins,
       mapping = aes(x = bill_length_mm, fill = sex))+
  geom_histogram(position = "identity", alpha = 0.5)+
  facet_grid(.~sex)

ggplot(data = penguins %>% filter(!is.na(sex)),
       mapping = aes(x = bill_length_mm, fill = sex))+
  geom_histogram(position = "identity", alpha = 0.5)+
  facet_grid(sex~species)



###########################################################
# year ~ body mass relationship

ggplot(data = penguins,mapping = aes(x= factor(year),
                                     y= body_mass_g,
                                     fill = sex))+
  geom_boxplot() +#boxplot expects x to be a catagory
  facet_grid(island~species)

###########################################################
# density plot

ggplot(data = penguins %>% filter(!is.na(sex)),
       mapping = aes(x = bill_length_mm, color = sex))+
  geom_density()



