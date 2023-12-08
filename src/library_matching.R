#####################################################
## Standard concentration curves & retention times ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

#Matching up retention times with library of single-standards
######################################################################################################################
#use filtered_df clean to calculate median retention times (accounts for slight shifts in average) per Molecule & mode

retention_times_median <- filtered_df_clean %>% 
  group_by(Molecule, mode) %>% 
  summarise(rt_median = median(retention_time))

#create filter df to filter out molecules which are anyway undetectable (dilutionfactor 1 only)
filter_molecules <- detection_limit %>% 
  group_by(mode) %>% 
  filter(dilutionfactor == 1) %>% 
  distinct(Molecule, mode)

#anti-join the two dfs
filtered_rt <- anti_join(retention_times_median, filter_molecules, by = c("Molecule", "mode"))

#separate pos & negative for subsequent library matching
rt_neg <- filtered_rt %>% 
  filter(mode == "neg")

rt_pos <- filtered_rt %>% 
  filter(mode == "pos")

library_rt_neg <- inner_join(rt_neg, library_neg, by = "Molecule")
library_rt_pos <- inner_join(rt_pos, library_pos, by = "Molecule")

#calculate differences in rt 
library_rt_neg <- library_rt_neg %>% 
  mutate(rt_diff = abs(rt_median - rt_neg))

library_rt_pos <- library_rt_pos %>% 
  mutate(rt_diff = abs(rt_median - rt_pos))

#plot histogram of RT shifts

ggplot(library_rt_neg, aes(rt_diff)) + 
  geom_histogram() +
  labs(title = "Retention times shift negative mixed mode 4min", x = "Retention time difference") +
  theme_classic()

ggplot(library_rt_pos, aes(rt_diff)) + 
  geom_histogram() +
  labs(title = "Retention times shift positive mixed mode 4min", x = "Retention time difference") +
  theme_classic()
