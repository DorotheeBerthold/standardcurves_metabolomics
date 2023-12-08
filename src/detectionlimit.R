#####################################################
## Standard concentration curves & retention times ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

#Detection limit and residual calculations
######################################################################################################################

detection_limit <- filtered_df_clean %>%
  group_by(Molecule, mode, dilution) %>%
  mutate(conc_nM = 10^(log_conc)) %>% 
  slice_max(order_by = dilutionfactor) %>% 
  ungroup() %>% 
  select(Molecule, mode, dilution, dilutionfactor, conc_nM)

detection_limit <- filtered_df_clean %>%
  group_by(Molecule, mode, dilution) %>%
  mutate(conc_nM = 10^(log_conc)) %>% 
  slice_max(order_by = dilutionfactor) %>% 
  ungroup() %>% 
  select(Molecule, mode, dilution, dilutionfactor, conc_nM)

#add R-squared to detection limit
detection_r2 <- inner_join(detection_limit, lm_filtered, by = c("Molecule", "mode", "dilution"))

detection_r2 <- detection_r2 %>% 
  filter(dilutionfactor != 1)

write_csv(detection_r2, "results/detectionlimit_4min.csv")

ggplot(detection_r2, aes(R_squared, fill = dilution)) + 
  geom_histogram(alpha = 0.5, binwidth = 0.002) +
  labs(title = "R-squared of standard line fit", x = "R_squared") +
  theme_classic()


#Detection limit 15min and 4min combined (Standard 1:3 dilution)
######################################################################################################################

detection_limit_4min_15min <- filtered_4min_15min_clean %>%
  group_by(Molecule, mode, column) %>%
  mutate(conc_nM = 10^(log_conc)) %>% 
  slice_max(order_by = dilutionfactor) %>% 
  ungroup() %>% 
  select(Molecule, mode, column, dilutionfactor, conc_nM)

#add R-squared to detection limit
detection_r2_4min_15min <- inner_join(detection_limit_4min_15min, lm_filtered_4min_15min, by = c("Molecule", "mode", "column"))

detection_r2_4min_15min <- detection_r2_4min_15min %>% 
  filter(dilutionfactor != 1)

write_csv(detection_r2_4min_15min, "results/detectionlimit_4min_15min.csv")


ggplot(detection_r2_4min_15min, aes(R_squared, fill = group)) + 
  geom_histogram(alpha = 0.5, binwidth = 0.002) +
  labs(title = "R-squared of standard line fit", x = "R_squared") +
  theme_classic()


