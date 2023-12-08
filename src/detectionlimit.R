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

write_csv(detection_limit, "results/detectionlimit_4min.csv")
