#####################################################
## Standard concentration curves & retention times ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

#Fit linear models for standard lines and R squares for determination of matrix effects
######################################################################################################################
#separate positive and negative mode measurements
quant_pos <- quant_conc |> 
  filter(mode == "pos")

quant_neg <- quant_conc |> 
  filter(mode == "neg")

#Fit different linear model for each mode
standard_curve_neg <- quant_neg %>%
  filter(log_area != 0) %>% 
  group_by(Molecule, background, dilution) %>%
  summarize(Slope = coef(lm(log_area ~ log_conc))[2],
            Intercept = coef(lm(log_area ~ log_conc))[1],
            R_squared = summary(lm(log_area ~ log_conc))$r.squared)
standard_curve_neg$R_squared <- round(standard_curve_neg$R_squared, 3)

standard_curve_pos <- quant_pos %>%
  filter(log_area != 0) %>% 
  group_by(Molecule, background, dilution) %>%
  summarize(Slope = coef(lm(log_area ~ log_conc))[2],
            Intercept = coef(lm(log_area ~ log_conc))[1],
            R_squared = summary(lm(log_area ~ log_conc))$r.squared)
standard_curve_pos$R_squared <- round(standard_curve_pos$R_squared, 3)

standard_curve_neg$group <- paste0(standard_curve_neg$background, "_", standard_curve_neg$dilution)
standard_curve_pos$group <- paste0(standard_curve_pos$background, "_", standard_curve_pos$dilution)

#Concentration curves plots
######################################################################################################################
#define colour vector for plotting

colour_groups <- c("matrix_Standard2" = "darkgreen", "matrix_Standard3" = "springgreen", 
                   "water_Standard2" = "navyblue", "water_Standard3" = "lightblue")


#plot different metabolites, separated by background, faceted by metabolite

pdf("plots/concentration_curves_metabolites_negativemode_4min.pdf")
for(i in 1:15){
  print(ggplot(quant_neg, aes(log_conc, log_area, color = group)) +
          geom_point() +
          geom_abline(data=standard_curve_neg, aes(intercept = Intercept, slope = Slope, color = group)) +
          facet_wrap_paginate(~Molecule, scales = "free_y", ncol = 3, nrow= 3, page = i) +
          theme_classic() +
          scale_color_manual(values = colour_groups, name = "Background & dilution", labels = c("matrix_1:2", "matrix_1:3", "water_1:2", "water_1:3")) +
          labs(title = "Negative mixed mode 4min", x = "log10 concentration [uM]", legend = "background + dilution"))
}
dev.off()

pdf("plots/concentration_curves_metabolites_positivemode_4min.pdf")
for(i in 1:15){
  print(ggplot(quant_pos, aes(log_conc, log_area, color = group)) +
          geom_point() +
          geom_abline(data=standard_curve_neg, aes(intercept = Intercept, slope = Slope, color = group)) +
          facet_wrap_paginate(~Molecule, scales = "free_y", ncol = 3, nrow= 3, page = i) +
          theme_classic() +
          scale_color_manual(values = colour_groups, name = "Background & dilution", labels = c("matrix_1:2", "matrix_1:3", "water_1:2", "water_1:3")) +
          labs(title = "Positive mixed mode 4min", x = "log10 concentration [uM]", legend = "background + dilution"))
}
dev.off()

#Retention time plots
######################################################################################################################

quant_neg <- quant_neg %>% 
  mutate(conc = 10^(log_conc))


quant_pos <- quant_pos %>% 
  mutate(conc = 10^(log_conc))


pdf("plots/retentiontime_negativemode_4min.pdf")
for(i in 1:15){
  print(ggplot(quant_neg, aes(conc, retention_time, color = group)) +
          geom_line() +
          facet_wrap_paginate(~Molecule, scales = "free_y", ncol = 3, nrow= 3, page = i) +
          theme_classic() +
          scale_color_manual(values = colour_groups, name = "Background & dilution", labels = c("matrix_1:2", "matrix_1:3", "water_1:2", "water_1:3")) +
          labs(title = "Retention times [negative mixed mode 4min]", x = "Concentration [uM]", y = "Retention time"))
}
dev.off()

pdf("plots/retentiontime_positiveemode_4min.pdf")
for(i in 1:15){
  print(ggplot(quant_pos, aes(conc, retention_time, color = group)) +
          geom_line() +
          facet_wrap_paginate(~Molecule, scales = "free_y", ncol = 3, nrow= 3, page = i) +
          theme_classic() +
          scale_color_manual(values = colour_groups, name = "Background & dilution", labels = c("matrix_1:2", "matrix_1:3", "water_1:2", "water_1:3")) +
          labs(title = "Retention times [positive mixed mode 4min]", x = "Concentration [uM]", y = "Retention time"))
}
dev.off()

