#####################################################
## Standard concentration curves & retention times ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################



#Compare 4min & 15min mixed mode
######################################################################################################################

filtered_4min_15min <- rbind(filtered_df3, filtered_df_15)
filtered_4min_15min <- filtered_4min_15min %>% 
  select(-Precursor)

#generate grouping column for plotting
filtered_4min_15min$group <- paste0(filtered_4min_15min$column, "_", filtered_4min_15min$mode)

colour_mode <- c("4min_neg" = "lightblue", "15min_neg" = "navyblue", 
                 "4min_pos" = "lightgreen", "15min_pos" = "darkgreen")


#run regression model on filtered df depleted of column dilution_change

filtered_4min_15min_clean <- filtered_4min_15min %>% 
  select(-dilution_change) %>% 
  filter(log_area != 0)

lm_filtered_4min_15min <- filtered_4min_15min_clean %>%
  group_by(Molecule, mode, column) %>%
  summarize(Slope = coef(lm(log_area ~ log_conc))[2],
            Intercept = coef(lm(log_area ~ log_conc))[1],
            R_squared = summary(lm(log_area ~ log_conc))$r.squared)
lm_filtered_4min_15min$group <- paste0(lm_filtered_4min_15min$column, "_", lm_filtered_4min_15min$mode)

# Calculate the number of unique molecules
num_unique_molecules <- length(unique(filtered_4min_15min$Molecule))

# Set the desired number of columns and rows in each page
ncol_per_page <- 3
nrow_per_page <- 4

# Calculate the total number of pages
total_pages <- ceiling(num_unique_molecules / (ncol_per_page * nrow_per_page))

pdf("plots/concentration_curves_4min_15min.pdf")
for(i in 1:total_pages){
  print(ggplot(filtered_4min_15min, aes(log_conc, log_area, color = group)) +
          geom_point() +
          geom_abline(data=lm_filtered_4min_15min, aes(intercept = Intercept, slope = Slope, color = group)) +
          facet_wrap_paginate(~Molecule, scales = "free_y", ncol = 3, nrow= 4, page = i) +
          theme_classic() +
          scale_color_manual(values = colour_mode, name = "column & mode") +
          labs(title = "Comparison 4min & 15min mixed mode", x = "log10 concentration [uM]"))
}
dev.off()
