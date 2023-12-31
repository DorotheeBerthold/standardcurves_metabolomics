#####################################################
## Standard concentration curves & retention times ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################

#Check for linearity in the dilution changes
######################################################################################################################
#for this, we are only considering the water dilutions
#also, we run a separate loop for Standard 1:2 & Standard 1:3 as there were troubles with repeating values in the
#nested loop

quant_conc_water <- quant_conc %>% 
  filter(background == "water", 
         dilution == "Standard2")

#result_df2 is the output of 1:2 dilution loop
result_df2 <- data.frame()

# Define the dilution factors
dilution_factors <- 1:24


# Loop through each combination of Molecule and mode and dilutionfactor
unique_combinations <- quant_conc_water |>  
  distinct(Molecule, mode)


# Calculate the dilution change

dilution_factor1 <- quant_conc_water %>%
  filter(dilutionfactor == 1)
dilution_factor1$dilution_change <- NA

for (i in 1:nrow(unique_combinations)) {
  molecule <- unique_combinations$Molecule[i]
  mode <- unique_combinations$mode[i]
  
  # Calculate dilution change for dilution factors 2 and above
  for (j in 2:length(dilution_factors)) {
    current_df <- dilution_factors[j]
    previous_df <- dilution_factors[j - 1]
    
    current_area <- quant_conc_water$Area[quant_conc_water$Molecule == molecule & quant_conc_water$mode == mode & quant_conc_water$dilutionfactor == current_df]
    previous_area <- quant_conc_water$Area[quant_conc_water$Molecule == molecule & quant_conc_water$mode == mode & quant_conc_water$dilutionfactor == previous_df]
    
    dilution_change <- current_area / previous_area
    dilution_change <- round(dilution_change, 4)  # Round to 4 decimal places
    
    result_row <- data.frame(Molecule = molecule, mode = mode, dilutionfactor = current_df, dilution_change = dilution_change)
    result_df2 <- rbind(result_df2, result_row)
  }
}

#only keep columns matching with result_df3
dilution_factor1 <- dilution_factor1[, names(result_df2)]

# Combine the data frames
result_df2 <- rbind(result_df2, dilution_factor1)
result_df2$dilution <- "Standard2"



#repeat for 1:3 Standard dilution --> result_df3

quant_conc_water <- quant_conc %>% 
  filter(background == "water", 
         dilution == "Standard3")

dilution_factors <- 1:16
result_df3 <- data.frame()

unique_combinations <- quant_conc_water |>  
  distinct(Molecule, mode)


dilution_factor1 <- quant_conc_water %>%
  filter(dilutionfactor == 1)
dilution_factor1$dilution_change <- NA

for (i in 1:nrow(unique_combinations)) {
  molecule <- unique_combinations$Molecule[i]
  mode <- unique_combinations$mode[i]
  
  # Calculate dilution change for dilution factors 2 and above
  for (j in 2:length(dilution_factors)) {
    current_df <- dilution_factors[j]
    previous_df <- dilution_factors[j - 1]
    
    current_area <- quant_conc_water$Area[quant_conc_water$Molecule == molecule & quant_conc_water$mode == mode & quant_conc_water$dilutionfactor == current_df]
    previous_area <- quant_conc_water$Area[quant_conc_water$Molecule == molecule & quant_conc_water$mode == mode & quant_conc_water$dilutionfactor == previous_df]
    
    dilution_change <- current_area / previous_area
    dilution_change <- round(dilution_change, 4)  # Round to 4 decimal places
    
    result_row <- data.frame(Molecule = molecule, mode = mode, dilutionfactor = current_df, dilution_change = dilution_change)
    result_df3 <- rbind(result_df3, result_row)
  }
}

#only keep columns matching with result_df3
dilution_factor1 <- dilution_factor1[, names(result_df3)]

# Combine the data frames
result_df3 <- rbind(result_df3, dilution_factor1)
result_df3$dilution <- "Standard3"

#bind result_df2 & result_df3 together before joining with quant_conc
result_joined <- rbind(result_df2, result_df3)

quant_conc_water <- quant_conc %>% 
  filter(background == "water")


#create different dfs after loop, attach to big df
quant_calc <- left_join(quant_conc_water, result_joined, by = c("Molecule", "mode", "dilution", "dilutionfactor"))

#Create a cut-off for the dilution_change based on histogram distribution
######################################################################################################################

#generate a subset without dilution factor 1 (as all NA) for hsitogram plotting
quant_calc_hist <- quant_calc |> 
  filter(dilution_change < 1.5)

ggplot(quant_calc_hist, aes(x = dilution_change, fill = dilution)) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 30) +
  geom_vline(data = quant_calc_hist %>% group_by(dilution) %>% summarise(mean_value = mean(dilution_change)),
             aes(xintercept = mean_value, color = dilution),
             linetype = "dashed", linewidth = 1) +
  labs(title = "Dilution Change Histogram",
       x = "Dilution Change",
       y = "Frequency") +
  theme_minimal()


#Generate separate filtered_df based on 1:2 and 1:3 dilution
#1:2 "perfect" dilution_change would be 0.5
filtered_df2 <- quant_calc %>%
  group_by(Molecule, mode, dilution) %>%
  filter(dilution_change != 0) %>% 
  filter(dilution == "Standard2") %>%
  do({
    keep_rows <- numeric()
    for (i in 2:24) {
      if (all(between(.$dilution_change[.$dilutionfactor == i], 0.5 - 0.25, 0.5 + 0.25))) {
        keep_rows <- c(keep_rows, i)
      } else {
        break  # Terminate the loop if the condition is not met
      }
    }
    filter(., dilutionfactor %in% keep_rows)
  })

dilution_factor1 <- quant_calc %>% 
  filter(dilution == "Standard2",
         dilutionfactor == 1)

filtered_df2 <- rbind(filtered_df2, dilution_factor1)
#1:3 "perfect" dilution_change would be 0.33
filtered_df3 <- quant_calc %>%
  group_by(Molecule, mode, dilution) %>%
  filter(dilution_change != 0) %>% 
  filter(dilution == "Standard3") %>%
  do({
    keep_rows <- numeric()
    for (i in 2:16) {
      if (all(between(.$dilution_change[.$dilutionfactor == i], 0.33 - 0.25, 0.33 + 0.25))) {
        keep_rows <- c(keep_rows, i)
      } else {
        break  # Terminate the loop if the condition is not met
      }
    }
    filter(., dilutionfactor %in% keep_rows)
  })

dilution_factor1 <- quant_calc %>% 
  filter(dilution == "Standard3",
         dilutionfactor == 1)

filtered_df3 <- rbind(filtered_df3, dilution_factor1)

#bind them together
filtered_df <- rbind(filtered_df2, filtered_df3)
filtered_df$group <- paste0(filtered_df$dilution, "_", filtered_df$mode)

colour_mode <- c("Standard2_neg" = "pink", "Standard3_neg" = "plum", 
                   "Standard2_pos" = "brown4", "Standard3_pos" = "brown2")


#run regression model on filtered df depleted of column dilution_change

filtered_df_clean <- filtered_df %>% 
  select(-dilution_change) %>% 
  filter(log_area != 0)

lm_filtered <- filtered_df_clean %>%
  group_by(Molecule, mode, dilution) %>%
  summarize(Slope = coef(lm(log_area ~ log_conc))[2],
            Intercept = coef(lm(log_area ~ log_conc))[1],
            R_squared = summary(lm(log_area ~ log_conc))$r.squared)
lm_filtered$group <- paste0(lm_filtered$dilution, "_", lm_filtered$mode)

# Calculate the number of unique molecules
num_unique_molecules <- length(unique(filtered_df$Molecule))

# Set the desired number of columns and rows in each page
ncol_per_page <- 3
nrow_per_page <- 4

# Calculate the total number of pages
total_pages <- ceiling(num_unique_molecules / (ncol_per_page * nrow_per_page))

pdf("plots/concentration_curves_corrected_4min.pdf")
for(i in 1:total_pages){
  print(ggplot(filtered_df_clean, aes(log_conc, log_area, color = group)) +
          geom_point() +
          geom_abline(data=lm_filtered, aes(intercept = Intercept, slope = Slope, color = group)) +
          facet_wrap_paginate(~Molecule, scales = "free_y", ncol = 3, nrow= 4, page = i) +
          theme_classic() +
          scale_color_manual(values = colour_mode, name = "Dilution & mode", labels = c("neg_1:2", "pos_1:2", "neg_1:3", "pos_1:3")) +
          labs(title = "Mixed mode 4min", x = "log10 concentration [uM]"))
}
dev.off()
