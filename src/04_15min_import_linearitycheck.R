#####################################################
## Standard concentration curves & retention times ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################



#Import analysis from Skyline together with Standard concentration curves & extract metadata from samplename
######################################################################################################################

quant_15 <- read.csv("skyline_DB017_15min.csv")
quant_15$Molecule <- tolower(quant_15$Molecule)
conc <- read.csv("standard_conc.csv")

#separate standard dilution, water, matrix & sample number into different columns
quant_15$string <- sub("\\.wiff2$", "", quant_15$File.Name)
quant_15$dilution <- str_left(quant_15$string, n= 9)
quant_15$background <- sub(".*_(water|matrix)_.*", "\\1", quant_15$string)
quant_15$mode <- sub(".*_(neg|pos)_.*", "\\1", quant_15$string)
quant_15$dilutionfactor <- as.numeric(sub(".*_(\\d+)$", "\\1", quant_15$string))
colnames(quant_15)[6:7] <- c("Area", "retention_time")
quant_15$Area <- as.numeric(quant_15$Area)
quant_15$retention_time <- as.numeric(quant_15$retention_time)

#logtransform Area
quant_15 <- quant_15 |> 
  mutate(log_area = log10(Area+1))

#filter out QC probes
qc <- c("QC_pos", "QC_pos2", "QC_neg")
quant_15 <- quant_15 |> 
  filter(!string %in% qc)

#match with standard concentrations
quant_15$log_conc <- conc3[match(quant_15$dilutionfactor, seq_along(conc3))]

#Check for linearity in the dilution changes
######################################################################################################################
quant_15_water <- quant_15 %>% 
  filter(background == "water")

dilution_factors <- 1:16
result_df_15 <- data.frame()

unique_combinations_15 <- quant_15_water |>  
  distinct(Molecule, mode)


dilution_factor15 <- quant_15_water %>%
  filter(dilutionfactor == 1)
dilution_factor15$dilution_change <- NA

for (i in 1:nrow(unique_combinations_15)) {
  molecule <- unique_combinations_15$Molecule[i]
  mode <- unique_combinations_15$mode[i]
  
  # Calculate dilution change for dilution factors 2 and above
  for (j in 2:length(dilution_factors)) {
    current_df <- dilution_factors[j]
    previous_df <- dilution_factors[j - 1]
    
    current_area <- quant_15_water$Area[quant_15_water$Molecule == molecule & quant_15_water$mode == mode & quant_15_water$dilutionfactor == current_df]
    previous_area <- quant_15_water$Area[quant_15_water$Molecule == molecule & quant_15_water$mode == mode & quant_15_water$dilutionfactor == previous_df]
    
    dilution_change <- current_area / previous_area
    dilution_change <- round(dilution_change, 4)  # Round to 4 decimal places
    
    result_row <- data.frame(Molecule = molecule, mode = mode, dilutionfactor = current_df, dilution_change = dilution_change)
    result_df_15 <- rbind(result_df_15, result_row)
  }
}

#only keep columns matching with result_df3
dilution_factor15 <- dilution_factor15[, names(result_df_15)]

# Combine the data frames
result_df_15 <- rbind(result_df_15, dilution_factor15)


#create different dfs after loop, attach to big df
quant_calc_15 <- left_join(quant_15_water, result_df_15, by = c("Molecule", "mode", "dilutionfactor"))

#Create a cut-off for the dilution_change based on histogram distribution
######################################################################################################################

#generate a subset without dilution factor 1 (as all NA) for histogram plotting
quant_calc_hist_15 <- quant_calc_15 |> 
  filter(dilution_change < 1.5)

ggplot(quant_calc_hist_15, aes(x = dilution_change, fill = mode)) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 30) +
  geom_vline(data = quant_calc_hist_15 %>% group_by(mode) %>% summarise(mean_value = mean(dilution_change)),
             aes(xintercept = mean_value, color = mode),
             linetype = "dashed", linewidth = 1) +
  labs(title = "Dilution Change Histogram",
       x = "Dilution Change",
       y = "Frequency") +
  theme_minimal()

filtered_df_15 <- quant_calc_15 %>%
  group_by(Molecule, mode) %>%
  filter(dilution_change != 0) %>% 
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

dilution_factor15 <- quant_calc_15 %>% 
  filter(dilutionfactor == 1)

filtered_df_15 <- rbind(filtered_df_15, dilution_factor15)
filtered_df_15$column <- "15min"
