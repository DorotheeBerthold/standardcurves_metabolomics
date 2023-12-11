#####################################################
## Standard concentration curves & retention times ##
## Dorothée L. Berthold, ETH Zürich                ##
#####################################################


#Import required packages and functions
source("gutPackages.R")
gutPackages()

#Import analysis from Skyline together with Standard concentration curves & extract metadata from samplename
######################################################################################################################
#quant: long format df where each sample is matched with Area and retention time to each of the 128 metabolites
#conc: df for known Standard concentrations for 1:2 dilution series (Standard2) and 1:3 dilution series (Standard3)
#each standard is measured with both water (extraction buffer) and matrix (diluted cecal content) background

quant <- read.csv("skyline_DB017.csv")
quant$Molecule <- tolower(quant$Molecule)
conc <- read.csv("standard_conc.csv")

#Load standard library for 4min containing single-standard retention times both in positive and negative mode
library_neg <- read_excel("library_MSMM.xlsx", sheet = 1, col_names = F)
colnames(library_neg) <- c("Molecule", "Precursor formula", "rt_neg")
library_neg$Molecule <- tolower(library_neg$Molecule)
library_pos <- read_excel("library_MSMM.xlsx", sheet = 2, col_names = F)
colnames(library_pos) <- c("Molecule", "Precursor formula", "rt_pos")
library_pos$Molecule <- tolower(library_pos$Molecule)

#separate standard dilution, water, matrix & sample number into different columns
quant$string <- sub("\\.wiff2$", "", quant$File.Name)
quant$dilution <- str_left(quant$string, n= 9)
quant$background <- sub(".*_(water|matrix)_.*", "\\1", quant$string)
quant$mode <- sub(".*_(neg|pos)_.*", "\\1", quant$string)
quant$dilutionfactor <- as.numeric(sub(".*_(\\d+)$", "\\1", quant$string))
colnames(quant)[6:7] <- c("Area", "retention_time")
quant$Area <- as.numeric(quant$Area)
quant$retention_time <- as.numeric(quant$retention_time)

#logtransform Area
quant <- quant |> 
  mutate(log_area = log10(Area+1))

#filter out QC probes
qc <- c("QC_pos", "QC_pos2", "QC_neg")
quant <- quant |> 
  filter(!string %in% qc)

#Divide into positive and negative mode and match with log-transformed concentrations
######################################################################################################################
#logtransform concentrations
conc <- conc %>% 
  mutate(log_conc2 = log10(concentration2),
         log_conc3 = log10(concentration3))

#create vectors for concentrations
conc3 <- conc$log_conc3
conc3 <- na.omit(conc3)
conc2 <- conc$log_conc2

#separate quant based on standard-dilution factor 
#quant2: Standard dilution 1:2
#quant3: Standard dilution 1:3

quant2 <- quant %>% 
  filter(dilution == "Standard2")

quant3 <- quant %>% 
  filter(dilution == "Standard3")

#match with standard concentrations
quant3$log_conc <- conc3[match(quant3$dilutionfactor, seq_along(conc3))]
quant2$log_conc <- conc2[match(quant2$dilutionfactor, seq_along(conc2))]

#rowbind the two dataframes again
quant_conc <- bind_rows(quant2, quant3)

#generate grouping column for plotting
quant_conc$group <- paste0(quant_conc$background, "_", quant_conc$dilution)

#separate positive and negative mode measurements
quant_pos <- quant_conc |> 
  filter(mode == "pos")

quant_neg <- quant_conc |> 
  filter(mode == "neg")
