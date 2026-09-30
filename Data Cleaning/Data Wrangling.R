# Kelley Sinning 9/15/2026
# Cleaning up all the manjor datasheets so they can be used in an analysis together!

# trying something

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

setwd("~/Library/CloudStorage/OneDrive-TheUniversityofMontana/Data/BVR/Data Cleaning")

diets <- read.csv("Diet_Data.csv")
bugs <- read.csv("BMI_SurberData.csv")
fish <- read.csv("Meta_Field_Data.csv")


# First, we will start with cleaning up the diets data sheet
diets <- diets %>%
  select(-Initials, -Date_entered, -IDer, -ID_date, -Occasion, -Diet_observation_number,
         -Measurement_mm,-Extra_counts, -Empty_case_mm, -PROOFED, -USE_SAMPLE, -Total_Sampled, -X,
         -NOTES) %>%
  mutate(across(measurement_1:measurement_10, ~ as.numeric(as.character(.x)))) %>%
  mutate(Length_mean_mm = rowMeans(select(., measurement_1:measurement_10), na.rm = TRUE)) %>%
  select(-(measurement_1:measurement_10)) %>%
  rename(Abundance = Total.Measured)



# Renaming  typos

diets <- diets %>%   
  mutate(Order = recode_values(Order,
                            c("Ephemeropterta", "Ephemerroptera", "Emphemeroptera", "Ehpemeroptera", "Ephemerillidae",
                              "Ephemeropterra", "Eohemeroptera", "Ephemerptera", "Epehemeroptera", "Ephemeroptera ") ~ "Ephemeroptera",
                            c("Dptera", "Dipteraa", "diptera", "Diptera ", "DIptera") ~ "Diptera",
                            c("Trichooptera", "Tricchoptera", "Trichopteraa", "Trichptera", "Trichoptera (empty)",
                              "Trichiptera", "Tirchoptera", " Trichoptera", "trichoptera", "Trichoptera ") ~ "Trichoptera",
                            c("Hemioptera") ~ "Hemiptera",
                            c("Plecotera") ~ "Plecoptera",
                            c("Miscellllanous", "Miscelllanous", "miscellanous", "Miscellanous", "Miscellanous",
                              "Misc.", "Miscellanous ") ~ "Miscellaneous",
                            c("Terrestrial ") ~ "Terrestrial",
                            c("Larval fish ", "Trout fry", "Fish Fry", "Trout", "Fry") ~ "Trout Fry",
                            c("empty", "EMPTY", "Empty ", "MT") ~ "Empty",
                            c("isopoda", "Isopoda") ~ "Asellidae",
                            c("insecta") ~ "Insecta",
                            c("Coleoptera ") ~ "Coleoptera",
                            c("Aquatic mite", "Hydrachidae") ~ "Hydrachnidae",
                            c("Trombiditonnes", "Trombidiformes") ~ "Ant",
                            c("Hirudinea") ~ "Leech",
                            c("gastropoda") ~ "Gastropoda",
                            c("Polychaeta") ~ "Oligochaeta",
                            default = Order))
                            
unique(diets$Order) # check this for more fixes as data is added    


diets <- diets %>%
  mutate(Family = recode_values(Family,
                                c("Baaetidae", "Baetidaae", "Baetidae ", "baetidae") ~ "Baetidae",
                                c("Lepdostomatidae", "Lepidostommatidae", "Lepidostomatdae", "Lepidstomatidae",
                                  "Lepdiostomatidae", "Lepidostomatidae (empty)", "lepidostomatidae") ~ "Lepidostomatidae",
                                c("Brachycetridae", "Brachycentratidae", "Bracychentratidae", "Brachycentridae ") ~ "Brachycentridae",
                                c("Rhyacolphilidae", "Rhycophilidae") ~ "Rhyacophilidae",
                                c("Perlodide", "Perlodd") ~ "Perlodidae",
                                c("Chiroonomidae", "Chronomidae", "Chironimidae", "Chiironomidae", "Chironomidaae",
                                  "Chironomidaee", "Chironmidae", "Chironomodae", "Chironomidae ") ~ "Chironomidae",
                                c("Ephemerllidae", "Ephemerillidae", "Epemerelllidae", "Ephemeriillidae",
                                  "Ephemerilllidae", "Ephermerellidae") ~ "Ephemerellidae",
                                c("Glossostomatidae", "Glossostoomatidae") ~ "Glossosomatidae",
                                c("Hydropsychiae", "Hydropyschidae", "Hydropsychiodae", "Hydropschidae", "Hydropsychidae ") ~ "Hydropsychidae",
                                c("Simulidae", "Simullidae", "Simuidae", "Simulliidae", "Simluidae", "SImulidae", "simuliidae") ~ "Simuliidae",
                                c("Tipuliidae", "Tipullidae", "Tiplulidae") ~ "Tipulidae",
                                c("Pelidae", "perlidae", "Perlidae ") ~ "Perlidae",
                                c("Spider", "(spider)", "Sipder", "spider") ~ "Araneae",
                                c("Anthericidae", "Athericidae ") ~ "Athericidae",
                                "Aselllidae" ~ "Asellidae",
                                c("Formiicidae", "Ant", "ant") ~ "Formicidae",
                                c("Hydrophilidae?", "Hyfrophilidae") ~ "Hydrophilidae",
                                c("Leptophlebidae", "Leptophlebididae", "Leptophlebidiae", "Leptoplebiidae") ~ "Leptophlebiidae",
                                "Fly" ~ "Diptera",
                                c("Snail", "Snail ") ~ "Gastropoda",
                                "Hydrachinidae" ~ "Hydrachnidae",
                                "Stratomyidae" ~ "Stratiomyidae",
                                "Empidiae" ~ "Empididae",
                                c("Trout", "Brown trout fry", "Fry", "Fish fry", "Sculpin", "Sculpin fry", "Sculpin Fry", "Trout fry",
                                  "Trout fry ", "FIsh fry", "Fish fry", "Fish fry", "Fish fry ", "Fish Fry", "FIsh Fry") ~ "Fry",
                                c("", "N/A", "-", "Nothing measurable", "(gut)", "nothing measurable") ~ NA_character_,
                                "Chloroperlidae " ~ "Chloroperlidae",
                                "Ophiogomphus" ~ "Gomphidae",
                                "Homoptera" ~ "Hemiptera",
                                "Nemouridae " ~ "Nemouridae",
                                "Tabanidae " ~ "Tabanidae",
                                default = Family))

unique(diets$Family) # check this for more fixes as data is added    

# For biomass
ADDING IN LENGTH-MASS EQUATIONS
SECPROD <- SECPROD %>%
  mutate(
    Biomass.mg = case_when(
      Genus == "Acentrella" ~ (0.00962 * (Length ^ 2.75)) * Abundance,
      Genus == "Acerpenna" ~ (0.0076 * (Length ^ 2.691)) * Abundance,
      Genus == "Acroneuria" ~ (0.0019 * (Length ^ 3.232)) * Abundance,
      Genus == "Allocapnia" ~ (0.004 * (Length ^ 2.487)) * Abundance,
      Genus == "Allognasta" ~ (0.0032 * (Length ^ 2.61)) * Abundance,
      Genus == "Alloperla" ~ (0.0062 * (Length ^ 2.724)) * Abundance,
      Genus == "Ameletus" ~ (0.0077 * (Length ^ 2.588)) * Abundance,
      Genus == "Amphinemura" ~ (0.004 * (Length ^ 2.975)) * Abundance,
      Genus == "Antocha" ~ (0.0041 * (Length ^ 2.4439)) * Abundance,
      Genus == "Atherix" ~ (0.0038 * (Length ^ 2.586)) * Abundance,
      Genus == "Attenella" ~ (0.00928 * (Length ^ 2.9)) * Abundance,
      Genus == "Baetis" ~ (0.0076 * (Length ^ 2.691)) * Abundance,
      Genus == "Baetidae" ~ (0.0076 * (Length ^ 2.691)) * Abundance,
      Genus == "Baetisca" ~ (0.0116 * (Length ^ 2.905)) * Abundance,
      Genus == "Boyeria" ~ (0.0082 * (Length ^ 2.813)) * Abundance,
      Genus == "Calopteryx" ~ (0.005 * (Length ^ 2.742)) * Abundance,
      Genus == "Cernotina" ~ (0.0071 * (Length ^ 2.531)) * Abundance, # Polycentropodidae family
      Genus == "Chauliodes" ~ (0.0062 * (Length ^ 2.691)) * Abundance,
      Genus == "Chelifera" ~ (0.004 * (Length ^ 2.655)) * Abundance,
      Genus == "Cheumatopsyche" ~ (0.0045 * (Length ^ 2.721)) * Abundance,
      Genus == "Chimarra" ~ (0.0044 * (Length ^ 2.652)) * Abundance,
      Genus == "Chironomini" ~ (0.0007 * (Length ^ 2.952)) * Abundance,
      Genus == "Collembola" ~ (0.0024 * (Length ^ 3.676)) * Abundance,
      Genus == "Cordulegaster" ~ (0.0067 * (Length ^ 2.782)) * Abundance,
      Genus == "Cyrnellus" ~ (0.0071 * (Length ^ 2.531)) * Abundance,
      Genus == "Dicranota" ~ (0.0027 * (Length ^ 2.637)) * Abundance,
      Genus == "Diplectrona" ~ (0.0049 * (Length ^ 2.62)) * Abundance,
      Genus == "Discocerina" ~ (0.00033 * (Length ^ 3.55)) * Abundance,
      Genus == "Dixa" ~ (0.0433 * (Length ^ 0)) * Abundance,
      Genus == "Dixella" ~ (0.0433 * (Length ^ 0)) * Abundance,
      Genus == "Dolophilodes" ~ (0.00408 * (Length ^ 2.82)) * Abundance,
      Genus == "Ectopria" ~ (0.0164 * (Length ^ 2.929)) * Abundance,
      Genus == "Eloeophila" ~ (0.0014 * (Length ^ 2.667)) * Abundance,
      Genus == "Ephemera" ~ (0.0021 * (Length ^ 2.737)) * Abundance,
      Genus == "Epeorus" ~ (0.0121 * (Length ^ 2.667)) * Abundance,
      Genus == "Eriopterini" ~ (0.0016 * (Length ^ 2.914)) * Abundance, # Limoniidae family
      Genus == "Eurylophella" ~ (0.008 * (Length ^ 2.663)) * Abundance,
      Genus == "Gerris" ~ (0.015 * (Length ^ 2.596)) * Abundance,
      Genus == "Glossosoma" ~ (0.0092 * (Length ^ 2.888)) * Abundance,
      Genus == "Goera" ~ (0.00156 * (Length ^ 2.75)) * Abundance,
      Genus == "Gomphus" ~ (00.044 * (Length ^ 3.124)) * Abundance,
      Genus == "Gyrinus" ~ (0.0531 * (Length ^ 2.586)) * Abundance, # Dineutes sp. from Benke
      Genus == "Hagenella" ~ (0.0054 * (Length ^ 2.811)) * Abundance, # Ptilostomis from Benke
      Genus == "Hemiptera" ~ (0.00836 * (Length ^ 3.075)) * Abundance,
      Genus == "Hetaerina" ~ (0.005 * (Length ^ 2.742)) * Abundance,
      Genus == "Helichus" ~ (0.0011 * (Length ^ 3.1)) * Abundance,
      Genus == "Hexatoma" ~ (0.0042 * (Length ^ 2.596)) * Abundance,
      Genus == "Hydatophylax" ~ (0.0049 * (Length ^ 2.85)) * Abundance,
      Genus == "Hydropsyche" ~ (0.0051 * (Length ^ 2.824)) * Abundance,
      Genus == "Isonychia" ~ (0.0031 * (Length ^ 3.167)) * Abundance,
      Genus == "Isoperla" ~ (0.01 * (Length ^ 2.658)) * Abundance,
      Genus == "Langessa" ~ (0.00271 * (Length ^ 2.959)) * Abundance, # Lepidoptera from Greg's sheet 
      Genus == "Lanthus" ~ (0.0097 * (Length ^ 2.895)) * Abundance,
      Genus == "Lepidostoma" ~ (0.0079 * (Length ^ 2.649)) * Abundance,
      Genus == "Leuctra" ~ (0.003 * (Length ^ 2.761)) * Abundance,
      Genus == "Limnephilidae" ~ (0.0049 * (Length ^ 2.85)) * Abundance,
      Genus == "Limnophila" ~ (0.0014 * (Length ^ 2.667)) * Abundance,
      Genus == "Limoniidae" ~ (0.0016 * (Length ^ 2.914)) * Abundance,
      Genus == "Lypodiversa" ~ (0.0039 * (Length ^ 2.873)) * Abundance, # Polycentropodidae family
      Genus == "Micrasema" ~ (0.0181 * (Length ^ 2.410)) * Abundance,
      Genus == "Microvelia" ~ (0.0083 * (Length ^ 2.777)) * Abundance,
      Genus == "Molophilus" ~ (0.0016 * (Length ^ 2.914)) * Abundance, 
      Genus == "Neocleon" ~ (0.0076 * (Length ^ 2.691)) * Abundance,# Baetis formula
      Genus == "Neophylax" ~ (0.0049 * (Length ^ 2.85)) * Abundance,
      Genus == "Neoplasta" ~ (0.004 * (Length ^ 2.655)) * Abundance,
      Genus == "Nigronia" ~ (0.0062 * (Length ^ 2.691)) * Abundance,
      Genus == "Optioservus" ~ (0.0039 * (Length ^ 2.96)) * Abundance,
      Genus == "Oreogeton" ~ (0.0033 * (Length ^ 2.392)) * Abundance,
      Genus == "Orthocladine" ~ (0.002 * (Length ^ 2.254)) * Abundance,
      Genus == "Oulimnius" ~ (0.0138 * (Length ^ 2.5548)) * Abundance,
      Genus == "Paracapnia" ~ (0.004 * (Length ^ 2.487)) * Abundance,
      Genus == "Paraleptophlebia" ~ (0.0038 * (Length ^ 2.918)) * Abundance,
      Genus == "Polycentropodidae" ~ (0.0071 * (Length ^ 2.531)) * Abundance,
      Genus == "Polycentropus" ~ (0.0071 * (Length ^ 2.531)) * Abundance,
      Genus == "Probezzia" ~ (0.0033 * (Length ^ 2.392)) * Abundance,
      Genus == "Prodaticus" ~ (0.1029 * (Length ^ 0)) * Abundance,
      Genus == "Prosimulium" ~ (0.0012 * (Length ^ 3.19)) * Abundance, 
      Genus == "Prostoia" ~ (0.004 * (Length ^ 2.975)) * Abundance,
      Genus == "Psephenus" ~ (0.0077 * (Length ^ 2.883)) * Abundance,
      Genus == "Pseudolimnophila" ~ (0.0014 * (Length ^ 2.667)) * Abundance,
      Genus == "Pteronarcys" ~ (0.0064 * (Length ^ 2.845)) * Abundance,
      Genus == "Pycnopsyche" ~ (0.0049 * (Length ^ 2.85)) * Abundance,
      Genus == "Psychodini" ~ (0.0007 * (Length ^ 2.592)) * Abundance, # Chironomini equation
      Genus == "Remenus" ~ (0.0119 * (Length ^ 2.695)) * Abundance,
      Genus == "Rhagovelia" ~ (0.0083 * (Length ^ 2.777)) * Abundance,
      Genus == "Rhyacophila" ~ (0.0099 * (Length ^ 2.48)) * Abundance,
      Genus == "Sialis" ~ (0.0031 * (Length ^ 2.801)) * Abundance,
      Genus == "Simulium" ~ (0.004 * (Length ^ 2.807)) * Abundance,
      Genus == "Stenelmis" ~ (0.0111 * (Length ^ 2.49)) * Abundance,
      Genus == "Stenonema" ~ (0.0128 * (Length ^ 2.616)) * Abundance, 
      Genus == "Stratiomyidae" ~ (0.005 * (Length ^ 2.591)) * Abundance,
      Genus == "Stylogomphus" ~ (0.0097 * (Length ^ 2.895)) * Abundance,
      Genus == "Tallaperla" ~ (0.0194 * (Length ^ 2.853)) * Abundance,
      Genus == "Taeniopteryx" ~ (0.0067 * (Length ^ 2.711)) * Abundance,
      Genus == "Tanypodinae" ~ (0.0038 * (Length ^ 0.0006)) * Abundance,
      Genus == "Tanytarsini" ~ (0.0008 * (Length ^ 2.728)) * Abundance,
      Genus == "Tipula" ~ (0.0064 * (Length ^ 2.443)) * Abundance, 
      Genus == "Triacanthagyna" ~ (0.0082 * (Length ^ 2.813)) * Abundance, # Using Boyeria from Benke
      Genus == "Wormaldia" ~ (0.0044 * (Length ^ 2.652)) * Abundance,
      Genus == "Zoraena" ~ (0.0067 * (Length ^ 2.782)) * Abundance, # Family Cordulegaster
      TRUE ~ NA_real_  # Assign NA for genera not specified
    ))

# checking for NAs, should be zero!
sum(is.na(SECPROD$Biomass.mg))

      
      # Finally, let's add a new Density column, then use it to correct biomass by area
      SECPROD <- SECPROD %>% 
        mutate(Density = Abundance / 0.0929) %>% # Making density column
        mutate(Biomass.g = Biomass.mg / 1000) %>% # Biomass was in mg bc of the length mass regressions, divide by 1000 to get to g
        mutate(Biomass.Area.Corrected = Biomass.g*Density) # Making biomass.area.corrected column
      # Saving as a CSV for geom_ridge code
      write.csv(SECPROD, "SEC_PROD.csv", row.names = FALSE)