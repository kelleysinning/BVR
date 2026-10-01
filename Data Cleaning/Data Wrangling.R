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
         -NOTES) %>% # With this move, we are disregarding empty cases, look deeper into this and maybe delete the whole row if this is a thing
  mutate(across(measurement_1:measurement_10, ~ as.numeric(as.character(.x)))) %>%
  mutate(Measurement_mean_mm= rowMeans(select(., measurement_1:measurement_10), na.rm = TRUE)) %>%
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
                            c("Trombiditonnes", "Trombidiformes", "Ant") ~ "Formicidae",
                            c("Hirudinea") ~ "Leech",
                            c("gastropoda") ~ "Gastropoda",
                            c("Polychaeta") ~ "Oligochaeta",
                            "Crangonyctidae" ~ "Amphipoda",
                            c("Turbellaria","Planariidae") ~ "Planarian",
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
                                "Planariidae" ~ "Planarian",
                                default = Family))

unique(diets$Family) # check this for more fixes as data is added    

# Making a new column for highest taxonomic resolution
diets <- diets %>%
  mutate(
    Family = na_if(trimws(Family), ""),
    Order  = na_if(trimws(Order), ""),
    Taxon  = coalesce(Family, Order)
  )

# Adding in length-mass equations; assigning biomass to highest taxonomic resolution
diets <- diets %>% # First, convert these to numeric
  mutate(
    Measurement_mean_mm = as.numeric(as.character(Measurement_mean_mm)),
    Abundance = as.numeric(as.character(Abundance))
  )

diets <- diets %>%
  mutate(
    Biomass.mg = case_when(
      Taxon == "Ephemeroptera" ~ (0.0066 * (Measurement_mean_mm ^ 2.88)) * Abundance,
      Taxon == "Trichoptera" ~ (0.0019 * (Measurement_mean_mm ^ 3.12)) * Abundance,
      Taxon == "Diptera" ~ (0.00096 * (Measurement_mean_mm ^ 3)) * Abundance,
      Taxon == "Plecoptera" ~ (0.0023 * (Measurement_mean_mm ^ 2.45)) * Abundance,
      Taxon == "Terrestrial" ~ (0.04142 * (Measurement_mean_mm ^ 2.213)) * Abundance, # dipteran adults
      Taxon == "Ostracoda" ~ (0.0484 * (Measurement_mean_mm ^ 1.943)) * Abundance,
      Taxon == "Asellidae" ~ (0.0072 * (Measurement_mean_mm ^ 2.785)) * Abundance,
      Taxon == "Coleoptera" ~ (0.0035 * (Measurement_mean_mm ^ 2.4033)) * Abundance,
      Taxon == "Gastropoda" ~ (0.172 * (Measurement_mean_mm ^ 1.688)) * Abundance,
      Taxon == "Turbellaria" ~ (0.0089 * (Measurement_mean_mm ^ 2.145)) * Abundance,
      Taxon == "Odonata" ~ (0.01399 * (Measurement_mean_mm ^ 2.78)) * Abundance,
      Taxon == "Oligochaeta" ~ (0.00241 * (Measurement_mean_mm ^ 1.875)) * Abundance,
      Taxon == "Hemiptera" ~ (0.00836 * (Measurement_mean_mm ^ 3.075)) * Abundance,
      Taxon == "Araneae" ~ (0.1044 * (Measurement_mean_mm ^ 2.296)) * Abundance,
      Taxon == "Leech" ~ (0.0071 * (Measurement_mean_mm ^ 2.531)) * Abundance, # come back
      Taxon == "Hymenoptera" ~ (0.01379 * (Measurement_mean_mm ^ 2.696)) * Abundance,
      Taxon == "Lepidoptera" ~ (0.00271 * (Measurement_mean_mm ^ 2.959)) * Abundance, # put this as larval but is it more likely to be adult?
      #Taxon == "Planariidae" ~ ( * (Measurement_mean_mm ^ )) * Abundance, # NA
      #Taxon == "Hydrachnidae" ~ ( * (Measurement_mean_mm ^ )) * Abundance, # NA
      #Taxon == "Amphipoda" ~ ( * (Measurement_mean_mm ^ )) * Abundance,# NA
      
      Taxon == "Baetidae" ~ (0.0076 * (Measurement_mean_mm ^ 2.691)) * Abundance,
      Taxon == "Lepidostomatidae" ~ (0.0079 * (Measurement_mean_mm ^ 2.649)) * Abundance,
      Taxon == "Brachycentridae" ~ (0.0024 * (Measurement_mean_mm ^ 3.676)) * Abundance, # check old code
      Taxon == "Rhyacophilidae" ~ (0.0024 * (Measurement_mean_mm ^ 3.676)) * Abundance,# check old code
      Taxon == "Perlodidae" ~ (0.01 * (Measurement_mean_mm ^ 2.658)) * Abundance,
      Taxon == "Chironomidae" ~ (0.0006 * (Measurement_mean_mm ^ 2.77)) * Abundance,
      Taxon == "Chloroperlidae" ~ (0.0062 * (Measurement_mean_mm ^ 2.724)) * Abundance,
      Taxon == "Tabanidae" ~ (0.005 * (Measurement_mean_mm ^ 2.591)) * Abundance,
      Taxon == "Ephemerellidae" ~ (0.00928 * (Measurement_mean_mm ^ 2.9)) * Abundance,
      Taxon == "Glossosomatidae" ~ (0.0024 * (Measurement_mean_mm ^ 2.616)) * Abundance, #old code
      Taxon == "Heptageniidae" ~ (0.0128 * (Measurement_mean_mm ^ 2.616)) * Abundance,
      Taxon == "Hydropsychidae" ~ (0.0049 * (Measurement_mean_mm ^ 2.62)) * Abundance,
      Taxon == "Simuliidae" ~ (0.0048 * (Measurement_mean_mm ^ 2.55)) * Abundance,
      Taxon == "Tipulidae" ~ (0.00392 * (Measurement_mean_mm ^ 2.4403)) * Abundance,
      Taxon == "Perlidae" ~ (0.003 * (Measurement_mean_mm ^ 3.232)) * Abundance,
      Taxon == "Araneae" ~ (0.1044 * (Measurement_mean_mm ^ 2.296)) * Abundance,
      Taxon == "Athericidae" ~ (0.0024 * (Measurement_mean_mm ^ 3.676)) * Abundance, #check old code
      Taxon == "Asellidae" ~ (0.0072 * (Measurement_mean_mm ^ 2.785)) * Abundance,
      Taxon == "Elmidae" ~ (0.0111 * (Measurement_mean_mm ^ 2.49)) * Abundance,
      Taxon == "Formicidae" ~ (0.00885 * (Measurement_mean_mm ^ 2.919)) * Abundance,
      Taxon == "Hemiptera" ~ (0.00836 * (Measurement_mean_mm ^ 3.075)) * Abundance,
      Taxon == "Hydrophilidae" ~ (0.0024 * (Measurement_mean_mm ^ 2.2)) * Abundance,
      Taxon == "Gomphidae" ~ (0.0044 * (Measurement_mean_mm ^ 3.124)) * Abundance,
      Taxon == "Dytiscidae" ~ (0.1029 * (Measurement_mean_mm ^ 0)) * Abundance,
      Taxon == "Capniidae" ~ (0.004 * (Measurement_mean_mm ^ 2.487)) * Abundance,
      Taxon == "Cicadellidae" ~ (0.02387 * (Measurement_mean_mm ^ 2.561)) * Abundance, # leaf hoppers
      Taxon == "Coleoptera" ~ (0.0035 * (Measurement_mean_mm ^ 2.4033)) * Abundance,
      Taxon == "Leptophlebiidae" ~ (0.0054 * (Measurement_mean_mm ^ 2.836)) * Abundance,
      Taxon == "Physidae" ~ (0.172 * (Measurement_mean_mm ^ 1.688)) * Abundance, # using gastropoda
      Taxon == "Hydroptilidae" ~ (0.01268 * (Measurement_mean_mm ^ 2.901)) * Abundance,
      Taxon == "Muscidae" ~ (0.00033 * (Measurement_mean_mm ^ 3.55)) * Abundance,
      Taxon == "Diptera" ~ (0.00096 * (Measurement_mean_mm ^ 3)) * Abundance, # using larva
      Taxon == "Gastropoda" ~ (0.172 * (Measurement_mean_mm ^ 1.688)) * Abundance,
      Taxon == "Terrestrial" ~ (0.04142 * (Measurement_mean_mm ^ 2.213)) * Abundance, # dipteran adults
      Taxon == "Oligochaeta" ~ (0.00241 * (Measurement_mean_mm ^ 1.875)) * Abundance,
      Taxon == "Aphididae" ~ (0.0598 * (Measurement_mean_mm ^ 1.724)) * Abundance,
      Taxon == "Nemouridae" ~ (0.004 * (Measurement_mean_mm ^ 2.975)) * Abundance,
      Taxon == "Pteronarcyidae" ~ (0.0064 * (Measurement_mean_mm ^ 2.845)) * Abundance,
      Taxon == "Limoniidae" ~ (0.00392 * (Measurement_mean_mm ^ 2.4403)) * Abundance, 
      Taxon == "Stratiomyidae" ~ (0.005 * (Measurement_mean_mm ^ 2.591)) * Abundance,
      Taxon == "Empididae" ~ (0.004 * (Measurement_mean_mm ^ 2.655)) * Abundance,
      #Taxon == "Crangonyctidae" ~ ( * (Measurement_mean_mm ^ )) * Abundance, # NA
      #Taxon == "Fry" ~ ( * (Measurement_mean_mm ^ )) * Abundance, # NA
      #Taxon == "Fish egg" ~ ( * (Measurement_mean_mm ^ )) * Abundance, # NA
      #Taxon == "Hydrachnidae" ~ ( * (Measurement_mean_mm ^ )) * Abundance, # NA
      #Taxon == "Cyclorrhapha" ~ ( * (Measurement_mean_mm ^ )) * Abundance, # NA
      #Taxon == "Notonectidae" ~ ( * (Measurement_mean_mm ^ )) * Abundance, # NA
      #Taxon == "Vespidae" ~ ( * (Measurement_mean_mm ^ )) * Abundance, #NA
      #Taxon == "Planariidae" ~ ( * (Measurement_mean_mm ^ )) * Abundance, # NA
      TRUE ~ NA_real_  # Assign NA for genera not specified
    ))

      
# Now, let's add a new Density column, then use it to correct biomass by area
diets <- diets %>% 
        mutate(Density = Abundance / 0.0929) %>% # Making density column based on 30 cm x 30 cm surber area --> .09 m^2
        mutate(Biomass.g = Biomass.mg / 1000) %>% # Biomass was in mg bc of the Measurement_mean_mm mass regressions, divide by 1000 to get to g
        mutate(Biomass.Area.Corrected = Biomass.g*Density) %>% # Making biomass.area.corrected column
        #filter(!is.na(Biomass.Area.Corrected)) %>% # Removing NAs bc that'd mean data was never entered, casualty of a messy dataset
        select(-Biomass.mg, -Biomass.g)

# How to clean this up? Should we?
# No biomass just means there wasn't an equation for it, which is a good bit of taxa
# No abundance means its empty, which is information

# Lovely!!
# Now, let's assign some coarse FFGs

diets <- diets %>%
  mutate(     FFG = case_when(
    # Collector-gatherers
    Taxon %in% c("Ephemeroptera", "Baetidae", "Leptophlebiidae", "Ephemerellidae",
                 "Diptera", "Chironomidae", "Stratiomyidae", "Elmidae",
                 "Ostracoda", "Oligochaeta") ~ "Collector-gatherer",
    
    # Collector-filterers
    Taxon %in% c("Hydropsychidae", "Brachycentridae", "Simuliidae") ~ "Collector-filterer",
    
    # Scrapers
    Taxon %in% c("Heptageniidae", "Glossosomatidae", "Hydroptilidae",
                 "Gastropoda", "Physidae") ~ "Scraper",
    
    # Shredders
    Taxon %in% c("Lepidostomatidae", "Tipulidae", "Limoniidae", "Capniidae",
                 "Nemouridae", "Pteronarcyidae", "Asellidae",
                 "Amphipoda", "Lepidoptera") ~ "Shredder",
    
    # Predators
    Taxon %in% c("Rhyacophilidae", "Perlodidae", "Chloroperlidae", "Perlidae",
                 "Plecoptera", "Tabanidae", "Athericidae", "Empididae", "Muscidae",
                 "Dytiscidae", "Hydrophilidae", "Coleoptera",
                 "Gomphidae", "Odonata", "Hemiptera",
                 "Turbellaria", "Leech", "Araneae") ~ "Predator",
    
    # Terrestrial inputs (not aquatic feeding groups)
    Taxon %in% c("Terrestrial", "Formicidae", "Hymenoptera", "Ant",
                 "Aphididae", "Cicadellidae") ~ "Terrestrial",
    
    # Taxa with no biomass equation (the commented-out lines above)
    Taxon == "Crangonyctidae" ~ "Shredder",     # amphipod; often gatherer/shredder
    Taxon == "Hydrachnidae"   ~ "Predator",     # water mites
    Taxon == "Cyclorrhapha"   ~ "Predator",     # mixed; adjust if mostly adults
    Taxon == "Notonectidae"   ~ "Predator",     # backswimmers
    Taxon == "Planariidae"    ~ "Predator",     # flatworms
    Taxon == "Vespidae"       ~ "Terrestrial",  # wasps
    Taxon == "Fry"            ~ "Predator",     # fish; change to NA if you'd rather exclude
    Taxon == "Fish egg"       ~ NA_character_,  # not a feeding group
    
    # Mixed at order level, needs a decision
    Taxon == "Trichoptera" ~ NA_character_,
    
    TRUE ~ NA_character_
  ))


# Ultimately, we'd want to join this with fish and say if it has a diet it also has a column for diversity of diet, 
# # of diet items, % scrapers

# Next, summarize by sample occasions left join diets with fish ID codes to assign species

