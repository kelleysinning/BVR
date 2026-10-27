# Kelley Sinning 9/15/2026
# Cleaning up all the manjor datasheets so they can be used in an analysis together!

# trying something

library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
install.packages("vegan")
library(vegan)
library(stringr)

setwd("~/Library/CloudStorage/OneDrive-TheUniversityofMontana/Data/BVR/Data Cleaning")

original_diets <- read.csv("Diet_Data.csv")
bugs <- read.csv("BMI_SurberData.csv")


# First, we will start with cleaning up the diets data sheet
meas_cols <- paste0("measurement_", 1:10)
keys <- c("Sample_ID", "Order", "Family", "Life_stage", "BL.HW")

diets <- original_diets %>%
  select(-Initials, -Date_entered, -IDer, -ID_date, -Occasion, -Diet_observation_number,
         -Measurement_mm, -Extra_counts, -Empty_case_mm, -PROOFED, -USE_SAMPLE, -Total_Sampled,
         -NOTES) %>% # removing extra columns
  mutate(`Life_stage` = if_else(Life_stage %in% c("l","L "), "L", Life_stage)) %>% # Convert lowercase "l" and "L " to capital "L" in Life Stage column
  filter(`Life_stage` %in% c("A", "L", "P")) %>% # Only keep Life Stages that are A, L, or P, ignoring cases and other random letters
  mutate(BL.HW = case_when(
    BL.HW %in% c("Bl", "bl", "BL") ~ "BL",
    BL.HW %in% c("hw", "Hw", "HW") ~ "HW",
    TRUE ~ BL.HW
  )) %>%
  filter(BL.HW %in% c("BL", "HW")) %>%
  mutate(across(all_of(meas_cols), ~ as.numeric(as.character(.x))), 
         Total.Measured = as.numeric(as.character(Total.Measured)),
          Sample_ID = str_replace(Sample_ID, "^C0525(\\d{3})", "C20525\\1")) # cleaning up this sample ID error

diets$Sample_date <- parse_date_time(as.character(diets$Sample_date), orders = c("mdy", "mdY")) %>% as.Date() # putting sample dates in chronological order
diets$Sample_date <- factor(diets$Sample_date, levels = sort(unique(diets$Sample_date)))

    # 1. Pool the measurements from duplicate rows, keeping the first 10 in row order
pooled <- diets %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(all_of(meas_cols), names_to = "slot", values_to = "value",
               values_drop_na = TRUE) %>%
  mutate(slot_n = as.integer(str_remove(slot, "measurement_"))) %>% # what number measurement was it? 1, or 10?
  arrange(row_id, slot_n) %>% 
  group_by(across(all_of(keys))) %>% 
  mutate(new_slot = row_number()) %>%
  filter(new_slot <= 10) %>%          # measurements beyond 10 are dropped here
  ungroup() %>%
  select(all_of(keys), new_slot, value) %>%
  pivot_wider(names_from = new_slot, values_from = value, names_prefix = "measurement_") # removing temporary "slot" and "row" 
# columns, but putting their values under measurements where they belong

    # make sure all 10 columns exist even if no group fills every slot
pooled[setdiff(meas_cols, names(pooled))] <- NA_real_
pooled <- pooled %>% select(all_of(keys), all_of(meas_cols))

    # 2. Collapse duplicates to one row; counts are summed so extras still count toward abundance
diets <- diets %>%
  group_by(across(all_of(keys))) %>%
  summarise(
    across(-c(all_of(meas_cols), Total.Measured), first),
    Total.Measured = sum(Total.Measured, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pooled, by = keys)

    # 3. Summing measurements into one column and renaming Abundance column
diets <- diets %>%
  mutate(Measurement_mean_mm = rowMeans(select(., all_of(meas_cols)), na.rm = TRUE)) %>%
  select(-all_of(meas_cols)) %>%
  rename(Abundance = Total.Measured)


# Now, fixing duplicate site IDs that occur because the sample location was entered incorrectly
diets <- diets %>%
  mutate(
    Sample_ID = toupper(trimws(Sample_ID)),
    Site_ID = str_extract(Sample_ID, "^(C1\\.5|C2|T1|T2|T3)")
  )

    # Check: any IDs that don't start with one of the five codes?
diets %>% filter(is.na(Site_ID)) %>% distinct(Sample_ID)

    # Check: how many rows disagree with the ID prefix?
diets %>% filter(Sample_location != Site_ID) %>% distinct(Sample_ID, Sample_location, Site_ID)

    # Keep only rows where the location matches the ID prefix
diets <- diets %>%
  mutate(Sample_location = Site_ID)

# Confirm each Sample_ID now has exactly one location
diets %>% distinct(Sample_ID, Sample_location) %>% count(Sample_ID) %>% filter(n > 1)



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
                                "Crangonyctidae" ~ "Amphipoda",
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
      BL.HW == "BL" & Taxon == "Ephemeroptera" ~ (0.0066 * (Measurement_mean_mm ^ 2.88)) * Abundance,
      BL.HW == "BL" & Taxon == "Trichoptera" ~ (0.0019 * (Measurement_mean_mm ^ 3.12)) * Abundance,
      BL.HW == "BL" & Taxon == "Diptera" ~ (0.00096 * (Measurement_mean_mm ^ 3)) * Abundance, # using larva
      BL.HW == "BL" & Taxon == "Plecoptera" ~ (0.0023 * (Measurement_mean_mm ^ 2.45)) * Abundance,
      BL.HW == "BL" & Taxon == "Terrestrial" ~ (0.04142 * (Measurement_mean_mm ^ 2.213)) * Abundance, # dipteran adults
      BL.HW == "BL" & Taxon == "Ostracoda" ~ (0.0484 * (Measurement_mean_mm ^ 1.943)) * Abundance,
      BL.HW == "BL" & Taxon == "Asellidae" ~ (0.0072 * (Measurement_mean_mm ^ 2.785)) * Abundance,
      BL.HW == "BL" & Taxon == "Coleoptera" ~ (0.0035 * (Measurement_mean_mm ^ 2.4033)) * Abundance,
      BL.HW == "BL" & Taxon == "Gastropoda" ~ (0.172 * (Measurement_mean_mm ^ 1.688)) * Abundance,
      BL.HW == "BL" & Taxon == "Turbellaria" ~ (0.0089 * (Measurement_mean_mm ^ 2.145)) * Abundance,
      BL.HW == "BL" & Taxon == "Odonata" ~ (0.01399 * (Measurement_mean_mm ^ 2.78)) * Abundance,
      BL.HW == "BL" & Taxon == "Oligochaeta" ~ (0.00241 * (Measurement_mean_mm ^ 1.875)) * Abundance,
      BL.HW == "BL" & Taxon == "Hemiptera" ~ (0.00836 * (Measurement_mean_mm ^ 3.075)) * Abundance,
      BL.HW == "BL" & Taxon == "Araneae" ~ (0.1044 * (Measurement_mean_mm ^ 2.296)) * Abundance,
      BL.HW == "BL" & Taxon == "Leech" ~ (0.0071 * (Measurement_mean_mm ^ 2.531)) * Abundance, # come back
      BL.HW == "BL" & Taxon == "Hymenoptera" ~ (0.01379 * (Measurement_mean_mm ^ 2.696)) * Abundance,
      BL.HW == "BL" & Taxon == "Lepidoptera" ~ (0.00271 * (Measurement_mean_mm ^ 2.959)) * Abundance, # larval vs adult?
      BL.HW == "BL" & Taxon == "Baetidae" ~ (0.0076 * (Measurement_mean_mm ^ 2.691)) * Abundance,
      BL.HW == "BL" & Taxon == "Lepidostomatidae" ~ (0.0079 * (Measurement_mean_mm ^ 2.649)) * Abundance,
      BL.HW == "BL" & Taxon == "Brachycentridae" ~ (0.0024 * (Measurement_mean_mm ^ 3.676)) * Abundance, # check old code
      BL.HW == "BL" & Taxon == "Rhyacophilidae" ~ (0.0024 * (Measurement_mean_mm ^ 3.676)) * Abundance, # check old code
      BL.HW == "BL" & Taxon == "Perlodidae" ~ (0.01 * (Measurement_mean_mm ^ 2.658)) * Abundance,
      BL.HW == "BL" & Taxon == "Chironomidae" ~ (0.0006 * (Measurement_mean_mm ^ 2.77)) * Abundance,
      BL.HW == "BL" & Taxon == "Chloroperlidae" ~ (0.0062 * (Measurement_mean_mm ^ 2.724)) * Abundance,
      BL.HW == "BL" & Taxon == "Tabanidae" ~ (0.005 * (Measurement_mean_mm ^ 2.591)) * Abundance,
      BL.HW == "BL" & Taxon == "Ephemerellidae" ~ (0.00928 * (Measurement_mean_mm ^ 2.9)) * Abundance,
      BL.HW == "BL" & Taxon == "Glossosomatidae" ~ (0.0024 * (Measurement_mean_mm ^ 2.616)) * Abundance, # old code
      BL.HW == "BL" & Taxon == "Heptageniidae" ~ (0.0128 * (Measurement_mean_mm ^ 2.616)) * Abundance,
      BL.HW == "BL" & Taxon == "Hydropsychidae" ~ (0.0049 * (Measurement_mean_mm ^ 2.62)) * Abundance,
      BL.HW == "BL" & Taxon == "Simuliidae" ~ (0.0048 * (Measurement_mean_mm ^ 2.55)) * Abundance,
      BL.HW == "BL" & Taxon == "Tipulidae" ~ (0.00392 * (Measurement_mean_mm ^ 2.4403)) * Abundance,
      BL.HW == "BL" & Taxon == "Perlidae" ~ (0.003 * (Measurement_mean_mm ^ 3.232)) * Abundance,
      BL.HW == "BL" & Taxon == "Athericidae" ~ (0.0024 * (Measurement_mean_mm ^ 3.676)) * Abundance, # check old code
      BL.HW == "BL" & Taxon == "Elmidae" ~ (0.0111 * (Measurement_mean_mm ^ 2.49)) * Abundance,
      BL.HW == "BL" & Taxon == "Formicidae" ~ (0.00885 * (Measurement_mean_mm ^ 2.919)) * Abundance,
      BL.HW == "BL" & Taxon == "Hydrophilidae" ~ (0.0024 * (Measurement_mean_mm ^ 2.2)) * Abundance,
      BL.HW == "BL" & Taxon == "Gomphidae" ~ (0.0044 * (Measurement_mean_mm ^ 3.124)) * Abundance,
      BL.HW == "BL" & Taxon == "Dytiscidae" ~ (0.1029 * (Measurement_mean_mm ^ 0)) * Abundance, # exponent 0, check this
      BL.HW == "BL" & Taxon == "Capniidae" ~ (0.004 * (Measurement_mean_mm ^ 2.487)) * Abundance,
      BL.HW == "BL" & Taxon == "Cicadellidae" ~ (0.02387 * (Measurement_mean_mm ^ 2.561)) * Abundance, # leaf hoppers
      BL.HW == "BL" & Taxon == "Leptophlebiidae" ~ (0.0054 * (Measurement_mean_mm ^ 2.836)) * Abundance,
      BL.HW == "BL" & Taxon == "Physidae" ~ (0.172 * (Measurement_mean_mm ^ 1.688)) * Abundance, # using gastropoda
      BL.HW == "BL" & Taxon == "Hydroptilidae" ~ (0.01268 * (Measurement_mean_mm ^ 2.901)) * Abundance,
      BL.HW == "BL" & Taxon == "Muscidae" ~ (0.00033 * (Measurement_mean_mm ^ 3.55)) * Abundance,
      BL.HW == "BL" & Taxon == "Aphididae" ~ (0.0598 * (Measurement_mean_mm ^ 1.724)) * Abundance,
      BL.HW == "BL" & Taxon == "Nemouridae" ~ (0.004 * (Measurement_mean_mm ^ 2.975)) * Abundance,
      BL.HW == "BL" & Taxon == "Pteronarcyidae" ~ (0.0064 * (Measurement_mean_mm ^ 2.845)) * Abundance,
      BL.HW == "BL" & Taxon == "Limoniidae" ~ (0.00392 * (Measurement_mean_mm ^ 2.4403)) * Abundance,
      BL.HW == "BL" & Taxon == "Stratiomyidae" ~ (0.005 * (Measurement_mean_mm ^ 2.591)) * Abundance,
      BL.HW == "BL" & Taxon == "Empididae" ~ (0.004 * (Measurement_mean_mm ^ 2.655)) * Abundance,
      BL.HW == "BL" & Taxon == "Amphipoda" ~ (0.0058 * (Measurement_mean_mm ^ 2.798)) * Abundance,
      # No coefficients (NA): Planariidae, Hydrachnidae, 
      # Fry, Fish egg, Cyclorrhapha, Notonectidae, Vespidae
      TRUE ~ NA_real_  # NA for other taxa, HW rows, or missing BL.HW
    ))

# Adding headwidth biomass where possible
diets <- diets %>%
  mutate(
    Biomass.mg = case_when(
      BL.HW == "HW" & Taxon == "Chironomidae"      ~ (2.7842 * (Measurement_mean_mm ^ 2.835)) * Abundance,
      BL.HW == "HW" & Taxon == "Baetidae"          ~ (0.815 * (Measurement_mean_mm ^ 3.349)) * Abundance,
      BL.HW == "HW" & Taxon == "Diptera"           ~ (2.7842 * (Measurement_mean_mm ^ 2.835)) * Abundance,
      BL.HW == "HW" & Taxon == "Ephemeroptera"     ~ (0.815 * (Measurement_mean_mm ^ 3.349)) * Abundance,
      BL.HW == "HW" & Taxon == "Trichoptera"       ~ (2.221 * (Measurement_mean_mm ^ 3.349)) * Abundance,
      BL.HW == "HW" & Taxon == "Simuliidae"        ~ (2.553 * (Measurement_mean_mm ^ 4.347)) * Abundance,
      BL.HW == "HW" & Taxon == "Heptageniidae"     ~ (0.060 * (Measurement_mean_mm ^ 4.111)) * Abundance,
      BL.HW == "HW" & Taxon == "Lepidostomatidae"  ~ (1.666 * (Measurement_mean_mm ^ 2.987)) * Abundance,
      BL.HW == "HW" & Taxon == "Hydropsychidae"    ~ (0.984 * (Measurement_mean_mm ^ 2.814)) * Abundance,
      BL.HW == "HW" & Taxon == "Ephemerellidae"    ~ (0.450 * (Measurement_mean_mm ^ 3.476)) * Abundance,
      BL.HW == "HW" & Taxon == "Brachycentridae"   ~ (2.221 * (Measurement_mean_mm ^ 3.349)) * Abundance,
      BL.HW == "HW" & Taxon == "Plecoptera"        ~ (0.3208 * (Measurement_mean_mm ^ 3.189)) * Abundance,
      BL.HW == "HW" & Taxon == "Rhyacophilidae"    ~ (1.750 * (Measurement_mean_mm ^ 3.522)) * Abundance,
      BL.HW == "HW" & Taxon == "Perlodidae"        ~ (0.5462 * (Measurement_mean_mm ^ 2.826)) * Abundance,
      BL.HW == "HW" & Taxon == "Perlidae"          ~ (0.3208 * (Measurement_mean_mm ^ 3.189)) * Abundance,
      BL.HW == "HW" & Taxon == "Asellidae"         ~ (0.6525 * (Measurement_mean_mm ^ 3.001)) * Abundance,
      BL.HW == "HW" & Taxon == "Elmidae"           ~ (1.4040 * (Measurement_mean_mm ^ 3.794)) * Abundance,
      BL.HW == "HW" & Taxon == "Coleoptera"        ~ (1.4040 * (Measurement_mean_mm ^ 3.794)) * Abundance,
      BL.HW == "HW" & Taxon == "Gomphidae"         ~ (0.8177 * (Measurement_mean_mm ^ 2.454)) * Abundance,
      BL.HW == "HW" & Taxon == "Amphipoda"         ~ (1.091 * (Measurement_mean_mm ^ 3.891)) * Abundance,
      # Not included (no BL equation either, or not a real taxon):
      #  Fry, Vespidae, Insecta, Empty
      TRUE ~ Biomass.mg
    ))
# These taxa had no HW equations in Benke et al. 1999
#BL.HW == "HW" & Taxon == "Glossosomatidae"   ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Tipulidae"         ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Formicidae"        ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Terrestrial"       ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Dytiscidae"        ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Leptophlebiidae"   ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Araneae"           ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Hydrophilidae"     ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Oligochaeta"       ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Cicadellidae"      ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Gastropoda"        ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Hymenoptera"       ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Hemiptera"         ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Hydroptilidae"     ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Leech"             ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Lepidoptera"       ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,
#BL.HW == "HW" & Taxon == "Ostracoda"         ~ (NA_real_ * (Measurement_mean_mm ^ NA_real_)) * Abundance,

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

# Looking at % FFG
diets <- diets %>%
  group_by(Sample_ID) %>%
  mutate(Percent.Scraper.Abundance = 100 * sum(Abundance[FFG == "Scraper"], na.rm = TRUE)
         / sum(Abundance, na.rm = TRUE),
         Percent.Scraper.Biomass = 100 * sum(Biomass.Area.Corrected[FFG == "Scraper"], na.rm = TRUE)
         / sum(Biomass.Area.Corrected, na.rm = TRUE)) %>%
  ungroup()

# A little bit more cleaning up
diets <- diets %>%
  select(-Order, -Family, - Density, -Measurement_mean_mm) %>%
  filter(Taxon != "Empty") %>%
  arrange(Sample_date) %>%
  mutate(Sample_date = factor(Sample_date, levels = sort(unique(Sample_date)))) %>%
  relocate(Sample_ID, Sample_location, Sample_date, Life_stage, BL.HW, Taxon, Abundance, Biomass.Area.Corrected, FFG)

str(diets)
# First, we want to crunch this down to one fish with summed abundance, biomass, FFG metrics, diversity

# Ultimately, we'd want to join this with fish and say if it has a diet it also has a column for diversity of diet, 
# # of diet items, % scrapers

# Next, summarize by sample occasions left join diets with fish ID codes to assign species

# CALCULATING VEGAN INDICIES------------------------------------------------------

community <- diets %>%
  group_by(Sample_ID, Taxon) %>% # not grouping by BL or HW rn
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Taxon, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("Sample_ID")

str(community)
# Richness (number of different taxa)
S <- specnumber(community)

# Shannon, Simpson, inverse Simpson
H     <- diversity(community, index = "shannon") # Higher values mean greater diversity, sensitive to rare species
simp  <- diversity(community, index = "simpson") # A higher value means lower diversity (higher dominance by one or a few species).
invS  <- diversity(community, index = "invsimpson") # Higher values mean greater diversity. It represents the "effective number of species"—the number of equally-abundant species that would give the same community diversity

# Evenness (Pielou's J)
J <- H / log(S) # how similar is the number of different species, or balance



div_table <- data.frame(
  Sample_ID = rownames(community),
  richness = S, shannon = H, simpson = simp,
  invsimpson = invS, evenness = J
)

# FINALIZING THE DIET TABLE-----------------------------------------------------

diets_final <- diets %>%
  group_by(Sample_ID, Sample_location) %>%
  summarise(
    Total.Biomass = sum(Biomass.Area.Corrected),
    Total.Abundance = sum(Abundance),
    Percent.Scraper.Abundance = mean(Percent.Scraper.Abundance),
    Percent.Scraper.Biomass = mean(Percent.Scraper.Biomass),
    .groups = "drop"
  ) %>%
  left_join(div_table, by = "Sample_ID")

  
# LINKING WITH FISH DATA--------------------------------------------------------
fish <- read.csv("Meta_Field_Data.csv")

fish <- fish %>%
  select(-Initials, -Date_entered, -Observation_number, -Occasion, -Recorder, -Pass,
         -Fin_clip, -Scales, -Handling_mortality, -Comments, -Recap, -PIT_ID, -Handler,
         -Proofed_by_initials, -Proofed_date, -Corrected_by_initials, -Corrected_date) %>%
  rename(Sample_ID = Individual_ID)

fish <- fish %>%
  left_join(diets_final, by = "Sample_ID") %>%
  select(-Sample_location)

fish_with_diets <- fish %>% # making a df with all fish diet info
  filter(Diet == "Y", !is.na(Total.Abundance))

plot(Weight~Fork_length , data = fish)

ggplot(fish_with_diets, aes(x = Species, y = Percent.Scraper.Abundance)) +
       geom_boxplot() 
       # + facet_wrap(~Sampling_date)

# Adding in relative weight (I should check these equations)
fish <- fish %>%
  mutate(Fork_length = as.numeric(Fork_length),
         Weight      = as.numeric(Weight)) %>%
  mutate(Ws = case_when(
    Species == "BNT" ~ 10^(-4.867 + 2.960*log10(Fork_length)), # Milewksi and Brown et al, 1994 but doesn't apply to small fish, min 120 cm, can fit a linear model for our own pooled data, take log10 of. each length and weight, each fish above and below are different
    Species == "MTS" ~ 10^(-5.202 + 3.16*log10(Fork_length))
  ),
  Wr = Weight/Ws * 100)



fish_with_diets <- fish_with_diets %>%
  mutate(Fork_length = as.numeric(Fork_length),
         Weight      = as.numeric(Weight)) %>%
  mutate(Ws = case_when(
    Species == "BNT" ~ 10^(-4.867 + 2.960*log10(Fork_length)),
    Species == "MTS" ~ 10^(-5.202 + 3.16*log10(Fork_length))
  ),
  Wr = Weight/Ws * 100)


# Left join with flow

# RETRIEVING USGS GAUGE DATA
install.packages("dataRetrieval", type = "source")
library(dataRetrieval)
packageVersion("dataRetrieval")  # should be >= 2.7.19


# https://doi-usgs.github.io/dataRetrieval/reference/read_waterdata_daily.html
# read_waterdata_daily is newer 
site <- "USGS-09057500"  # USGS site number (Blue River below Green Mountain)
parameter_code <- "00060"    # Parameter code for discharge (ft³/s)
statistic_id <- "00003"      # Statistic code for daily mean
start_date <- "2021-04-01" # Date before first sampling event ever
end_date <- "2026-10-10" # Recent data

# Retrieve daily discharge data
discharge_data <- read_waterdata_daily(
  monitoring_location_id = site,
  parameter_code = parameter_code,
  statistic_id = statistic_id,
  time = c(start_date, end_date),
  skipGeometry = TRUE,
)

discharge_data <- discharge_data %>%
  select(time, value)

discharge_data <- discharge_data %>%
  rename(Discharge_cfs = value,
         Sampling_date = time
  )

discharge_data <- as.data.frame(discharge_data)
str(fish)

fish <- fish %>%
  mutate(Sampling_date = mdy(Sampling_date))%>%
  left_join(discharge_data, by = "Sampling_date") 

fish_with_diets <- fish_with_diets %>%
  mutate(Sampling_date = mdy(Sampling_date))%>%
  left_join(discharge_data, by = "Sampling_date") 





# Okay, now I want to add in surber community------------------------------------
library(dplyr)

original_bugs<- bugs


meas_cols <- paste0("measurement_", 1:10)

bugs <- bugs %>%
  dplyr::select(-initials, -date_entered, -sorter_ID, -occasion_ID, -collection_date, -sorting_date,
         -grids_sampled, -empty_case, -NOTES, -X, -PROOFED.BY, -PROOF.DATE) %>% # removing extra columns
  filter(!grepl("Eagle River|CR|ER", Site, ignore.case = TRUE)) %>%
  mutate(life_stage = case_when(
    life_stage %in% c("L ", "L?") ~ "L",
    life_stage %in% c("P?") ~ "P",
    life_stage %in% c("A?") ~ "A",
    TRUE ~ life_stage)) %>%
  filter(`life_stage` %in% c("A", "L", "P")) %>% # Only keep Life Stages that are A, L, or P, ignoring random spaces
  mutate(BL_HW = case_when(
    BL_HW %in% c("Bl", "BW") ~ "BL",
    TRUE ~ BL_HW
  )) %>%
  filter(BL_HW %in% c("BL", "HW")) # ignoring ""

# also want to filter out eagle and colorado river
unique(bugs$sample_location)

bugs$Sample_date <- parse_date_time(as.character(bugs$Sample_date), orders = c("mdy", "mdY")) %>% as.Date() # putting sample dates in chronological order
bugs$Sample_date <- factor(bugs$Sample_date, levels = sort(unique(bugs$Sample_date)))

# 1. Pool the measurements from duplicate rows, keeping the first 10 in row order
pooled <- bugs %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(all_of(meas_cols), names_to = "slot", values_to = "value",
               values_drop_na = TRUE) %>%
  mutate(slot_n = as.integer(str_remove(slot, "measurement_"))) %>% # what number measurement was it? 1, or 10?
  arrange(row_id, slot_n) %>% 
  group_by(across(all_of(keys))) %>% 
  mutate(new_slot = row_number()) %>%
  filter(new_slot <= 10) %>%          # measurements beyond 10 are dropped here
  ungroup() %>%
  select(all_of(keys), new_slot, value) %>%
  pivot_wider(names_from = new_slot, values_from = value, names_prefix = "measurement_") # removing temporary "slot" and "row" 
# columns, but putting their values under measurements where they belong

# make sure all 10 columns exist even if no group fills every slot
pooled[setdiff(meas_cols, names(pooled))] <- NA_real_
pooled <- pooled %>% select(all_of(keys), all_of(meas_cols))

# 2. Collapse duplicates to one row; counts are summed so extras still count toward abundance
bugs <- bugs %>%
  group_by(across(all_of(keys))) %>%
  summarise(
    across(-c(all_of(meas_cols), Total.Measured), first),
    Total.Measured = sum(Total.Measured, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pooled, by = keys)

# 3. Summing measurements into one column and renaming Abundance column
bugs <- bugs %>%
  mutate(Measurement_mean_mm = rowMeans(select(., all_of(meas_cols)), na.rm = TRUE)) %>%
  select(-all_of(meas_cols)) %>%
  rename(Abundance = Total.Measured)
