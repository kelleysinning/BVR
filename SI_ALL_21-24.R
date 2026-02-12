# BVR SI Bug and Fish Data for Kelley Sinning from May 2021 to present
# Focusing on macro and fish overlap

# Load packages
library(tidyverse)
library(SIBER)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)
library(car)
library(kableExtra)
library(lubridate)
library(vegan)
library(ggplot2)
library(dplyr)


setwd("~/Library/CloudStorage/OneDrive-Colostate/Data/BVR")

# IMPORT SIA DATA
SIA_ALL <- read_csv("Meta_SIA_Data.csv")

Sampling_dates <- read_csv("Sampling_dates.csv")

Sampling_dates <- Sampling_dates %>%
  mutate(Sampling_date = mdy(Sampling_date))

# Putting sampling occassions in order
SIA_ALL$Occasion <- factor(SIA_ALL$Occasion, levels = c("MAY_2021", "AUG_2021", "OCT_2021",
                                                          "MAY_2022", "AUG_2022", "OCT_2022",
                                                          "MAY_2023", "AUG_2023", "OCT_2023",
                                                          "MAY_2024", "AUG_2024", "OCT_2024")) # include "Unknown if you want to see the algae point


# Adding a sampling date column so it isn't just month/year 
# better for linking with specific hydrograph patterns
SIA_ALL <- SIA_ALL %>%
  left_join(Sampling_dates, by = c("Occasion", "Location"))

# NMDS----------------------------------------------------------------------------
head(SIA_ALL)


SIA_ALL <- SIA_ALL %>%
  select(d13C, d15N, Species, Occasion, Location) %>%
  drop_na()

str(SIA_ALL)

SIA_ALL_clean <- SIA_ALL %>% # Eliminating a bunch of family bug variation in permanova
  select(d13C, d15N, Species, Occasion, Location) %>%
  mutate(Species = if_else(Species == "FRY", "Fry", Species)) %>%
  mutate(
    Group = case_when(
      Species %in% c("BNT", "MTS", "Fry", "Fish Eggs") ~ Species,
      TRUE ~ "Macro"
    )
  ) %>%
  drop_na()


# Taking average so it is a lot cleaner
SIA_avg <- SIA_ALL_clean %>%
  group_by(Location, Occasion, Group) %>% # change from group to species if using SIA_ALL
  summarise(
    d13C = mean(d13C, na.rm = TRUE),
    d15N = mean(d15N, na.rm = TRUE),
    .groups = "drop"
  )


SIA_matrix <- SIA_avg %>% # Use SIA_ALL if you want chaos
  select(d13C, d15N)  
  #scale()                  # important: standardize isotopes

set.seed(123)  # reproducibility

NMDS_SIA <- metaMDS(
  SIA_matrix,
  distance = "euclidean",  # appropriate for isotope space
  k = 2,                   # 2D solution
  trymax = 10,
  autotransform = FALSE
)

NMDS_SIA


scores_nmds <- as.data.frame(scores(NMDS_SIA))


scores_nmds <- scores_nmds %>%
  bind_cols(SIA_avg %>% select(Group, Occasion, Location))%>%
  drop_na() 



scores_nmds$Occasion <- factor(scores_nmds$Occasion, levels = c("MAY_2021", "MAY_2022","MAY_2023","MAY_2024",
                                                                        "AUG_2021", "AUG_2022","AUG_2023","AUG_2024",
                                                                        "OCT_2021", "OCT_2022","OCT_2023","OCT_2024"))


library(rcartocolor)
display_carto_all()
mycolors <- carto_pal(9, "Geyser")
mycolors

species_colors <- c(
  "BNT" = "#008080",
  "MTS"   = "#DE8A5A",
  "Macro" = "#92B69E",
  "Algae"        = "#798234",
  "Fish Eggs"  = "#6C8EBF",
  "Fry" = "#DC869A"
)


ggplot(scores_nmds_clean, aes(x = NMDS1, y = NMDS2, color = Group)) +
  facet_wrap(~Occasion, scales = "free_y") +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = Group),
               geom = "polygon",
               alpha = 0.2,
               color = NA,
               level = 0.95) +
  scale_color_manual(values = species_colors) +
  scale_fill_manual(values = species_colors) +
  theme_classic() +
  labs(
    title = "NMDS of Stable Isotope Space",
    subtitle = paste("Stress =", round(NMDS_SIA$stress, 3)),
    x = "NMDS1",
    y = "NMDS2"
  )


# PERMANOVA

dispersion <- betadisper(
  dist(SIA_matrix),
  SIA_avg$Group
)

anova(dispersion)

boxplot(dispersion)

set.seed(123)

permanova <- adonis2(
  SIA_matrix ~ Group * Occasion,
  data = SIA_avg,
  method = "euclidean",
  permutations = 999
)

permanova # Groups differ in multivariate isotopic structure.

permanova <- adonis2(
  SIA_matrix ~ Group * Occasion,
  data = SIA_avg,
  method = "euclidean",
  permutations = 999,
  by = "terms"
)

permanova



# attempting at permanova for all occasions

library(vegan)
library(dplyr)

# List of all unique occasions
occasions <- unique(SIA_avg$Occasion)

# Create an empty list to store results
permanova_results <- list()

# Loop through each Occasion
for (occ in occasions) {
  # Subset the data for this occasion
  data_occ <- SIA_avg %>% filter(Occasion == occ)
  
  # Create the isotope matrix (columns d13C and d15N)
  SIA_matrix_occ <- data_occ %>% select(d13C, d15N)
  
  # Optional: scale if needed
  # SIA_matrix_occ <- scale(SIA_matrix_occ)
  
  # Run PERMANOVA using Group as factor
  perm <- adonis2(SIA_matrix_occ ~ Group, data = data_occ, method = "euclidean")
  
  # Save results in the list
  permanova_results[[occ]] <- perm
}

# Check results for one occasion
permanova_results[["MAY_2021"]]

# Or print all results nicely
for (occ in names(permanova_results)) {
  cat("\nPERMANOVA results for", occ, ":\n")
  print(permanova_results[[occ]])
}


