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
                                                          "MAY_2024", "AUG_2024", "OCT_2024"))


# Adding a sampling date column so it isn't just month/year 
# better for linking with specific hydrograph patterns
SIA_ALL <- SIA_ALL %>%
  left_join(Sampling_dates, by = c("Occasion", "Location"))

# NMDS----------------------------------------------------------------------------
head(SIA_ALL)

SIA_ALL <- SIA_ALL %>%
  select(d13C, d15N, Species) %>%
  drop_na()


SIA_matrix <- SIA_ALL %>%
  select(d13C, d15N)  # add more isotopes here if you have them
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
  bind_cols(SIA_ALL %>% select(Species))

ggplot(scores_nmds, aes(x = NMDS1, y = NMDS2, color = Taxon)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = Taxon),
               geom = "polygon",
               alpha = 0.2,
               color = NA) +
  theme_classic() +
  labs(
    title = "NMDS of Stable Isotope Space",
    subtitle = paste("Stress =", round(nmds_iso$stress, 3)),
    x = "NMDS1",
    y = "NMDS2"
  )

