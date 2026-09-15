# Kelley Sinning 9/15/2026
# Cleaning up all the manjor datasheets so they can be used in an analysis together!

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
         -Measurement_mm,-Extra_counts, -Empty_case_mm, -PROOFED, -USE_SAMPLE, -Total_Sampled, -X)
