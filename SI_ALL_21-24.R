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
problems(SIA_ALL)
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



# Eliminating a bunch of family bug variation by not differentiating by family
SIA_ALL <- SIA_ALL %>% 
  mutate(Species = if_else(Species == "FRY", "Fry", Species)) %>%
  filter(Species != "Fish Eggs") %>%   # remove fish eggs them here
  mutate(
    Group = case_when(
      Species %in% c("BNT", "MTS", "Fry") ~ Species,
      TRUE ~ "Macro"
    )
  ) %>%
  drop_na()

# Keeping family differentiation for SIBER
SIA_ALL_messy <- SIA_ALL %>% 
  mutate(Species = if_else(Species == "FRY", "Fry", Species)) %>%
  filter(Species != "Fish Eggs") %>%   # remove fish eggs them here
  drop_na()


# INORGANIC CARBON-------------------------
SIA_ALL <- SIA_ALL %>%
  mutate(CN = as.numeric(CN))


ggplot(SIA_ALL, aes(x = CN, y = d13C)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = T, color = "#70A494") +
  facet_grid(Group ~ Occasion,  scales = "free_x")+
  scale_x_continuous(n.breaks = 4) +
  theme(
    panel.grid = element_blank()
  )

# LIPID EFFECTS-------------------------

# C:N vs d13C SCATTER PLOT
SIA_ALL <- SIA_ALL %>%
  mutate(CN = as.numeric(CN))


ggplot(SIA_ALL, aes(x = CN, y = d13C)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = T, color = "#70A494") +
  facet_grid(Group ~ Occasion,  scales = "free_x")+
  scale_x_continuous(n.breaks = 4) +
  theme(
    panel.grid = element_blank()
  )

# C:N BOXPLOT
ggplot(SIA_ALL, aes(x=Species, group = Group, y=CN, colour=Group)) +
  geom_boxplot() +
  geom_point(position="jitter", size = 0.75, shape = 1)+
  #ylim(3,5) +
  geom_hline(aes(yintercept=4), colour="#70A494") +
  facet_grid(. ~ Occasion) +
  theme_minimal()

# PROPORTION OF TAXA UNDER 3.5 C:N THRESHOLD
SIA_ALL %>%
  summarize(proportion = sum(CN < 3.5) / n()) %>% 
  kable(caption = "% C:N < 3.5", 
        booktabs= T, 
        digits=2) # .61 under 3.5

# MEAN C:N 
SIA_ALL %>%
  group_by(Group) %>%
  summarise('Mean C:N' = mean(CN)) %>% 
  kable(caption = "MEAN LIPID CONTENT (C:N)", 
        booktabs= T, 
        digits=2) 



# NICHE OVERLAP FOR EACH OCCASION-----------------------------------------------

# FORMAT FOR SIBER

data_SIA_SIBER <- SIA_ALL %>%
  select(d13C, d15N, Group, Occasion) %>%
  rename("iso1" = "d13C",
         "iso2" = "d15N",
         "group" = "Group",
         "community" = "Occasion") %>%
  as.data.frame()


# Prepare data in chronological order
data_SIA_SIBER <- SIA_ALL %>%
  select(d13C, d15N, Group, Occasion) %>%
  rename(
    iso1 = d13C,
    iso2 = d15N,
    group = Group,
    community = Occasion
  ) %>%
  mutate(
    # convert community to factor ordered by date
    community = factor(community, levels = unique(community))
  ) %>%
  as.data.frame()


data_SIA_SIBER <- data_SIA_SIBER %>%
  group_by(community, group) %>%
  filter(n() >= 3) %>%   # minimum for SIBER to work
  ungroup()

str(data_SIA_SIBER)
# Create Siber object
SIA_SIBER_OBJECT <- createSiberObject(data_SIA_SIBER)


# Putting stuff in order
data_SIA_SIBER$community <- factor(data_SIA_SIBER$community, levels = c("MAY_2021", "AUG_2021", "OCT_2021",
                                                                        "MAY_2022", "AUG_2022", "OCT_2022",
                                                                        "MAY_2023", "AUG_2023", "OCT_2023",
                                                                        "MAY_2024", "AUG_2024", "OCT_2024"))

data_SIA_SIBER$community <- factor(data_SIA_SIBER$community, levels = c("MAY_2021", "MAY_2022","MAY_2023","MAY_2024",
                                                                        "AUG_2021", "AUG_2022","AUG_2023","AUG_2024",
                                                                        "OCT_2021", "OCT_2022","OCT_2023","OCT_2024"))


data_SIA_SIBER$group <- factor(data_SIA_SIBER$group, levels = c("BNT", "MTS", "Macro", "Fry"))

# SIBER PLOT

ggplot(data_SIA_SIBER, aes(x = iso1, y = iso2, colour = group)) +
  geom_point(size = 1.5, shape = 1) +
  stat_ellips