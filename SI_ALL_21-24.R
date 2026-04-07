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
  stat_ellipse(position="identity", geom = "polygon", aes(fill = group), level=0.4, linewidth = 1, alpha = 0.25) +
  stat_ellipse(position="identity", level=0.95, linewidth=0.5, linetype = 2) +
  facet_wrap(. ~ community, ncol = 4, nrow = 3, scales = "free") +
  #facet_wrap(. ~ community, ncol = 3, nrow = 4, scales = "free") + Use this for first order option
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  scale_color_manual(values = c("#008080", "#DE8A5A","#70A494", "#B4C8A8")) +
  scale_fill_manual(values = c("#008080", "#DE8A5A", "#70A494", "#B4C8A8")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.99, 0.01),
    legend.justification = c("right", "bottom"),
    legend.text = element_text(size = 8), # smaller legend text
    legend.key.size = unit(0.4, "cm"),    # smaller legend boxes
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA)
  )



# With macro family differentiation THIS LOOKS INSANE LOL-----------------------------
data_SIA_SIBER_messy <- SIA_ALL_messy %>%
  select(d13C, d15N, Species, Occasion) %>%
  rename("iso1" = "d13C",
         "iso2" = "d15N",
         "group" = "Species",
         "community" = "Occasion") %>%
  as.data.frame()


# Prepare data in chronological order
data_SIA_SIBER_messy <- SIA_ALL_messy %>%
  select(d13C, d15N, Species, Occasion) %>%
  rename(
    iso1 = d13C,
    iso2 = d15N,
    group = Species,
    community = Occasion
  ) %>%
  mutate(
    # convert community to factor ordered by date
    community = factor(community, levels = unique(community))
  ) %>%
  as.data.frame()


data_SIA_SIBER_messy <- data_SIA_SIBER_messy %>%
  group_by(community, group) %>%
  filter(n() >= 3) %>%   # minimum for SIBER to work
  ungroup()

str(data_SIA_SIBER_messy)
# Create Siber object
SIA_SIBER_OBJECT <- createSiberObject(data_SIA_SIBER_messy)


# Putting stuff in order
data_SIA_SIBER_messy$community <- factor(data_SIA_SIBER_messy$community, levels = c("MAY_2021", "AUG_2021", "OCT_2021",
                                                                        "MAY_2022", "AUG_2022", "OCT_2022",
                                                                        "MAY_2023", "AUG_2023", "OCT_2023",
                                                                        "MAY_2024", "AUG_2024", "OCT_2024"))

data_SIA_SIBER_messy$community <- factor(data_SIA_SIBER_messy$community, levels = c("MAY_2021", "MAY_2022","MAY_2023","MAY_2024",
                                                                        "AUG_2021", "AUG_2022","AUG_2023","AUG_2024",
                                                                        "OCT_2021", "OCT_2022","OCT_2023","OCT_2024"))


unique(data_SIA_SIBER_messy$group)

# SIBER PLOT


ggplot(data_SIA_SIBER_messy, aes(x = iso1, y = iso2, colour = group)) +
  geom_point(size = 1.5, shape = 1) +
  stat_ellipse(position="identity", geom = "polygon", aes(fill = group), level=0.4, linewidth = 1, alpha = 0.25) +
  stat_ellipse(position="identity", level=0.95, linewidth=0.5, linetype = 2) +
  facet_wrap(. ~ community, ncol = 4, nrow = 3, scales = "free") +
  #facet_wrap(. ~ community, ncol = 3, nrow = 4, scales = "free") + Use this for first order option
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
 # scale_color_manual(values = rand_cols) +
  #scale_fill_manual(values = rand_cols) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.99, 0.01),
    legend.justification = c("right", "bottom"),
    legend.text = element_text(size = 8), # smaller legend text
    legend.key.size = unit(0.4, "cm"),    # smaller legend boxes
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA)
  )






# SUM STATS FOR EACH GROUP: TA, SEA & SEAc ----------------------------------
# TA = Total Area = convex hull around all the data 
# SEA = Standard Ellipse Area = 40% of the data 
# SEAc = Standard Ellipse Area Corrected (small sample size corrected [< 20-30 individuals])
group.ML <- groupMetricsML(SIA_SIBER_OBJECT)

NICHE_WIDTHS <-group.ML %>%
  as.data.frame() %>%
  tibble::rownames_to_column("NICHE_AREA") %>%
  pivot_longer(!NICHE_AREA, names_to = "GROUP", values_to = "NICHE_AREA_VALUES") %>%
  separate(GROUP, into = c("COMMUNITY", "SPECIES"), sep = "\\.") %>%
  separate(COMMUNITY, into = c("MONTH", "YEAR"), sep = "\\_") %>%
  select(YEAR, MONTH, SPECIES, NICHE_AREA, NICHE_AREA_VALUES) %>%
  pivot_wider(names_from = NICHE_AREA, values_from = NICHE_AREA_VALUES) %>%
  mutate(
    YEAR = factor(YEAR, levels = c("2021", "2022", "2023", "2024")),
    MONTH = factor(MONTH, levels = c("MAY", "AUG", "OCT")))

kable(NICHE_WIDTHS, caption = "ISOTOPIC NICHE AREAS", digits = 1, booktabs = T) 


# PLOT NICHE WIDTHS
# This is SEA but if you change to SEAc it's basically the same result
ggplot(NICHE_WIDTHS, aes(x = MONTH, y = SEA, group = SPECIES)) +
  geom_point(aes(color = SPECIES)) +
  geom_line(aes(color = SPECIES)) +
  facet_wrap(~YEAR) +
  theme_minimal() +
  ylab("SEA") +
  scale_color_manual(values = c("#008080", "#DE8A5A","#70A494", "#B4C8A8")) +
  scale_fill_manual(values = c("#008080", "#DE8A5A","#70A494", "#B4C8A8")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_blank()) 

rm(NICHE_WIDTHS) # remove



# Out of the total isotopic space occupied by these two species, around about 
# __ is overlapping
# Setting p.interval = 0.40 will result in an ellipse that contains approximately 40% of the data (SEA).
# Setting p.interval = 0.95 will result in an ellipse that contains approximately 95% of the data.

######### MAY 2021 ######### 
MAY_2021_AREAS <- maxLikOverlap("MAY_2021.BNT", "MAY_2021.MTS", SIA_SIBER_OBJECT, 
                                p.interval = 0.40, n = 100)
# Overlap relative to the total area used by that population -- so in this case, 
# the overlap as a proportion of the non-overlapping area of the two ellipses, would be:
MAY_2021_OVERLAP <- MAY_2021_AREAS[3] / (MAY_2021_AREAS[2] + 
                                           MAY_2021_AREAS[1] -
                                           MAY_2021_AREAS[3])

######### AUG 2021 ######### 
AUG_2021_AREAS <- maxLikOverlap("AUG_2021.BNT", "AUG_2021.MTS", SIA_SIBER_OBJECT, 
                                p.interval = 0.40, n = 100)
AUG_2021_OVERLAP <- AUG_2021_AREAS[3] / (AUG_2021_AREAS[2] + 
                                           AUG_2021_AREAS[1] -
                                           AUG_2021_AREAS[3])

######### OCT 2021 ######### 

OCT_2021_AREAS <- maxLikOverlap("OCT_2021.BNT", "OCT_2021.MTS", SIA_SIBER_OBJECT, 
                                p.interval = 0.40, n = 100)
OCT_2021_OVERLAP <- OCT_2021_AREAS[3] / (OCT_2021_AREAS[2] + 
                                           OCT_2021_AREAS[1] -
                                           OCT_2021_AREAS[3])
######### MAY 2022 ######### 
MAY_2022_AREAS <- maxLikOverlap("MAY_2022.BNT", "MAY_2022.MTS", SIA_SIBER_OBJECT, 
                                p.interval = 0.40, n = 100)
MAY_2022_OVERLAP <- MAY_2022_AREAS[3] / (MAY_2022_AREAS[2] + 
                                           MAY_2022_AREAS[1] -
                                           MAY_2022_AREAS[3])
######### AUGUST 2022 ######### 
AUG_2022_AREAS <- maxLikOverlap("AUG_2022.BNT", "AUG_2022.MTS", SIA_SIBER_OBJECT, 
                                p.interval = 0.40, n = 100)
AUG_2022_OVERLAP <- AUG_2022_AREAS[3] / (AUG_2022_AREAS[2] + 
                                           AUG_2022_AREAS[1] -
                                           AUG_2022_AREAS[3])
######### OCTOBER 2022 ######### 
OCT_2022_AREAS <- maxLikOverlap("OCT_2022.BNT", "OCT_2022.MTS", SIA_SIBER_OBJECT, 
                                p.interval = 0.40, n = 100)
OCT_2022_OVERLAP <- OCT_2022_AREAS[3] / (OCT_2022_AREAS[2] + 
                                           OCT_2022_AREAS[1] -
                                           OCT_2022_AREAS[3])

######### MAY 2023 ######### 
MAY_2023_AREAS <- maxLikOverlap("MAY_2023.BNT", "MAY_2023.MTS", SIA_SIBER_OBJECT, 
                                p.interval = 0.40, n = 100)
MAY_2023_OVERLAP <- MAY_2023_AREAS[3] / (MAY_2023_AREAS[2] + 
                                           MAY_2023_AREAS[1] -
                                           MAY_2023_AREAS[3])

######### AUG 2023 ######### 
AUG_2023_AREAS <- maxLikOverlap("AUG_2023.BNT", "AUG_2023.MTS", SIA_SIBER_OBJECT, 
                                p.interval = 0.40, n = 100)
AUG_2023_OVERLAP <- AUG_2023_AREAS[3] / (AUG_2023_AREAS[2] + 
                                           AUG_2023_AREAS[1] -
                                           AUG_2023_AREAS[3])

######### OCT 2023 ######### 

OCT_2023_AREAS <- maxLikOverlap("OCT_2023.BNT", "OCT_2023.MTS", SIA_SIBER_OBJECT, 
                                p.interval = 0.40, n = 100)
OCT_2023_OVERLAP <- OCT_2023_AREAS[3] / (OCT_2023_AREAS[2] + 
                                           OCT_2023_AREAS[1] -
                                           OCT_2023_AREAS[3])
######### MAY 2024 ######### 
MAY_2024_AREAS <- maxLikOverlap("MAY_2024.BNT", "MAY_2024.MTS", SIA_SIBER_OBJECT, 
                                p.interval = 0.40, n = 100)
MAY_2024_OVERLAP <- MAY_2024_AREAS[3] / (MAY_2024_AREAS[2] + 
                                           MAY_2024_AREAS[1] -
                                           MAY_2024_AREAS[3])
######### AUGUST 2024 ######### 
AUG_2024_AREAS <- maxLikOverlap("AUG_2024.BNT", "AUG_2024.MTS", SIA_SIBER_OBJECT, 
                                p.interval = 0.40, n = 100)
AUG_2024_OVERLAP <- AUG_2024_AREAS[3] / (AUG_2024_AREAS[2] + 
                                           AUG_2024_AREAS[1] -
                                           AUG_2024_AREAS[3])
######### OCTOBER 2024 ######### 
OCT_2024_AREAS <- maxLikOverlap("OCT_2024.BNT", "OCT_2024.MTS", SIA_SIBER_OBJECT, 
                                p.interval = 0.40, n = 100)
OCT_2024_OVERLAP <- OCT_2024_AREAS[3] / (OCT_2024_AREAS[2] + 
                                           OCT_2024_AREAS[1] -
                                           OCT_2024_AREAS[3])



# PRINT RESULTS
MAY_2021_OVERLAP <- MAY_2021_OVERLAP %>%
  as.data.frame() %>%
  mutate(YEAR = 2021,
         SEASON = "MAY") %>%
  rename("OVERLAP" = ".") %>%
  rownames_to_column(var = "RowNumber") %>%
  select(-RowNumber)

AUG_2021_OVERLAP <- AUG_2021_OVERLAP %>%
  as.data.frame() %>%
  mutate(YEAR = 2021,
         SEASON = "AUGUST") %>%
  rename("OVERLAP" = ".") %>%
  rownames_to_column(var = "RowNumber") %>%
  select(-RowNumber)

OCT_2021_OVERLAP <- OCT_2021_OVERLAP %>%
  as.data.frame() %>%
  mutate(YEAR = 2021,
         SEASON = "OCTOBER") %>%
  rename("OVERLAP" = ".") %>%
  rownames_to_column(var = "RowNumber") %>%
  select(-RowNumber)

MAY_2022_OVERLAP <- MAY_2022_OVERLAP %>%
  as.data.frame() %>%
  mutate(YEAR = 2022,
         SEASON = "MAY") %>%
  rename("OVERLAP" = ".") %>%
  rownames_to_column(var = "RowNumber") %>%
  select(-RowNumber)

AUG_2022_OVERLAP <- AUG_2022_OVERLAP %>%
  as.data.frame() %>%
  mutate(YEAR = 2022,
         SEASON = "AUGUST") %>%
  rename("OVERLAP" = ".") %>%
  rownames_to_column(var = "RowNumber") %>%
  select(-RowNumber) 

OCT_2022_OVERLAP <- OCT_2022_OVERLAP %>%
  as.data.frame() %>%
  mutate(YEAR = 2022,
         SEASON = "OCTOBER") %>%
  rename("OVERLAP" = ".") %>%
  rownames_to_column(var = "RowNumber") %>%
  select(-RowNumber)

MAY_2023_OVERLAP <- MAY_2023_OVERLAP %>%
  as.data.frame() %>%
  mutate(YEAR = 2023,
         SEASON = "MAY") %>%
  rename("OVERLAP" = ".") %>%
  rownames_to_column(var = "RowNumber") %>%
  select(-RowNumber)

AUG_2023_OVERLAP <- AUG_2023_OVERLAP %>%
  as.data.frame() %>%
  mutate(YEAR = 2023,
         SEASON = "AUGUST") %>%
  rename("OVERLAP" = ".") %>%
  rownames_to_column(var = "RowNumber") %>%
  select(-RowNumber)

OCT_2023_OVERLAP <- OCT_2023_OVERLAP %>%
  as.data.frame() %>%
  mutate(YEAR = 2023,
         SEASON = "OCTOBER") %>%
  rename("OVERLAP" = ".") %>%
  rownames_to_column(var = "RowNumber") %>%
  select(-RowNumber)

MAY_2024_OVERLAP <- MAY_2024_OVERLAP %>%
  as.data.frame() %>%
  mutate(YEAR = 2024,
         SEASON = "MAY") %>%
  rename("OVERLAP" = ".") %>%
  rownames_to_column(var = "RowNumber") %>%
  select(-RowNumber)

AUG_2024_OVERLAP <- AUG_2024_OVERLAP %>%
  as.data.frame() %>%
  mutate(YEAR = 2024,
         SEASON = "AUGUST") %>%
  rename("OVERLAP" = ".") %>%
  rownames_to_column(var = "RowNumber") %>%
  select(-RowNumber) 

OCT_2024_OVERLAP <- OCT_2024_OVERLAP %>%
  as.data.frame() %>%
  mutate(YEAR = 2024,
         SEASON = "OCTOBER") %>%
  rename("OVERLAP" = ".") %>%
  rownames_to_column(var = "RowNumber") %>%
  select(-RowNumber)

OVERLAP_DF <- rbind(MAY_2021_OVERLAP, AUG_2021_OVERLAP, OCT_2021_OVERLAP, 
                    MAY_2022_OVERLAP, AUG_2022_OVERLAP, OCT_2022_OVERLAP,
                    MAY_2023_OVERLAP, AUG_2023_OVERLAP, OCT_2023_OVERLAP, 
                    MAY_2024_OVERLAP, AUG_2024_OVERLAP, OCT_2024_OVERLAP) %>%
  select(YEAR, SEASON, OVERLAP) %>%
  mutate(
    YEAR = factor(YEAR, levels = c("2021", "2022", "2023", "2024")),
    SEASON = factor(SEASON, levels = c("MAY", "AUGUST", "OCTOBER")))

# TABLE
kable(OVERLAP_DF, caption = "ISOTOPIC NICHE OVERLAP (SEA)", digits = 2, booktabs = T) 

# PLOT
ggplot(OVERLAP_DF, aes(x = SEASON, y = OVERLAP, group = YEAR)) +
  geom_point(aes(color = YEAR)) +
  geom_line(aes(color = YEAR)) +
  theme_minimal() +
  ylab("Isotpoic Niche Overlap (SEA)") +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_blank()) 





# NMDS----------------------------------------------------------------------------
head(SIA_ALL)


SIA_ALL <- SIA_ALL %>%
  select(d13C, d15N, Species, Occasion, Location) %>%
  drop_na()

str(SIA_ALL)

SIA_ALL_clean <- SIA_ALL %>% # Eliminating a bunch of family bug variation in permanova
  select(Occasion, Location, Species, d13C, d15N) %>%
  mutate(Species = if_else(Species == "FRY", "Fry", Species)) #%>%
  #mutate(
   # Group = case_when(
    #  Species %in% c("BNT", "MTS", "Fry", "Fish Eggs") ~ Species,
     # TRUE ~ "Macro"
    #)
  #) %>%
  #drop_na()


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


