# BVR SI Fish Data for Kelley Sinning from May 2021 to present
# Focusing on isotopic niche overlap between MTS and BNT

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

setwd("~/Library/CloudStorage/OneDrive-Colostate/Data/BVR")

# IMPORT SIA DATA
SIA_ALL <- read_csv("Meta_SIA_Data.csv")

Sampling_dates <- read_csv("Sampling_dates.csv")

Sampling_dates <- Sampling_dates %>%
  mutate(Sampling_date = mdy(Sampling_date))


problems()


str(SIA_ALL)

# Putting sampling occassions in order
SIA_ALL$Occasion <- factor(SIA_ALL$Occasion, levels = c("MAY_2021", "AUG_2021", "OCT_2021",
                                                 "MAY_2022", "AUG_2022", "OCT_2022",
                                                 "MAY_2023", "AUG_2023", "OCT_2023",
                                                 "MAY_2024", "AUG_2024", "OCT_2024"))


# Adding a sampling date column so it isn't just month/year 
# better for linking with specific hydrograph patterns
SIA_ALL <- SIA_ALL %>%
  left_join(Sampling_dates, by = c("Occasion", "Location"))





## PLOTTING----------------------------------------------------------------------

# INORGANIC CARBON

# First we want to test whether there is evidence for an effect of inorganic carbon
# on the carbon isotope values. If this was a factor we would expect to see a 
# positive relationship between %C and d13C. Why?

# If you just use SIA_FISH the below graphs are hideous bc the bugs make it really
# hard to differentiate facets, so here I filtered for just BNT and MTS
SIA_JUST_FISH <- SIA_ALL %>%
  filter(Species %in% c("BNT", "MTS"))

SIA_JUST_FISH$Occasion <- factor(SIA_JUST_FISH$Occasion, levels = c("MAY_2021", "AUG_2021", "OCT_2021",
                                                                    "MAY_2022", "AUG_2022", "OCT_2022",
                                                                    "MAY_2023", "AUG_2023", "OCT_2023",
                                                                    "MAY_2024", "AUG_2024", "OCT_2024"))


ggplot(SIA_JUST_FISH, aes(x=C, y=d13C)) +
  geom_point() +
  geom_smooth(method = "lm", se=T,formula=y~x) +
  #facet_grid(Species ~ Month, scales = "free_x") 
  facet_grid(Species ~  Occasion, scales = "free_x") 
  #facet_grid(Species ~  Year, scales = "free_x")



# LIPID EFFECTS

# Lipids contain carbon but no nitrogen AND are 13C depleted relative to protein 
# (i.e. muscle). So, lipid rich samples will be depleted in 13C (more negative delta 13C), in extreme cases 
#this may bias our interpretation of their position in the food web.

# Alternatively we can identify whether any individuals / populations have C:N 
# values above the 3.5 threshold. If C:N value are mostly below 4, it is unlikely 
# that lipids will influence the d13C values of most species.

# C:N vs d13C SCATTER PLOT
SIA_JUST_FISH <- SIA_JUST_FISH %>%
  mutate(CN = as.numeric(CN))


ggplot(SIA_JUST_FISH, aes(x = CN, y = d13C)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = T, color = "#70A494") +
  facet_grid(Species ~ Occasion,  scales = "free_x")+
  scale_x_continuous(n.breaks = 4) +
  theme(
    panel.grid = element_blank()
  )


# C:N BOXPLOT
ggplot(SIA_JUST_FISH, aes(x=Species, group = Species, y=CN, colour=Species)) +
  geom_boxplot() +
  geom_point(position="jitter", size = 0.75, shape = 1)+
  #ylim(3,5) +
  geom_hline(aes(yintercept=4), colour="#70A494") +
  facet_grid(. ~ Occasion) +
  theme_minimal()

# PROPORTION OF FISH UNDER 3.5 C:N THRESHOLD
SIA_JUST_FISH %>%
  summarize(proportion = sum(CN < 3.5) / n()) %>% 
  kable(caption = "% C:N < 3.5", 
        booktabs= T, 
        digits=2) # I think we are in the safe

# MEAN C:N FOR BROWN TROUT & MOTTLED SCULPIN
SIA_FISH %>%
  group_by(Species) %>%
  summarise('Mean C:N' = mean(CN)) %>% 
  kable(caption = "MEAN LIPID CONTENT (C:N)", 
        booktabs= T, 
        digits=2) 





  
  
  
# NICHE OVERLAP FOR EACH OCCASION-----------------------------------------------

# FORMAT FOR SIBERx

data_SIA_SIBER <- SIA_FISH %>%
  select(d13C, d15N, Species, Occasion) %>%
  rename("iso1" = "d13C",
         "iso2" = "d15N",
         "group" = "Species",
         "community" = "Occasion") %>%
  filter(group %in% c("BNT", "MTS")) %>%
as.data.frame()

# Prepare data in chronological order
data_SIA_SIBER <- SIA_FISH %>%
  select(d13C, d15N, Species, Occasion) %>%
  rename(
    iso1 = d13C,
    iso2 = d15N,
    group = Species,
    community = Occasion
  ) %>%
  filter(group %in% c("BNT", "MTS")) %>%
  mutate(
    # convert community to factor ordered by date
    community = factor(community, levels = unique(community))
  ) %>%
  as.data.frame()


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





# SIBER PLOT

ggplot(data_SIA_SIBER, aes(x = iso1, y = iso2, colour = group)) +
  geom_point(size = 1.5, shape = 1) +
  stat_ellipse(position="identity", geom = "polygon", aes(fill = group), level=0.4, linewidth = 1, alpha = 0.25) +
  stat_ellipse(position="identity", level=0.95, linewidth=0.5, linetype = 2) +
  facet_wrap(. ~ community, ncol = 4, nrow = 3, scales = "free") +
  #facet_wrap(. ~ community, ncol = 3, nrow = 4, scales = "free") + Use this for first order option
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  scale_color_manual(values = c("#008080", "#DE8A5A")) +
  scale_fill_manual(values = c("#008080", "#DE8A5A")) +
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


# SUM STATS FOR EACH GROUP: TA, SEA & SEAc -----------
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
  scale_color_manual(values = c("#008080", "#DE8A5A")) +
  scale_fill_manual(values = c("#008080", "#DE8A5A")) +
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







# NICHE OVERLAP FOR EACH YEAR-----------------------------------------------


# Prepare data in chronological order
data_SIA_SIBER <- SIA_FISH %>%
  select(d13C, d15N, Species, Year) %>%
  rename(
    iso1 = d13C,
    iso2 = d15N,
    group = Species,
    community = Year
  ) %>%
  filter(group %in% c("BNT", "MTS")) %>%
  mutate(
    # convert community to factor ordered by date
    community = factor(community, levels = unique(community))
  ) %>%
  as.data.frame()

# CREATE SIBER OBJECT
SIA_SIBER_OBJECT <- createSiberObject(data_SIA_SIBER)

# SIBER PLOT
ggplot(data_SIA_SIBER, aes(x = iso1, y = iso2, colour = group)) +
  geom_point(size = 1.5, shape = 1) +
  stat_ellipse(position="identity", geom = "polygon", aes(fill = group), level=0.4, linewidth = 1, alpha = 0.25) +
  stat_ellipse(position="identity", level=0.95, linewidth=0.5, linetype = 2) +
  facet_wrap(. ~ community) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom")
  )



# NICHE OVERLAP FOR EACH SITE-----------------------------------------------

# Prepare data in chronological order
data_SIA_SIBER <- SIA_FISH %>%
  select(d13C, d15N, Species, Location) %>%
  rename(
    iso1 = d13C,
    iso2 = d15N,
    group = Species,
    community = Location
  ) %>%
  filter(group %in% c("BNT", "MTS")) %>%
  mutate(
    # convert community to factor ordered by date
    community = factor(community, levels = unique(community))
  ) %>%
  as.data.frame()

# CREATE SIBER OBJECT
SIA_SIBER_OBJECT <- createSiberObject(data_SIA_SIBER)

# SIBER PLOT
ggplot(data_SIA_SIBER, aes(x = iso1, y = iso2, colour = group)) +
  geom_point(size = 1.5, shape = 1) +
  stat_ellipse(position="identity", geom = "polygon", aes(fill = group), level=0.4, linewidth = 1, alpha = 0.25) +
  stat_ellipse(position="identity", level=0.95, linewidth=0.5, linetype = 2) +
  facet_wrap(. ~ community) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  scale_color_manual(values = c("#E69F00", "#0072B2")) +
  scale_fill_manual(values = c("#E69F00", "#0072B2")) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom")
  )


### PLAYING AROUND
NICHE_WIDTHS_BNT <- NICHE_WIDTHS %>% 
  filter(SPECIES == "BNT")

NICHE_WIDTHS_MTS <- NICHE_WIDTHS %>% 
  filter(SPECIES == "MTS")



model <- lmer(SEA ~ MONTH + (1|YEAR),
              data = NICHE_WIDTHS_BNT) 
summary(model)

model <- lmer(SEA ~ MONTH + (1|YEAR),
              data = NICHE_WIDTHS_MTS) 
summary(model)



model <- lm(SEA ~ SPECIES,
            data = NICHE_WIDTHS) 
summary(model)  


anova_model <- aov(SEA ~ MONTH, data = NICHE_WIDTHS_BNT)
summary(anova_model)
TukeyHSD(anova_model)


anova_model <- aov(SEA ~ MONTH, data = NICHE_WIDTHS)
summary(anova_model)
TukeyHSD(anova_model)

### HYDROGRAPH------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(tidyverse)
library(purrr)


# RETRIEVING USGS GAUGE DATA
install.packages("dataRetrieval", type = "source")
library(dataRetrieval)
packageVersion("dataRetrieval")  # should be >= 2.7.19


# https://doi-usgs.github.io/dataRetrieval/reference/read_waterdata_daily.html
# read_waterdata_daily is newer 
site <- "USGS-09057500"           # USGS site number (Blue River below Green Mountain)
parameter_code <- "00060"    # Parameter code for discharge (ft³/s)
statistic_id <- "00003"      # Statistic code for daily mean
start_date <- "2021-01-01"
end_date <- "2024-12-29"

# Retrieve daily discharge data
discharge_data <- read_waterdata_daily(
  monitoring_location_id = site,
  parameter_code = parameter_code,
  statistic_id = statistic_id,
  time = c(start_date, end_date)
)


# So using readNWISdv...it's older but works
site <- "09057500"          # USGS site number
parameterCd <- "00060"      # Discharge (cfs)
statCd <- "00003"           # Daily mean
startDate <- "2020-01-01"
endDate <- "2024-12-29"

# Retrieve daily values
discharge_data <- readNWISdv(
  siteNumbers = site,
  parameterCd = parameterCd,
  statCd = statCd,
  startDate = startDate,
  endDate = endDate
)

# Renaming columns

discharge_data <- discharge_data %>%
  rename(Discharge_cfs = X_00060_00003,
         Qualifier = X_00060_00003_cd
  )
# A, P = data qualifiers (approved, provisional) for Qualifier


ggplot(discharge_data, aes(x = Date, y = Discharge_cfs)) +
  geom_line() +
  geom_vline(
    data = Sampling_dates,
    aes(xintercept = Sampling_date),
    color = "#70A494",
    linetype = "solid"
  ) +
  labs(
    x = "Date",
    y = "Discharge (cfs)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

class(discharge_data$Date)
class(Sampling_dates$Sampling_date)


# Function to compute 30-day stats for each sample, this works
get_30day_stats <- function(Sampling_date, discharge_data) {
  start_date <- Sampling_date - days(30)
  subset <- discharge_data %>%
    filter(Date > start_date & Date <= Sampling_date)
  
  mean_cfs <- mean(subset$Discharge_cfs, na.rm = TRUE)
  cv_cfs <- sd(subset$Discharge_cfs, na.rm = TRUE) / mean_cfs
  
  return(c(mean_30d = mean_cfs, cv_30d = cv_cfs))
}

# Apply function to each sample
stats_matrix <- t(sapply(Sampling_dates$Sampling_date, get_30day_stats, discharge_data = discharge_data))

# Combine with Sampling_dates of each SI occasion so we know the mean 30d discharge prior to sampling
Sampling_datess <- cbind(Sampling_dates, stats_matrix)



# IS DISCHARGE RELATED TO ISOTOPE STUFF?---------------------------------------

#discharge_data <- discharge_data %>%
 # rename(Sampling_date = Date)


#SIA_FISH_DISCHARGE <- SIA_FISH %>%
 # left_join(discharge_data %>% select("Sampling_date", "Discharge_cfs"), by = c("Sampling_date"))

# Putting the df in chronological order
NICHE_WIDTHS <- NICHE_WIDTHS %>%
  mutate(
    MONTH_NUM = match(MONTH, c("MAY", "AUG", "OCT")),   # convert month to number
    YEAR = as.numeric(as.character(YEAR))               # ensure YEAR is numeric
  ) %>%
  arrange(YEAR, MONTH_NUM) %>%                          # sort by year then month
  select(-MONTH_NUM)                                    # drop temporary column


# Adding a column for sample dates so it can be left joined with mean_30d discharge prior to first
# sampling day of each occasion to see if there's a relationship with niche widths and flow
NICHE_WIDTHS <- NICHE_WIDTHS %>%
  mutate(Sampling_date = (c("2021-05-24", "2021-05-24","2021-08-02","2021-08-02", "2021-10-10","2021-10-10",
         "2022-05-23", "2022-05-23", "2022-08-01", "2022-08-01", "2022-10-10", "2022-10-10", "2023-05-15",
         "2023-05-15", "2023-08-07", "2023-08-07", "2023-10-09", "2023-10-09", "2024-05-13", "2024-05-13",
         "2024-07-29", "2024-07-29", "2024-10-07", "2024-10-07")))

NICHE_WIDTHS_BNT <- NICHE_WIDTHS_BNT %>%
  mutate(Sampling_date = (c( "2021-05-24","2021-08-02", "2021-10-10",
                             "2022-05-23", "2022-08-01", "2022-10-10", 
                            "2023-05-15", "2023-08-07", "2023-10-09", 
                            "2024-05-13", "2024-07-29", "2024-10-07")))

NICHE_WIDTHS_MTS <- NICHE_WIDTHS_MTS %>%
  mutate(Sampling_date = (c( "2021-05-24","2021-08-02", "2021-10-10",
                             "2022-05-23", "2022-08-01", "2022-10-10", 
                             "2023-05-15", "2023-08-07", "2023-10-09", 
                             "2024-05-13", "2024-07-29", "2024-10-07")))
                          
                          
NICHE_WIDTHS$Sampling_date   <- as.Date(NICHE_WIDTHS$Sampling_date)
NICHE_WIDTHS_BNT$Sampling_date   <- as.Date(NICHE_WIDTHS_BNT$Sampling_date)
NICHE_WIDTHS_MTS$Sampling_date   <- as.Date(NICHE_WIDTHS_MTS$Sampling_date)
Sampling_datess$Sampling_date <- as.Date(Sampling_datess$Sampling_date)


NICHE_WIDTHS_DISCHARGE <- NICHE_WIDTHS %>%
  left_join(Sampling_datess %>% select("Sampling_date", "mean_30d"), by = c("Sampling_date"))

NICHE_WIDTHS_BNT_DISCHARGE <- NICHE_WIDTHS_BNT %>%
  left_join(Sampling_datess %>% select("Sampling_date", "mean_30d"), by = c("Sampling_date"))

NICHE_WIDTHS_MTS_DISCHARGE <- NICHE_WIDTHS_MTS %>%
  left_join(Sampling_datess %>% select("Sampling_date", "mean_30d"), by = c("Sampling_date"))

model <- lmer(SEA ~ mean_30d + (1|MONTH) + (1|YEAR),
              data = NICHE_WIDTHS_DISCHARGE) 
summary(model)


model <- lm(SEA ~ mean_30d, data = NICHE_WIDTHS_DISCHARGE)
summary(model)


model <- lm(SEA ~ mean_30d, data = NICHE_WIDTHS_BNT_DISCHARGE)
summary(model)


model <- lm(SEA ~ mean_30d, data = NICHE_WIDTHS_MTS_DISCHARGE)
summary(model)


# It seems like the 30 days leading up to the first day of sampling was not a 
# significant factor in niche widths of BNT and MTS


# How about overlap?
# Putting the df in chronological order
OVERLAP_DF <- OVERLAP_DF %>%
  mutate(
    MONTH_NUM = match(SEASON, c("MAY", "AUG", "OCT")),   # convert month to number
    YEAR = as.numeric(as.character(YEAR))               # ensure YEAR is numeric
  ) %>%
  arrange(YEAR, MONTH_NUM) %>%                          # sort by year then month
  select(-MONTH_NUM)                                    # drop temporary column


OVERLAP_DF <- OVERLAP_DF %>%
  mutate(Sampling_date = (c( "2021-05-24","2021-08-02", "2021-10-10",
                             "2022-05-23", "2022-08-01", "2022-10-10", 
                             "2023-05-15", "2023-08-07", "2023-10-09", 
                             "2024-05-13", "2024-07-29", "2024-10-07")))

OVERLAP_DF$Sampling_date   <- as.Date(OVERLAP_DF$Sampling_date)

OVERLAP_DF_DISCHARGE <- OVERLAP_DF %>%
  left_join(Sampling_datess %>% select("Sampling_date", "mean_30d"), by = c("Sampling_date"))

model <- lm(OVERLAP ~ mean_30d, data = OVERLAP_DF_DISCHARGE)
summary(model)

# nope

library(readr)


# for sliding window
write_csv(discharge_data, "SI_discharge.csv")
write_csv(NICHE_WIDTHS, "NICHE_WIDTHS_21to24.csv")
write_csv(OVERLAP_DF, "OVERLAP.csv")
