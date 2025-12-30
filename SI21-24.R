# BVR SI Data for Kelley Sinning from May 2021 to present

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
SIA_FISH <- read_csv("Meta_SIA_Data.csv")

Sampling_dates <- read_csv("Sampling_dates.csv")

Sampling_dates <- Sampling_dates %>%
  mutate(Sampling_date = mdy(Sampling_date))


problems()


str(SIA_FISH)

# Putting sampling occassions in order
SIA_FISH$Occasion <- factor(SIA_FISH$Occasion, levels = c("MAY_2021", "AUG_2021", "OCT_2021",
                                                 "MAY_2022", "AUG_2022", "OCT_2022",
                                                 "MAY_2023", "AUG_2023", "OCT_2023",
                                                 "MAY_2024", "AUG_2024", "OCT_2024"))


# Adding a sampling date column so it isn't just month/year 
# better for linking with specific hydrograph patterns
SIA_FISH <- SIA_FISH %>%
  left_join(Sampling_dates, by = c("Occasion", "Location"))





## PLOTTING----------------------------------------------------------------------

# INORGANIC CARBON

# First we want to test whether there is evidence for an effect of inorganic carbon
# on the carbon isotope values. If this was a factor we would expect to see a 
# positive relationship between %C and d13C. Why?

# If you just use SIA_FISH the below graphs are hideous bc the bugs make it really
# hard to differentiate facets, so here I filtered for just BNT and MTS
SIA_JUST_FISH <- SIA_FISH %>%
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

# FORMAT FOR SIBER
data_SIA_SIBER <- SIA_FISH %>%
  select(d13C, d15N, Species, Occasion) %>%
  rename("iso1" = "d13C",
         "iso2" = "d15N",
         "group" = "Species",
         "community" = "Occasion") %>%
  filter(Species == "BNT", Species == "MTS")
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

# CREATE SIBER OBJECT
SIA_SIBER_OBJECT <- createSiberObject(data_SIA_SIBER)

# Putting stuff in order
data_SIA_SIBER$community <- factor(data_SIA_SIBER$community, levels = c("MAY_2021", "AUG_2021", "OCT_2021",
                                                                        "MAY_2022", "AUG_2022", "OCT_2022",
                                                                        "MAY_2023", "AUG_2023", "OCT_2023",
                                                                        "MAY_2024", "AUG_2024", "OCT_2024"))

# SIBER PLOT

ggplot(data_SIA_SIBER, aes(x = iso1, y = iso2, colour = group)) +
  geom_point(size = 1.5, shape = 1) +
  stat_ellipse(position="identity", geom = "polygon", aes(fill = group), level=0.4, linewidth = 1, alpha = 0.25) +
  stat_ellipse(position="identity", level=0.95, linewidth=0.5, linetype = 2) +
  facet_wrap(. ~ community, ncol = 3, nrow = 4) +
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


# PLOT
ggplot(NICHE_WIDTHS, aes(x = MONTH, y = SEA, group = SPECIES)) +
  geom_point(aes(color = SPECIES)) +
  geom_line(aes(color = SPECIES)) +
  facet_wrap(~YEAR) +
  theme_minimal() +
  ylab("SEA") +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_blank()) 

# REMOVE
rm(NICHE_WIDTHS)

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



# work on faceting by month and year









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
# read_waterdata_daily is newer but not working here for some reason.....
site <- "09057500"           # USGS site number (Blue River below Green Mountain)
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
startDate <- "2021-01-01"
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
