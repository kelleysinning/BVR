# Code for Kelley Sinning Ph.D.

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
site <- "USGS-09057500"  # USGS site number (Blue River below Green Mountain)
parameter_code <- "00060"    # Parameter code for discharge (ft³/s)
statistic_id <- "00003"      # Statistic code for daily mean
start_date <- "2021-04-01" # Date before first sampling event ever
end_date <- "2026-08-24" # Date following most recent benthotorching/sampling event

# Retrieve daily discharge data
discharge_data <- read_waterdata_daily(
  monitoring_location_id = site,
  parameter_code = parameter_code,
  statistic_id = statistic_id,
  time = c(start_date, end_date),
  skipGeometry = TRUE,
)

discharge_data <- discharge_data %>%
  select(time, value, approval_status)

discharge_data <- discharge_data %>%
  rename(Discharge_cfs = value,
         Date = time
  )


# So using readNWISdv...it's older but works
# More results than read_waterdata_dailybecause it includes provisional data
site <- "09057500"          # USGS site number
parameterCd <- "00060"      # Discharge (cfs)
statCd <- "00003"           # Daily mean
startDate <- "2021-04-01"
endDate <- "2026-08-024"

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

# BRINGING IN ALGAE DATA-------------------------------------------------------
# to merge with discharge

setwd("~/Library/CloudStorage/OneDrive-TheUniversityofMontana/Data/BVR")
didymo_benthotorch <- read.csv("ALL_Bentho_Core.csv")
didymo_benthotorch <- didymo_benthotorch %>%
  filter(!is.na(Sampling_date)) # removing NA columns that arose from comments in the csv

str(didymo_benthotorch)
# Making sure date formats are the same
didymo_benthotorch$Sampling_date <- as.Date(didymo_benthotorch$Sampling_date)  
discharge_data$time <- as.Date(discharge_data$time)

# Making sure columns are numeric
str(didymo_benthotorch)
didymo_benthotorch <- didymo_benthotorch %>%
  mutate(across(c(Cyano, Green, Diatoms, Chlorophyll.A, Velocity),
                as.numeric))

# HYDROGRAPH WITH DATES OF BENTHOTORCH AND VELOCITIES---------------------------------------
library(dplyr)

sort(unique(didymo_benthotorch$Sampling_date))
sample_dates <- c( # Dates that Benthotorch have been taken
  "2023-01-15", "2023-03-21", "2023-05-12", "2023-05-24", "2023-06-27", "2023-07-31",
  "2023-10-04", "2024-05-13", "2024-05-14", "2024-05-15", "2024-05-16", "2024-05-17",
  "2024-07-24", "2024-07-29", "2024-07-30", "2024-07-31", "2024-08-01", "2024-08-02",
  "2024-08-13", "2024-09-26", "2024-10-26", "2024-11-14", "2025-01-16", "2025-03-11",
  "2025-05-19", "2025-05-20", "2025-05-21", "2025-05-22", "2025-05-23", "2025-06-10",
  "2025-06-18", "2025-07-30", "2025-08-27", "2025-09-23", "2025-10-20", "2025-10-22", 
  "2025-10-21", "2025-10-19", "2025-11-14", "2025-12-15", "2026-02-14", "2026-03-07", 
  "2026-04-11", "2026-05-25", "2026-05-26", "2026-05-27", "2026-05-28", "2026-05-29",
  "2026-08-03", "2026-08-04", "2026-08-05", "2026-08-06", "2026-08-07"
)

velocity_dates <- c("2024-07-24", "2024-08-13", "2024-09-26", "2024-10-26", "2024-11-14", "2025-01-16",
                    "2025-03-11", "2025-05-19", "2025-05-20", "2025-05-21", "2025-05-22", "2025-05-23",
                    "2025-06-10", "2025-07-30", "2025-08-27", "2025-09-23", "2025-10-20","2025-10-22", 
                    "2025-10-21", "2025-10-19", "2025-11-14", "2025-12-15", "2026-02-14", "2026-03-07", 
                    "2026-04-11", "2026-05-25", "2026-05-26", "2026-05-27", "2026-05-28", "2026-05-29",
                    "2026-08-03", "2026-08-04", "2026-08-05", "2026-08-06", "2026-08-07"
                    ) # These additional dates have paired velocity



# Convert to Date
sample_dates <- as.Date(sample_dates)
velocity_dates <- as.Date(velocity_dates)


# Adding a new column for it it was sampled or not
discharge_data <- discharge_data %>%
  mutate(sampled = if_else(Date %in% sample_dates, "yes", "no"),
         velocity = if_else(Date %in% velocity_dates, "yes", "no"))


# This filters for all benthotorch data
discharge_data <- discharge_data %>%
  filter(Date > as.Date("2023-01-01") &
           Date < as.Date("2026-08-24"))


# Filter just the sampling dates for vertical lines
sample_dates <- discharge_data %>%
  filter(sampled == "yes") %>%
  pull(Date)

velocity_dates <- discharge_data %>%
  filter(velocity == "yes") %>%
  pull(Date)



ggplot(discharge_data, aes(x = Date, y = Discharge_cfs)) +
  geom_line() +
  geom_vline(xintercept = sample_dates, color = "#70A494", linetype = "solid") +
  geom_vline(xintercept = velocity_dates, color = "#DE8A5A", linetype = "solid") +
  labs(
       x = "Date",
       y = "Discharge (cfs)") +
  theme_minimal()+
  theme(
    panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank()   # remove minor grid lines
  )

# Here, you can see months better
ggplot(discharge_data, aes(x = Date, y = Discharge_cfs)) +
  geom_line() +
  geom_vline(xintercept = sample_dates, color = "#70A494") +
  geom_vline(xintercept = velocity_dates, color = "#DE8A5A", linetype = "solid") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(x = "Date", y = "Discharge (cfs)") +
  theme_minimal()

### DIDYMO OVER TIME W/ HYDROGRAPH -----------------------------------------------

# Pivot your algae columns to long format

didymo_over_time <- didymo_benthotorch %>%
  pivot_longer(
    cols = c(Cyano, Green, Diatoms),   # the algae columns
    names_to = "Algae_Type",
    values_to = "Concentration"
  ) %>%
  filter(Algae_Type == "Diatoms")

didymo_over_time <- didymo_over_time %>%
  mutate(
    Concentration = as.numeric(Concentration)
  )%>%
  filter(!is.na(Concentration)) # removing NA and making conc. numeric

didymo_over_time$Sampling_date <- as.character(didymo_over_time$Sampling_date)

str(didymo_over_time)

ggplot(didymo_over_time, aes(x = Sampling_date, y = Concentration, fill = Algae_Type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  
  labs(
    x = "Date",
    y = "Algae Concentration",
    fill = "Algae Type"
  ) +
  scale_fill_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"), 
    name = "Algae Type") + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.position = "top"
  ) # make sure date is a character to run this right



# Overlaying with hydrology
didymo_over_time <- didymo_over_time %>%
  mutate(Sampling_date = as.Date(Sampling_date))

scale_factor <- 0.01

ggplot() +
  ## Hydrograph (scaled)
  geom_line(
    data = discharge_data,
    aes(x = Date, y = Discharge_cfs * scale_factor),
    color = "grey40",
    linewidth = 0.8
  ) +
  
  ## Boxplots (grouped by sampling date)
  geom_boxplot(
    data = didymo_over_time,
    aes(x = Sampling_date,
        y = Concentration,
        group = Sampling_date),
    width = 3,
    fill = "#2887A1",
    color = "#008080",
    alpha = 0.7
  ) +
  
  ## Primary axis: concentration
  scale_y_continuous(
    name = expression("Diatom concentration (" * mu * "g " * "chl-a" / cm^2 * ")"),
    sec.axis = sec_axis(
      transform = ~ . / scale_factor,  # <-- must provide function
      name = "Discharge (CFS)"
    )
  ) +
  
  scale_x_date(name = "Date") +
  theme_bw()



# Overlaying with SI
Sampling_dates <- read_csv("Sampling_dates.csv")

Sampling_dates <- Sampling_dates %>%
  mutate(Sampling_date = mdy(Sampling_date))

ggplot() +
  # SI
  geom_vline(
    data = Sampling_dates,
    aes(xintercept = Sampling_date, color = Occasion), 
    linetype = "dashed",
    linewidth = 0.3
  ) +
  scale_color_manual(values = c( # For better visuals on poster
    MAY_2021 = "#B4C8A8",
    MAY_2022 = "#B4C8A8",
    MAY_2023 = "#B4C8A8",
    MAY_2024 = "#B4C8A8",
    MAY_2025 = "#B4C8A8",
    AUG_2021 = "#EDBB8A",
    AUG_2022 = "#EDBB8A",
    AUG_2023 = "#EDBB8A",
    AUG_2024 = "#EDBB8A",
    AUG_2025 = "#EDBB8A",
    OCT_2021 = "#DE8A5A",
    OCT_2022 = "#DE8A5A",
    OCT_2023 = "#DE8A5A",
    OCT_2024 = "#DE8A5A",
    OCT_2025 = "#DE8A5A")) +
  
  ## Hydrograph (scaled)
  geom_line(
    data = discharge_data,
    aes(x = Date, y = Discharge_cfs * scale_factor),
    color = "grey40",
    linewidth = 0.8
  ) +
  
  ## Boxplots (grouped by sampling date)
  geom_boxplot(
    data = didymo_over_time,
    aes(x = Sampling_date,
        y = Concentration,
        group = Sampling_date),
    width = 3,
    fill = "#2887A1",
    color = "#008080",
    alpha = 0.7
  ) +
  
  ## Axes
  scale_y_continuous(
    name = expression("Diatom concentration (" * mu * "g " * "chl-a" / cm^2 * ")"),
    sec.axis = sec_axis(
      transform = ~ . / scale_factor,  # <-- must provide function
      name = "Discharge (CFS)"
    )
  ) +
  
  scale_x_date(name = "Date") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.grid.minor = element_blank()
    #axis.line = element_line(color = "black")
  )

library(rcartocolor)
install.packages(rcartocolor)
mycolors <- carto_pal(7, "Earth")
mycolors


# DISCHARGE (30 DAY MEAN AVERAGE PRIOR TO SAMPLING) RELATIONSHIPS --------------------
# AVERAGED ACROSS REPLICATES

# Function to compute 30-day stats for each sample
get_30day_stats <- function(Sampling_date, discharge_data) {
  start_date <- Sampling_date - days(30)
  subset <- discharge_data %>%
    filter(Date > start_date & Date <= Sampling_date)
  
  mean_cfs <- mean(subset$Discharge_cfs, na.rm = TRUE)
  cv_cfs <- sd(subset$Discharge_cfs, na.rm = TRUE) / mean_cfs
  
  return(c(mean_30d = mean_cfs, cv_30d = cv_cfs))
}

# Apply function to each sample
stats_matrix <- t(sapply(didymo_benthotorch$Sampling_date, get_30day_stats, discharge_data = discharge_data))

# Combine with didymo_benthotorch
didymo_benthotorch <- cbind(didymo_benthotorch, stats_matrix)

# Check results
head(didymo_benthotorch)


# DISCHARGE RELATIONSHIPS ACROSS BENTHOTORCH REPS------------------------------
avg_didymo_benthotorch <- didymo_benthotorch %>%
  select(Sampling_date, Site, Sample.Type, 
         Cyano, Green, Diatoms, Chlorophyll.A,
         mean_30d, cv_30d) %>%
  group_by(Sampling_date, Site, Sample.Type) %>%
  summarise(
    across(c(Cyano, Green, Diatoms, Chlorophyll.A, mean_30d, cv_30d), 
           \(x) mean(x, na.rm = TRUE)),   # modern dplyr style
    Replicate_number = n(),   # how many replicates went into the average
    .groups = "drop"
  )

# Pivot your algae columns to long format
algae_long <- avg_didymo_benthotorch %>%
  pivot_longer(
    cols = c(Cyano, Green, Diatoms),   # the algae columns
    names_to = "Algae_Type",
    values_to = "Concentration"
  )

algae_long$Algae_Type <- factor(algae_long$Algae_Type, levels = c("Green", "Cyano", "Diatoms"))

# To look at boxplots before algae was averaged across replicates run this*****
# and switch out algae_long for nonavg_algae_long
nonavg_algae_long <- didymo_benthotorch %>%
  pivot_longer(
  cols = c(Cyano, Green, Diatoms),   # the algae columns
  names_to = "Algae_Type",
  values_to = "Concentration"
 )

nonavg_algae_long$Algae_Type <- factor(nonavg_algae_long$Algae_Type, levels = c("Green", "Cyano", "Diatoms"))

# Plotting algae over time and over discharge
library(rcartocolor)
display_carto_all()
carto_pal(7, "Geyser")

# Algae over time at each site
ggplot(nonavg_algae_long, aes(x = Sampling_date, y = Concentration, color = Algae_Type)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(aes(group = Algae_Type), method = "lm", se = FALSE, size = 1) +  # one line per algae type
  facet_wrap(~Site) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "Sampling Date",
    y = "Algae Concentration",
    color = "Algae Type",
  ) +
  scale_colour_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"), 
    name = "Algae Type") +  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# SAA with discharge fun
ggplot(nonavg_algae_long, aes(x = Sampling_date, y = Concentration, color = Algae_Type)) +
  geom_point(aes(size = mean_30d, alpha = 0.7)) +  # size represents mean discharge
  geom_smooth(aes(group = Algae_Type), method = "lm", se = FALSE, size = 1) +  
  facet_wrap(~Site) +
  scale_size_continuous(name = "Mean 30-day Discharge (cfs)") +
  labs(
    x = "Sampling Date",
    y = "Algae Concentration",
    color = "Algae Type",
  ) +
  scale_colour_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"), 
    name = "Algae Type") + 
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Algae response to discharge at each site
ggplot(nonavg_algae_long, aes(x = mean_30d, y = Concentration, color = Algae_Type)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(aes(group = Algae_Type), method = "lm", se = FALSE, size = 1) +  
  facet_wrap(~Site) +
  labs(
    x = "Mean 30-day Discharge (cfs)",
    y = "Algae Concentration",
    color = "Algae Type",
    shape = "Site"
  ) +
  scale_colour_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"), 
    name = "Algae Type") + 
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )


# Algae at each site with boxplots
ggplot(nonavg_algae_long, aes(x = Site, y = Concentration, fill = Algae_Type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # side-by-side boxes per site
  labs(
    x = "Site",
    y = "Algae Concentration",
    fill = "Algae Type"
  ) +
  scale_fill_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"), 
    name = "Algae Type") + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # tilt x labels if many sites
    legend.position = "top"
  )



# Boxplots across discharge...first have to make discharge categorical
bins <- nonavg_algae_long %>%
  mutate(discharge_bin = cut(mean_30d,
                             breaks = c(0, 50, 100, 200, 500, 1000, Inf),
                             labels = c("0–50", "50–100", "100–200", "200–500", "500–1000", "1000+")))

ggplot(bins, aes(x = discharge_bin, y = Concentration, fill = Algae_Type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # side-by-side boxes per site
  labs(
    x = "30-day mean discharge (binned)",
    y = "Algae Concentration",
    fill = "Algae Type"
  ) +
  scale_fill_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"), 
    name = "Algae Type") + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # tilt x labels if many sites
    legend.position = "top"
  )

# VELOCITY RELATIONSHIPS ACROSS BENTHOTORCH REPS--------------------------------------

#Filtering just for data with velocity
didymo_benthotorch_velocity <- didymo_benthotorch %>%
  filter(!is.na(Velocity)) 

sort(unique(didymo_benthotorch_velocity$Sampling_date))

# This is for when replicates are averaged
avg_didymo_benthotorch_velocity <- didymo_benthotorch_velocity %>%
  select(Sampling_date, Site, Sample.Type,  
         Cyano, Green, Diatoms, Chlorophyll.A, Velocity,
         mean_30d, cv_30d) %>%
  group_by(Sampling_date, Site, Sample.Type) %>%
  summarise(
    across(c(Cyano, Green, Diatoms, Chlorophyll.A, Velocity, mean_30d, cv_30d), 
           \(x) mean(x, na.rm = TRUE)),  
    Replicate_number = n(),   # how many replicates went into the average
    .groups = "drop"
  ) 

# When you want averaged replicates
velocity_long_avg <- avg_didymo_benthotorch_velocity %>%
  pivot_longer(
    cols = c(Cyano, Green, Diatoms),   # the algae columns
    names_to = "Algae_Type",
    values_to = "Concentration"
  )

velocity_long_avg$Algae_Type <- factor(velocity_long_avg$Algae_Type, levels = c("Green", "Cyano", "Diatoms"))

# To look at raw data before algae was averaged across replicates run this*****
# and switch out velocity_long for velocity_long_avg
velocity_long <- didymo_benthotorch_velocity %>%
  pivot_longer(
    cols = c(Cyano, Green, Diatoms),   # the algae columns
    names_to = "Algae_Type",
    values_to = "Concentration"
  ) 

velocity_long$Algae_Type <- factor(velocity_long$Algae_Type, levels = c("Green", "Cyano", "Diatoms"))
velocity_long <- velocity_long %>%
  mutate(Month_Year = format(Sampling_date, "%Y-%m"))


# Algae over time for each site
ggplot(velocity_long, aes(x = Sampling_date, y = Concentration, color = Algae_Type)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(aes(group = Algae_Type), method = "lm", se = FALSE, size = 1) +  # one line per algae type
  facet_wrap(~Site) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = "Sampling Date",
    y = "Algae Concentration",
    color = "Algae Type",
  ) +
  scale_colour_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"), 
    name = "Algae Type") +  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# SAA with velocity fun
ggplot(velocity_long, aes(x = Sampling_date, y = Concentration, color = Algae_Type)) +
  geom_point(aes(size = Velocity, alpha = 0.7)) +  # size represents mean Velocity
  geom_smooth(aes(group = Algae_Type), method = "lm", se = FALSE, size = 1) +  
  facet_wrap(~Site) +
  scale_size_continuous(name = "Velocity") +
  labs(
    x = "Sampling Date",
    y = "Algae Concentration",
    color = "Algae Type",
  ) +
  scale_colour_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"), 
    name = "Algae Type") + 
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Algae response to velocity at each site
ggplot(velocity_long, aes(x = Velocity, y = Concentration, color = Algae_Type)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(aes(group = Algae_Type), method = "lm", se = FALSE, size = 1) +
  facet_wrap(~Site, scales = "free") +
  labs(
    x = "Velocity",
    y = "Algae Concentration",
    color = "Algae Type",
    shape = "Site"
  ) +
  scale_colour_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"), 
    name = "Algae Type") + 
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )


# Algae at each site with boxplots
ggplot(velocity_long, aes(x = Site, y = Concentration, fill = Algae_Type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  # side-by-side boxes per site
  labs(
    x = "Site",
    y = "Algae Concentration",
    fill = "Algae Type"
  ) +
  scale_fill_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"), 
    name = "Algae Type") + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # tilt x labels if many sites
    legend.position = "top"
  )


# Categorizing by bins
vbins <- velocity_long %>% 
  mutate(velocity_bin = cut(Velocity,
                            breaks = c(0, 0.2,  0.4, 0.6, 1, Inf),
                            labels = c("0–0.2", "0.2–0.4", "0.4–0.6", 
                                       "0.6–0.1", "Above 1")))

# To look at diatom-velocity relationships
diatoms <- velocity_long %>%
  filter(Algae_Type == "Diatoms")

vbins <- diatoms %>%
  filter(Velocity >= 0) %>%   # removes negative velocities
  mutate(velocity_bin = cut(Velocity,
                            breaks = c(-Inf, 0.2, 0.3, 0.4, 0.6, 1, Inf),
                            labels = c("0–0.2", "0.2–0.3", 
                                       "0.3–0.4", "0.4–0.6", "0.6–1", "Above 1")))

# Algae across velocity "bins"
ggplot(vbins, aes(x = velocity_bin, y = Concentration, fill = Algae_Type)) + # Can run vbins with "diatoms" to only see diatoms
  geom_boxplot(position = position_dodge(width = 0.8)) +  # side-by-side boxes per site
  labs(
    x = "Velocity",
    y = "Algae Concentration",
    fill = "Algae Type"
  ) +
  scale_fill_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"), 
    name = "Algae Type") + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # tilt x labels if many sites
    legend.position = "top"
  )

# Linear diatom-velocity relationship
ggplot(diatoms, aes(x = Velocity, y = Concentration, color = Algae_Type)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, lwd = 1.2) +  # linear trend line
  labs(
    x = "Velocity",
    y = "Diatom Chl-a Concentration (µg/cm²)",
    color = "Algae Type"
  ) +
  scale_color_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1")
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )


# Algae response to velocity at each site at each sampling period
velocity_long %>%
  group_split(Site) %>%
  walk(~ {
    site_name <- unique(.x$Site)
    
    site_velocity <- ggplot(.x, aes(x = Velocity, y = Concentration, color = Algae_Type)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_smooth(aes(group = Algae_Type), method = "lm", se = FALSE, size = 1) +
      facet_wrap(~Month_Year, scales = "free") +   # only facet by time within site
      labs(
        x = "Velocity",
        y = "Algae Concentration",
        color = "Algae Type",
        title = paste("Site:", site_name)
      ) +
      scale_colour_manual(
        values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"),
        name = "Algae Type"
      ) +
      theme_bw(base_size = 14) +
      theme(
        legend.position = "top",
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()
      )
    
    print(site_velocity) 
    
  })


# Algae at each site with boxplots
ggplot(velocity_long_avg, aes(x = Site, y = Concentration, fill = Algae_Type)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +  
  labs(
    x = "Site",
    y = "Algae Concentration",
    fill = "Algae Type"
  ) +
  scale_fill_manual(
    values = c("Green" = "#70A494", "Cyano" = "#DE8A5A", "Diatoms" = "#2887A1"), 
    name = "Algae Type") + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.position = "top"
  )


# PLAYING WITH STATS-------------------------------------------------------------
# Not normal
shapiro.test(velocity_long$Concentration)


sites_x_algae <- velocity_long %>%
  group_by(Site) %>%
  rstatix::wilcox_test(Concentration ~ Algae_Type)
sites_x_algae # comparing algae type differences within sites

algae_x_sites <- velocity_long %>%
  group_by(Algae_Type) %>%
  rstatix::wilcox_test(Concentration ~ Site)
algae_x_sites # comparing algae type differences between sites



# Mixed models--------------

install.packages("lme4")      
install.packages("lmerTest")  
library(lme4)
library(lmerTest)

model <- lmer(Concentration ~ Velocity * Algae_Type + (1|Site) + (1|Month_Year),
              data = velocity_long)

summary(model)
# Significant relationship with diatoms and velocity in relation to the reference,
# which is green, so not really what we are interested in


# Just diatoms
# Velocity
diatoms <- velocity_long %>%
  filter(Algae_Type == "Diatoms")

model <- lmer(Concentration ~ Velocity  + (1|Site) + (1|Month_Year),
              data = diatoms) 
summary(model) # significant negative relationship: as velocity increases,
# diatom concentration decreases.

# without random effects
model <- lm(Concentration ~ Velocity,
              data = diatoms) 
summary(model)  


# Discharge
model <- lmer(Concentration ~ mean_30d  + (1|Site) + (1|Month_Year),
              data = diatoms) 
summary(model)
# No relationship with mean 30 d discharge before sampling
anova(model)


# COLORS
library(rcartocolor)
install.packages(rcartocolor)
mycolors <- carto_pal(7, "Earth")
mycolors
