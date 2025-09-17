# Code for Kelley Sinning Ph.D.
install.packages("dataRetrieval", type = "source")
library(dataRetrieval)
packageVersion("dataRetrieval")  # should be >= 2.7.19


# https://doi-usgs.github.io/dataRetrieval/reference/read_waterdata_daily.html
# read_waterdata_daily is newer but not working here for some reason.....
site <- "09057500"           # USGS site number (Blue River below Green Mountain)
parameter_code <- "00060"    # Parameter code for discharge (ft³/s)
statistic_id <- "00003"      # Statistic code for daily mean
start_date <- "2022-12-01"
end_date <- "2025-01-29"

# Retrieve daily discharge data
discharge_data <- read_waterdata_daily(
  monitoring_location_id = site,
  parameter_code = parameter_code,
  statistic_id = statistic_id,
  time = c(start_date, end_date)
)


# So using readNWISdv
site <- "09057500"          # USGS site number
parameterCd <- "00060"      # Discharge (cfs)
statCd <- "00003"           # Daily mean
startDate <- "2022-12-01"
endDate <- "2025-08-31"

# Retrieve daily values
discharge_data <- readNWISdv(
  siteNumbers = site,
  parameterCd = parameterCd,
  statCd = statCd,
  startDate = startDate,
  endDate = endDate
)

# Renaming columns

library(dplyr)

discharge_data <- discharge_data %>%
  rename(Discharge_cfs = X_00060_00003,
         Qualifier = X_00060_00003_cd
  )
 # A, P = data qualifiers (approved, provisional) for Qualifier

# Now, bringing in my algae data -----------------------------------------------
setwd("/Users/kelleysinning/Library/CloudStorage/OneDrive-Colostate/Data/BVR")
didymo_benthotorch <- read.csv("ALL_bentho_core.csv")
didymo_benthotorch <- didymo_benthotorch %>%
  filter(!is.na(Sampling_date)) # removing NA columns that arose from comments in the csv

# Making sure date formats are the same
didymo_benthotorch$Sampling_date <- as.Date(didymo_benthotorch$Sampling_date)  
discharge_data$Date <- as.Date(discharge_data$Date)

# Making sure columns are numeric
str(didymo_benthotorch)
didymo_benthotorch <- didymo_benthotorch %>%
  mutate(across(c(Cyano, Green, Diatoms, Chlorophyll.A, Velocity),
                as.numeric))


library(lubridate)

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


# Averaging across replicates

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



library(ggplot2)
library(dplyr)
library(tidyr)

# Pivot your algae columns to long format

algae_long <- avg_didymo_benthotorch %>%
  pivot_longer(
    cols = c(Cyano, Green, Diatoms),   # the algae columns
    names_to = "Algae_Type",
    values_to = "Concentration"
  )

# Plotting
library(rcartocolor)
display_carto_all()
carto_pal(7, "Temps")

# First, ordering things how I like
algae_long$Algae_Type <- factor(algae_long$Algae_Type, levels = c("Green", "Cyano", "Diatoms"))

# Algae over time at each site
ggplot(algae_long, aes(x = Sampling_date, y = Concentration, color = Algae_Type)) +
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
ggplot(algae_long, aes(x = Sampling_date, y = Concentration, color = Algae_Type)) +
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
ggplot(algae_long, aes(x = mean_30d, y = Concentration, color = Algae_Type)) +
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
    panel.grid.minor = element_bla