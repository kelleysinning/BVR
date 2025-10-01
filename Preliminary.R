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


# So using readNWISdv...it's older but works
site <- "09057500"          # USGS site number
parameterCd <- "00060"      # Discharge (cfs)
statCd <- "00003"           # Daily mean
startDate <- "2022-10-01"
endDate <- "2025-09-29"

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
    panel.grid.minor = element_blank()
  )


# Algae at each site with boxplots
ggplot(algae_long, aes(x = Site, y = Concentration, fill = Algae_Type)) +
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

# To look at boxplots before algae was averaged across replicates run this # code below
# and switch out algae_long for nonavg_algae_long

  # nonavg_algae_long <- didymo_benthotorch %>%
  # pivot_longer(
  #  cols = c(Cyano, Green, Diatoms),   # the algae columns
  #  names_to = "Algae_Type",
  # values_to = "Concentration"
  # )

bins <- algae_long %>%
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

# VELOCITY RELATIONSHIPS WITH AVERAGED REPS--------------------------------------

#Filtering just for data with velocity
didymo_benthotorch_velocity <- didymo_benthotorch %>%
  filter(!is.na(Velocity)) 

# This is for when replicates were avergaed
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

velocity_long_avg <- avg_didymo_benthotorch_velocity %>%
  pivot_longer(
    cols = c(Cyano, Green, Diatoms),   # the algae columns
    names_to = "Algae_Type",
    values_to = "Concentration"
  )

# Ordering things how I like
velocity_long_avg$Algae_Type <- factor(velocity_long_avg$Algae_Type, levels = c("Green", "Cyano", "Diatoms"))


# Plotting

# Algae over time for each site
ggplot(velocity_long_avg, aes(x = Sampling_date, y = Concentration, color = Algae_Type)) +
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
ggplot(velocity_long_avg, aes(x = Sampling_date, y = Concentration, color = Algae_Type)) +
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
ggplot(velocity_long_avg, aes(x = Velocity, y = Concentration, color = Algae_Type)) +
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
ggplot(velocity_long_avg, aes(x = Site, y = Concentration, fill = Algae_Type)) +
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


# Boxplots across velocity...first have to make velocity categorical-------

# This is still for averaged replicates, use velocity_long below to see non avg
# playing around with different blocks

vbins <- velocity_long_avg %>%
  mutate(velocity_bin = cut(Velocity,
                            breaks = c(0, 0.2, 0.3, 0.4, 0.6, 1, Inf),
                            labels = c("0–0.2", "0.2–0.3", 
                                       "0.3–0.4", "0.4–0.6", "0.6–1", "Above 1")))

vbins <- velocity_long_avg %>%
  mutate(velocity_bin = cut(Velocity,
                            breaks = c(0, 0.2,  0.4, 0.6, 1, Inf),
                            labels = c("0–0.2", "0.2–0.4", "0.4–0.6", 
                                       "0.6–0.1", "Above 1")))

# this is for unaveraged, raw data for diatoms only
vbins <- diatoms %>%
  filter(Velocity >= 0) %>%   # removes negative velocities
  mutate(velocity_bin = cut(Velocity,
                            breaks = c(-Inf, 0.2, 0.3, 0.4, 0.6, 1, Inf),
                            labels = c("0–0.2", "0.2–0.3", 
                                       "0.3–0.4", "0.4–0.6", "0.6–1", "Above 1")))

ggplot(vbins, aes(x = velocity_bin, y = Concentration, fill = Algae_Type)) +
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


install.packages("rstatix")
library(rstatix)

# Not normal
shapiro.test(velocity_long$Concentration)


# Non-parametric  pairwise comparisons (within bins)
wilcox_test <- vbins %>%
  group_by(velocity_bin) %>%
  wilcox_test(Concentration ~ Algae_Type)





# NOT AVERAGING REPLICATES------------------------------------------------------
# trying to do one box per site and per sampling

velocity_long <- didymo_benthotorch_velocity %>%
  pivot_longer(
    cols = c(Cyano, Green, Diatoms),   # the algae columns
    names_to = "Algae_Type",
    values_to = "Concentration"
  ) # This is what we want to use


velocity_long <- velocity_long %>%
  mutate(Month_Year = format(Sampling_date, "%Y-%m"))


# Ordering things how I like
velocity_long$Algae_Type <- factor(velocity_long$Algae_Type, levels = c("Green", "Cyano", "Diatoms"))


# Algae response to velocity at each site
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
# This goes with algae x sites below
ggplot(velocity_long, aes(x = Site, y = Concentration, fill = Algae_Type)) +
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


# Discharge
model <- lmer(Concentration ~ mean_30d  + (1|Site) + (1|Month_Year),
              data = diatoms) 
summary(model)
# No relationship with mean 30 d discharge before sampling
anova(model)




# MOVING WINDOW ----------------------------------------------------------------
# https://cran.r-project.org/web/packages/climwin/vignettes/climwin.html
install.packages("climwin")
library(climwin)

# Filtering all available diatom data, not just the ones with velocity
diatoms_discharge <- didymo_benthotorch %>%
  pivot_longer(
    cols = c(Cyano, Green, Diatoms),   # the algae columns
    names_to = "Algae_Type",
    values_to = "Concentration"
  ) # This is what we want to use

# Now, filtering just for algae
diatoms_discharge <- diatoms_discharge %>%
  filter(Algae_Type == "Diatoms")

# Remove rows with NA concentrations or missing dates
diatoms_discharge <- diatoms_discharge %>%
  filter(!is.na(Concentration), !is.na(Sampling_date))


# Convert dates to Date class
discharge_data$Date <- as.Date(discharge_data$Date)
diatoms_discharge$Sampling_date <- as.Date(diatoms_discharge$Sampling_date)

discharge_data <- discharge_data %>%
  filter(!is.na(Discharge_cfs), !is.na(Date)) # removing anything weird


# Add year and day-of-year columns to biological dataset
diatoms_discharge$refday <- as.numeric(format(diatoms_discharge$Sampling_date, "%j"))
diatoms_discharge$year   <- as.numeric(format(diatoms_discharge$Sampling_date, "%Y"))



MeanWin <- slidingwin(xvar = list(climate = discharge_data$Discharge_cfs),
                      cdate = discharge_data$Date,
                      bdate = diatoms_discharge$Sampling_date,
                      baseline = lm(Concentration ~ 1, data = diatoms_discharge),
                      cinterval = "day",
                      range = c(30, 0),
                      #type = "absolute", refday = c(20, 9), #not exactly sure what ref date we should use
                      type = "relative",
                      stat = "mean", #can also do max
                      func = "lin")


head(MeanWin[[1]]$Dataset)

MeanWin[[1]]$BestModel
MeanOutput <- MeanWin[[1]]$Dataset

plotdelta(dataset = MeanOutput) # Red shows strong regions/windows
plotweights(dataset = MeanOutput)
plotbetas(dataset = MeanOutput)
plotwin(dataset = MeanOutput)



MeanWin_single <- singlewin(xvar = list(climate = discharge_data$Discharge_cfs),
                        cdate = discharge_data$Date,
                        bdate = diatoms_discharge$Sampling_date,
                        baseline = lm(Concentration ~ 1, data = diatoms_discharge),
                        cinterval = "day",
                        range = c(30, 0),
                        type = "absolute", refday = c(20, 8), #not exactly sure what ref date we should use
                        stat = "mean", #can also do max
                        func = "lin")


plotbest(dataset = MeanOutput,
         bestmodel = MeanWin_single$BestModel, 
         bestmodeldata = MeanWin_single$BestModelData)

# Accounting for overfitting-----------------
MeanWin_rand <- randwin(repeats = 5, 
                        xvar = list(climate = discharge_data$Discharge_cfs),
                        cdate = discharge_data$Date,
                        bdate = diatoms_discharge$Sampling_date,
                        baseline = lm(Concentration ~ 1, data = diatoms_discharge),
                        cinterval = "day",
                        range = c(30, 0),
                        type = "absolute", refday = c(20, 05), #not exactly sure what ref date we should use
                        stat = "mean", #can also do max
                        func = "lin")

datasetrand = MeanWin_rand[[1]]
pvalue(MeanOutput <- MeanWin[[1]]$Dataset, datasetrand = MeanWin_rand[[1]], metric = "C", sample.size = 47)
# p < 0.05 original slidingwin result is unlikely to be an issue of overfitting

plothist(dataset = MassOutput, datasetrand = MassRand)
