# Sliding window analysis using climwin package
# 12/1/2025 Kanno Lab Meeting demonstration

# https://cran.r-project.org/web/packages/climwin/vignettes/climwin.html
install.packages("climwin")
library(climwin)

install.packages("quantreg")
library(quantreg)

# Set your working directory to your personal laptop
setwd("~/Library/CloudStorage/OneDrive-Colostate/Data/BVR")

# CSVs
# From USGS Green Mountain Gauge--mean discharge from 10/1/2022-9/29/2025
discharge_data <- read.csv("discharge_data.csv") # Qualifier A means it was approved, P = provisional

# Diatom conc. from benthotorch from 1/15/2023-8/27/2025
# ALL available diatom data so far
diatoms_discharge <- read.csv("diatoms_discharge.csv")

# Convert dates to Date class
discharge_data$Date <- as.Date(discharge_data$Date, format = "%m/%d/%Y")
diatoms_discharge$Sampling_date <- as.Date(diatoms_discharge$Sampling_date, format = "%m/%d/%Y")
discharge_data$Date <- as.Date(discharge_data$Date)
diatoms_discharge$Sampling_date <- as.Date(diatoms_discharge$Sampling_date)

# FIRST SLIDING WINDOW EXAMPLE

SlidWin <- slidingwin(xvar = list(climate = discharge_data$Discharge_cfs),
                      cdate = discharge_data$Date, #cdate = date variable for climate dataset
                      bdate = diatoms_discharge$Sampling_date, #bdate = date variable for biological dataset
                      #baseline = lm(Concentration ~ 1, data = diatoms_discharge), # no random effects
                      baseline = lme4::lmer(Concentration ~ 1 + (1|Site), data = diatoms_discharge),
                      cinterval = "day", # can also be weeks or months
                      range = c(30, 0), # can customize this to as many days out from biological sample day as you want
                      type = "relative", # bc we have no specific day in mind as a reference day
                      stat = c("max","mean"), # can choose if you want means or maxes or both
                      func = "lin") # can also do quadratic (“quad”), cubic (“cub”), 
                      # logarithmic (“log”) and inverse (“inv”)


head(SlidWin[[1]]$Dataset) # a sneak peak

SlidWin[[1]]$BestModel # BestModel: is a model object showing the relationship 
# between climate and biological date within the strongest climate window.

#When climate = 0, the predicted value of yvar is 0.5209
# For each 1-unit increase in climate (e.g., 1 cfs of discharge), 
# the model predicts an increase of 0.001069 units in diatom concentration.


# Now, plotting results!---------------------------------------
# a new object so we can plot examples from MeanWin
Output <- SlidWin[[1]]$Dataset

# distribution of ΔAICc values across all tested climate windows
# plotdelta shows you which windows are most predictive
# Red shows strong regions/windows, most predictive
plotdelta(dataset = Output)


# the higher a model weight value the greater confidence we have that this model 
# is the true ‘best’ model
# you can be we can be 95% confident that the best climate window falls within 
# 4% of the total fitted models
plotweights(dataset = Output)

# betalinear is interesting, it shows the relationships between the climate and biological data
# Blue shows positive relationships: if climate data goes up then biological data goes up
# Red shows negative relationship: if climate data goes up then biological data goes down or vice versa
plotbetas(dataset = Output)
# This is unexpected...tells us that the relationship between discharge and diatoms is low
# during the most predictive window
# We can interpret this as us targeting stable, low flows for sampling 
# This is where you can sleuth and see that there is also a decently predictive window
# that does have a strong negative relationship whcih is more what we'd expect


# shows boxplots of the start and end point of all climate windows that make up the 
# 95% confidence set. The values above each boxplot represent the median start and 
# end time for these models.
plotwin(dataset = Output)


# Accounting for overfitting-----------------
#By running slidingwin on this new randomized dataset we can determine the likelihood 
# of finding our original result by random chance
SlidWin_rand <- randwin(repeats = 5, 
                        xvar = list(climate = discharge_data$Discharge_cfs),
                        cdate = discharge_data$Date,
                        bdate = diatoms_discharge$Sampling_date,
                        #baseline = lm(Concentration ~ 1, data = diatoms_discharge),
                        baseline = lme4::lmer(Concentration ~ 1 + (1|Site), data = diatoms_discharge),
                        cinterval = "day",
                        range = c(30, 0),
                        type = "relative", 
                        stat = c("mean", "max"), 
                        func = "lin")

# Naming the model output as an object so we can plug it into plothist
datasetrand = SlidWin_rand[[1]]
SlidWin_rand[[1]]

pvalue(dataset = SlidWin[[1]]$Dataset, datasetrand = SlidWin_rand[[1]], metric = "C", sample.size = 3)
# Sample size is number of years in the dataset
# p < 0.05 original slidingwin result is unlikely to be an issue of overfitting


# This gives a histogram of the ΔAICc
# values extracted from randwin and a dashed line to show the ΔAICc
# of the observed result from slidingwin.
plothist(dataset = Output, datasetrand = datasetrand)




# we can plot the predictions of our best model over the biological data---------
SlidWin_single <- singlewin(xvar = list(climate = discharge_data$Discharge_cfs),
                            cdate = discharge_data$Date,
                            bdate = diatoms_discharge$Sampling_date,
                            #baseline = lm(Concentration ~ 1, data = diatoms_discharge),
                            baseline = lme4::lmer(Concentration ~ 1 + (1|Site), data = diatoms_discharge),
                            cinterval = "day",
                            range = c(3, 2), # change to single window range
                            type = "relative", 
                            stat = c("max"), 
                            func = "lin")

plotbest(dataset = Output,
         bestmodel = SlidWin_single$BestModel, 
         bestmodeldata = SlidWin_single$BestModelData)



# Plotting everything
plotall(dataset = Output,
        datasetrand = datasetrand,
        bestmodel = SlidWin_single$BestModel, 
        bestmodeldata = SlidWin_single$BestModelData)
