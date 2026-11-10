# SEM code for Fall 2026

install.packages("jagsUI")
library(jagsUI)
install.packages("rstan")
library(rstan)
install.packages("rtools")
library(rtools)
install.packages("vioplot")
library(vioplot)
install.packages("lavaan")
library(lavaan)
install.packages("piecewiseSEM")
library(piecewiseSEM)
install.packages("blavaan")
library(blavaan)


data <- read.csv("bugs_and_fish_with_diets_Unsummed.csv")

flow <- data$Discharge_cfs
food.qual <- data$Percent.Scraper.Biomass
diet.qual <- data$Percent.Diet.Scraper.Biomass
fish.condish <- data$Wr

# To write a path analysis we'll use 'lavaan', a package that allows for the construction of
# structural equation models. Lavaan can obtain ML parameter estimates for path coefficients
# produce cool figures, etc., but we're going to stick in BUGS-code land because lavaan
# CAN'T handle: 1) link functions, 2) joint analysis with other (e.g., capture-mark-recapture)
# data types, or 3) missing data.
data <- data.frame(flow, food.qual, diet.qual, fish.condish)
path1 <- sem(model =
              'food.qual ~ 1 + a1 * flow  
              diet.qual ~ 1 + b1 * food.qual
              fish.condish ~ 1 + c1 * diet.qual
              
              foodxfish := b1 * c1
              fishxflow := a1 * b1 * c1',   # this allows us to estimate the indirect effect of wolves on veg
             # also note we can specify parameter names!
             data = data)

summary(path1)

print(summary(path1))

install.packages("semPlot")   # once
library(semPlot)

semPaths(path1,
         what = "est",            # show estimated coefficients on arrows
         whatLabels = "est",
         layout = "tree2",        # left-to-right chain; try "tree", "spring", "circle"
         rotation = 2,            # 2 = left to right (default 1 is top to bottom)
         edge.label.cex = 1.1,
         sizeMan = 10,            # box size
         style = "lisrel",
         intercepts = FALSE,      # hides the intercept triangles, which clutter the plot
         residuals = FALSE,       # hides the residual loops
         nCharNodes = 0,          # don't abbreviate variable names
         color = list(man = c("lightblue", "lightgreen", "lightgreen", "plum")),
         edge.color = "gray30")

# Now, let's explore blavaan!
# note you can specify the number of chains, burn-in, even priors, but the real
# value of this (from my perspective) is that it will spit out a working JAGS
# or Stan model if you specify mcmcfile = T
# https://rdrr.io/cran/blavaan/man/bsem.html
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(runjags)
m.blavaan <- bsem(model = '
                w ~ 1 + h # here, the indirect effect is captured
                u ~ 1 + w + h', data = data, 
                  target = 'jags',
                  mcmcfile = T)
summary(m.blavaan)







