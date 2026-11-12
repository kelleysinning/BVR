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
install.packages("semPlot")   
library(semPlot)


bugs_and_fish <- read.csv("bugs_and_fish_with_diets_Unsummed.csv")

flow <- bugs_and_fish$Discharge_cfs
food.qual <- bugs_and_fish$Percent.Scraper.Biomass
diet.qual <- bugs_and_fish$Percent.Diet.Scraper.Biomass
fish.condish <- bugs_and_fish$Wr

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
              
              foodxfish := b1 * c1 # indirect effect of fish on food
              fishxflow := a1 * b1 * c1', # indirect effect of flow on fish
             data = data)

summary(path1)
# flow has a negative effect on food quality (e.g., percent scraper biomass in surbers)
# surber food quality has no effect on diet food quality, so basically fish are not necessarily
# selecting for food that is available in the macroinvertebrate community?
# diet quality has a weak negative effect on fish body condition, but this isn't significant...weird though
# Even if there are a lot of scrapers in a diet, if the assimilation of fatty acids is weakened then this 
# just may not be a strong predictor of body condition...body condition could be better explained from
# temp or something else

# no indirect effect of fish on food resources (e.g., there is a weak top-down effect)
# no indirect effect of flow on fish body condition
library(rcartocolor)
display_carto_all()
carto_pal(7, "Geyser")

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
         color = list(man = c("#008080", "#B4C8A8", "#EDBB8A", "#CA562C")),
         edge.color = "gray30", fade = FALSE)


# Now, let's explore blavaan!
# note you can specify the number of chains, burn-in, even priors, but the real
# value of this (from my perspective) is that it will spit out a working JAGS
# or Stan model if you specify mcmcfile = T
# https://rdrr.io/cran/blavaan/man/bsem.html

# Running this again with wet vs dry yearp-------------------------------------------

library(dplyr)
library(lubridate)


bugs_and_fish <- bugs_and_fish %>%
  mutate(
    year_type = case_when(
      year(Sampling_date.x) %in% c(2021, 2022, 2025) ~ "dry",
      year(Sampling_date.x) %in% c(2023, 2024)       ~ "wet",
      TRUE                                ~ NA_character_
    ),
    year_type = factor(year_type, levels = c("dry", "wet"))
  )

flow <- bugs_and_fish$Discharge_cfs
food.qual <- bugs_and_fish$Percent.Scraper.Biomass
diet.qual <- bugs_and_fish$Percent.Diet.Scraper.Biomass
fish.condish <- bugs_and_fish$Wr
year <- bugs_and_fish$year_type

library(lavaan)

data <- data.frame(flow, food.qual, diet.qual, fish.condish,
                   year = factor(year_type))   # e.g. "wet" / "dry"

# Check the group order lavaan will use (labels below must follow it)
levels(data$year)

model_mg <- '
  food.qual    ~ c(a1_dry, a1_wet) * flow
  diet.qual    ~ c(b1_dry, b1_wet) * food.qual
  fish.condish ~ c(c1_dry, c1_wet) * diet.qual

  # intercepts per group
  food.qual    ~ 1
  diet.qual    ~ 1
  fish.condish ~ 1

  # indirect effects per group
  foodxfish_dry := b1_dry * c1_dry
  foodxfish_wet := b1_wet * c1_wet
  flowxfish_dry := a1_dry * b1_dry * c1_dry
  flowxfish_wet := a1_wet * b1_wet * c1_wet

  # differences between groups
  diff_a    := a1_wet - a1_dry
  diff_b    := b1_wet - b1_dry
  diff_c    := c1_wet - c1_dry
  diff_flow := flowxfish_wet - flowxfish_dry
'

fit_mg <- sem(model_mg, data = data, group = "year")
summary(fit_mg, ci = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(runjags)

data <- data.frame(flow, food.qual, diet.qual, fish.condish)

path.blavaan <- bsem(model =
               'food.qual ~ 1 + a1 * flow  
              diet.qual ~ 1 + b1 * food.qual
              fish.condish ~ 1 + c1 * diet.qual
              
              foodxfish := b1 * c1 # indirect effect of fish on food
              fishxflow := a1 * b1 * c1', # indirect effect of flow on fish
              data = data,
              target = 'jags',
              mcmcfile = T)

summary(path.blavaan)

library(jagsUI)
sink("sem.jag")

# Plugging the blavaan code back in
cat("model { 
  for(i in 1:N) {
    food.qual[i] ~ dnorm(mu[i,1], 1/psi[1,1,g[i]])
    diet.qual[i] ~ dnorm(mu[i,2], 1/psi[2,2,g[i]])
    fish.condish[i] ~ dnorm(mu[i,3], 1/psi[3,3,g[i]])

  }

  # mu definitions
  for(i in 1:N) {
    mu[i,1] <- alpha[1,1,g[i]] + beta[1,4,g[i]]*flow[i]
    mu[i,2] <- alpha[2,1,g[i]] + beta[2,1,g[i]]*food.qual[i]
    mu[i,3] <- alpha[3,1,g[i]] + beta[3,2,g[i]]*diet.qual[i]

  }

  # Assignments from parameter vector & equality constraints
  alpha[1,1,1] <- parvec[1]
  beta[1,4,1] <- parvec[2]
  alpha[2,1,1] <- parvec[3]
  beta[2,1,1] <- parvec[4]
  alpha[3,1,1] <- parvec[5]
  beta[3,2,1] <- parvec[6]
  psi[1,1,1] <- pow(parvec[7],-1)
  psi[2,2,1] <- pow(parvec[8],-1)
  psi[3,3,1] <- pow(parvec[9],-1)
  psi[4,4,1] <- 61713.0979378224
  alpha[4,1,1] <- 332.78362745098
  def[1,1,1] <- beta[2,1,1] * beta[3,2,1]
  def[2,1,1] <- beta[1,4,1] * beta[2,1,1] * beta[3,2,1]
  # Priors
  parvec[1] ~ dnorm(0,1e-3)
  parvec[2] ~ dnorm(0,1e-2)
  parvec[3] ~ dnorm(0,1e-3)
  parvec[4] ~ dnorm(0,1e-2)
  parvec[5] ~ dnorm(0,1e-3)
  parvec[6] ~ dnorm(0,1e-2)
  parvec[7] ~ dgamma(1,.5)
  parvec[8] ~ dgamma(1,.5)
  parvec[9] ~ dgamma(1,.5)
}
")
sink()


library(jagsUI)


jags.data <- list(
  N            = nrow(data),
  g            = rep(1, nrow(data)),
  flow         = as.numeric(scale(data$flow)),
  food.qual    = data$food.qual,
  diet.qual    = data$diet.qual,
  fish.condish = data$fish.condish
)

params <- c("food.qual", "diet.qual", "fishcondish")

inits <- function() list(
  parvec = c(rnorm(6), rgamma(3, 1, 0.5))
)

fit <- jags(data = jags.data,
            inits = inits,
            parameters.to.save = params,
            model.file = "sem.jag",
            n.chains = 3,
            n.adapt = 1000,
            n.iter = 25000,
            n.burnin = 5000,
            n.thin = 5,
            parallel = TRUE)

print(fit)
traceplot(fit)



