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
  filter(!Species %in% c("Fry", "Fish Eggs")) %>%   # remove fish eggs them here
  mutate(
    Group = case_when(
      Species %in% c("BNT", "MTS") ~ Species,
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


# Example usage
occasions <- c("MAY_2021", "AUG_2021", "OCT_2021",
               "MAY_2022", "AUG_2022", "OCT_2022",
               "MAY_2023", "AUG_2023", "OCT_2023",
               "MAY_2024", "AUG_2024", "OCT_2024")


all_communities <- rownames(SIA_SIBER_OBJECT$sample.sizes)  # MAY_2021, AUG_2021, etc.
all_groups <- colnames(SIA_SIBER_OBJECT$sample.sizes)        # BNT, Macro, MTS


overlap_results <- list()

for(comm in all_communities) {
  
  # Only include groups that have data in this community (no NA in sample.sizes)
  valid_groups <- all_groups[!is.na(SIA_SIBER_OBJECT$sample.sizes[comm, ])]
  
  if(length(valid_groups) < 2) next  # skip if fewer than 2 groups
  
  community_groups <- paste(comm, valid_groups, sep = ".")
  pairs <- combn(community_groups, 2, simplify = FALSE)
  
  comm_results <- lapply(pairs, function(p) {
    maxLikOverlap(p[1], p[2], SIA_SIBER_OBJECT, p.interval = 0.40, n = 100)
  })
  
  names(comm_results) <- sapply(pairs, paste, collapse = " vs ")
  overlap_results[[comm]] <- comm_results
}

results_df <- do.call(rbind, lapply(names(overlap_results), function(comm) {
  do.call(rbind, lapply(names(overlap_results[[comm]]), function(pair) {
    res <- overlap_results[[comm]][[pair]]
    data.frame(
      community = comm,
      pair      = pair,
      overlap   = res["overlap"],
      area1     = res["area.1"],
      area2     = res["area.2"],
      row.names = NULL
    )
  }))
}))



results_df <- do.call(rbind, lapply(names(overlap_results), function(comm) {
  do.call(rbind, lapply(names(overlap_results[[comm]]), function(pair) {
    res <- overlap_results[[comm]][[pair]]
    
    area1   <- res[1]
    area2   <- res[2]
    overlap <- res[3]
    
    prop_overlap <- overlap / (area1 + area2 - overlap)  # matches your prior formula
    
    data.frame(
      community    = comm,
      pair         = pair,
      area1        = area1,
      area2        = area2,
      overlap_area = overlap,
      prop_overlap = prop_overlap,
      row.names    = NULL
    )
  }))
}))

print(results_df)



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


