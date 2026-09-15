Updated 8/25/26 with decriptors of R code and csvs


# DIATOM-VELOCITY

ALL_Bentho_Core.csv: benthotorch metasheet from 2023-present with paired velocity beginning in 2024 

Benthotorch-Velocity.R: Script with velocity and discharge relationships with the benthotorch data

Diatom_Velocity_RI.Rmd:ratio estimation project from Brien Gerber's class

diatom_velocity_RI.csv:the CSV I plugged into above R script, just a diatom-velocity datasheet

Sampling_dates: sampling dates associated with each site from 2021-2026, used in 
Benthotorch-Velocity.R to overlay hydrograph with SI and benthotorching


# ISOTOPIC NICHE

Meta_SIA_Data.csv:all up to date SIA data for fish and bugs

NICHE_WIDTHS_21to24.csv:niche widths for MTS and BNT across sampling seasons and years from 2021-2024

OVERLAP.csv: niche overlap for MTS and BNT across sampling seasons and years from 2021-2024

Sampling_dates: sampling dates associated with each site from 2021-2026, used in SI_ALL_21-24.R 
and SI_Fish_21-24.R for sampling occasions 

SI_ALL_21-24.R: looking at niche overlap for bugs and fish

SI_Fish_21-24.R: looking at niche overlap for BNT and MTS only


# CLIMWIN

ClimWin Demonstration.Rmd: Sliding window demonstration for lab meeting with diatom and discharge, updated 8/25/2026

ClimWin_demonstration.R: Same as above but in R script

SI_climwin.R: Sliding Window analysis for niche overlap and discharge--> revisit this sometime

SI_discharge.csv: discharge begining in Jan 2020-Dec 2024, to give plenty of window before the first SI sampling event

didymo_over_time.csv: all diatom data from 1/2023-Present, doesn't include velocity data, used in Cl