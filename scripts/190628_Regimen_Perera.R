# Caffine PopPK Model Simulation - Regimen
# ------------------------------------------------------------------------------
# Simulate dosing regimen for caffeine plasma concentration simulation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare workspace
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../calc_caff/")

# Load package libraries
  library(dplyr)	# dplyr required for mrgsolve
  library(mrgsolve)	 # Metrum differential equation solver for pharmacometrics
  library(ggplot2)  # graphical visualisation
  # library(MASS)  # mvrnorm

# Source external scripts
  # source("scripts/functions_utility.R")  # functions utility
  # source("model/caffpk_perera_MA.R")  # PopPK model script 
  source("model/caffpk_perera_noMA.R")  # PopPK model script 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Number of individuals
  nid <- 12  # Number of individuals
  ID <- 1:nid  # Sequence of individual ID's
  
# Source population
  source("scripts/190628_Population_MA.R")
  
# Create simulation input dataset
# Define time points
  conc_times <- seq(from = 0, 24, by = 0.5)  # 1 day of half hourly data
  dose_times <- 0  # 1 dose
  dose_amt <- 200
  
# Duplicate population data for number of times and alter dose times to include
#   dosage events for mrgsolve
  input_tb <- tibble::tibble(
    ID = rep(ID, each = length(conc_times)),
    time = rep(conc_times, nid),
    amt = dplyr::if_else(time %in% dose_times, dose_amt, 0),
    cmt = 1,
    evid = dplyr::if_else(amt != 0, 1, 0),
    rate = 0
  )
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulate standard dosing regimen
# Pipe dataset to model
  output_tb <- mod %>%
    mrgsolve::data_set(input_tb) %>%  # set input data (observations/events)
    mrgsolve::idata_set(pop_tb) %>%  # set individual data (sets tumour size)
    mrgsolve::carry_out(amt, evid, rate, cmt) %>%  # copy to simulated output
    mrgsolve::mrgsim() %>%  # simulate using mrgsolve
    tibble::as_tibble()
  
# Save data as .RDS for later use
  # saveRDS(output_tb, "output/pkdata_perera.rds")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot patient data
# Set ggplot2 theme
  # theme_bw2 <- theme_set(theme_bw(base_size = 14))
  # theme_update(plot.title = element_text(hjust = 0.5))

# Plot individual patient caffeine concentrations
  # p <- NULL
  # p <- ggplot(data = output_flat_tb)
  # p <- p + geom_line(aes(x = time, y = IPRECA), colour = "red", size = 1)
  # p <- p + geom_point(aes(x = time, y = DVCA), colour = "blue", shape = 1, alpha = 0.5)
  # p <- p + labs(x = "Time (days)", y = "Caffeine Concentration (mg/L)")
  # p <- p + coord_cartesian(xlim = c(0, 24), ylim = NULL)
  # p <- p + facet_wrap(~ID)
  # p

# Plot individual patient paraxanthine concentrations
  # p <- NULL
  # p <- ggplot(data = output_flat_tb)
  # p <- p + geom_line(aes(x = time, y = IPREPX), colour = "red", size = 1)
  # p <- p + geom_point(aes(x = time, y = DVPX), colour = "blue", shape = 1, alpha = 0.5)
  # p <- p + labs(x = "Time (days)", y = "Paraxanthine Concentration (mg/L)")
  # p <- p + coord_cartesian(xlim = c(0, 24), ylim = NULL)
  # p <- p + facet_wrap(~ID)
  # p
  
# Plot both individual patient caffeine and paraxanthine concentrations
  # p <- NULL
  # p <- ggplot(data = output_flat_tb)
  # p <- p + geom_line(aes(x = time, y = IPRECA), colour = "red", size = 1)
  # p <- p + geom_line(aes(x = time, y = IPREPX), colour = "blue", size = 1)
  # p <- p + labs(x = "Time (days)", y = "Caffeine Concentration (mg/L)")
  # p <- p + coord_cartesian(xlim = c(0, 24), ylim = NULL)
  # p <- p + facet_wrap(~ID)
  # p  
  