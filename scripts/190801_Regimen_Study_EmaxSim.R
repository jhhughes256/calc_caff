# Caffine PopPK Model Simulation - Regimen
# ------------------------------------------------------------------------------
# Simulate dosing regimen for caffeine plasma concentration simulation.
# Try to replicate the simulation component of the VPC in Perera et al.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Prepare workspace
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
  setwd("E:/Hughes/Git/calc_caff")

# Load package libraries
  library(dplyr)	# dplyr required for mrgsolve
  library(mrgsolve)	 # Metrum differential equation solver for pharmacometrics
  library(ggplot2)  # graphical visualisation
  # library(readxl)  # reading in observed data
  # library(tidyr)  # for nesting and unnesting data
  # library(purrr)  # iterative functions
  # library(MASS)  # mvrnorm in `190730_Population_Study.R`

# Source external scripts
  source("scripts/functions_utility.R")  # functions utility
  source("model/caffpk_perera_MA.R")  # PopPK model script
  source("scripts/190730_Data_Preparation.R")  # Observed calcium data

# Read in cabone model
  cabone_mod <- mrgsolve::mread("model/cabone_caffpkpd.cpp")  # cabone model

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Number of individuals
  nsim <- 1000
  trt <- with(dplyr::filter(pk_tb, !duplicated(ID)), table(TRT))
  uid <- sum(trt)  # Number of observed individuals
  nid <- uid*nsim  # Number of simulated individuals
  ID <- 1:nid  # Sequence of individual ID's

# Source population and add relevant info to population data
  source("scripts/190730_Population_Study.R")
  pop_tb <- dplyr::mutate(pop_tb, ID_ORIG = rep(unique(pk_tb$ID_ORIG), nsim))
  pop_tb <- dplyr::filter(pk_tb, !duplicated(ID)) %>%
    dplyr::select(ID_ORIG, P_0, GFR_0) %>%
    dplyr::inner_join(pop_tb, by = "ID_ORIG") %>%
    dplyr::arrange(ID) %>%
    dplyr::mutate(PDEMAX = 0.446)  # standard error (hessian): sigma: 0.760

# Create simulation input dataset
# Define time points
  conc_times <- seq(0, 48, by = 0.5)
  conc_tb <- tibble::tibble(
    ID = rep(ID, each = length(conc_times)),
    time = rep(conc_times, times = nid),
    amt = 0, evid = 0, rate = 0, cmt = 1
  )

# Create input and run model all in one pipe (to save on memory usage)
# Begin benchmark
  tictoc::tic()
# Create concentration dose dataset from observed data
  purrr::map_dfr(0:(nsim - 1)*uid, function(x) {
      dplyr::mutate(pk_tb, ID = x + ID)
    }) %>%
    dplyr::filter(!is.na(time)) %>%
    dplyr::select(ID, time, amt, evid, rate, cmt, SIM, ID_ORIG, TRT) %>%
    dplyr::bind_rows(conc_tb) %>%
    dplyr::arrange(ID, time, desc(evid)) %>%
# Pipe dataset into cabone model
    {mrgsolve::data_set(cabone_mod, .)} %>%  # set input data
    mrgsolve::idata_set(pop_tb) %>%  # set individual data (ETAs/ICs)
    mrgsolve::carry_out(amt, evid, rate, cmt, SIM, ID_ORIG, TRT) %>%  # output
    mrgsolve::mrgsim() %>%  # simulate using mrgsolve
    tibble::as_tibble() %>%
    readr::write_rds("output/EmaxSim.rds")
# Finish benchmark (simulation expected to take ~30 hours)
  tictoc::toc()
