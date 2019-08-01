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
# setwd("C:/.../calc_caff/")

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
  cabone_mod <- mrgsolve::mread("model/cabone_caffpk.cpp")  # cabone model
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Number of individuals
  nsim <- 1
  trt <- with(dplyr::filter(pk_tb, !duplicated(ID)), table(TRT))
  uid <- sum(trt)  # Number of observed individuals
  nid <- uid*nsim  # Number of simulated individuals
  ID <- 1:nid  # Sequence of individual ID's
  
# Identify caffeine and plasma ID's
  id_caff <- dplyr::filter(pk_tb, TRT == "Caffeine" & !duplicated(ID)) %>% 
    dplyr::pull(ID)
  id_plac <- dplyr::filter(pk_tb, TRT == "Placebo" & !duplicated(ID)) %>% 
    dplyr::pull(ID)
  
# Source population and add relevant info to population data
  source("scripts/190730_Population_Study.R")
  pop_tb <- dplyr::mutate_at(pop_tb, paste0("ETA", 1:9), function(x) 0)
  pop_tb <- dplyr::filter(pk_tb, !duplicated(ID)) %>%
    dplyr::select(ID, P_0) %>%
    merge(pop_tb, by = "ID") %>%
    tibble::as_tibble()

# Create simulation input dataset
# Define time points
  conc_times <- seq(0, 48, by = 0.5)
  conc_tb <- tibble::tibble(
    ID = rep(ID, each = length(conc_times)),
    time = rep(conc_times, times = nid),
    amt = 0, evid = 0, rate = 0, cmt = 1
  )
  
# Create concentration dose dataset from observed data
  input_tb <- pk_tb %>%
    dplyr::filter(!is.na(time)) %>%
    dplyr::select(ID, time, amt, evid, rate, cmt) %>%
    dplyr::bind_rows(conc_tb) %>%
    dplyr::arrange(ID, time, desc(evid))
  
# Simulate standard dosing regimen
# Pipe dataset to model
  output_tb <- cabone_mod %>%
    mrgsolve::data_set(input_tb) %>%  # set input data (observations/events)
    mrgsolve::idata_set(pop_tb) %>%  # set individual data (ETAs/ICs)
    mrgsolve::carry_out(amt, evid, rate, cmt, SIM) %>%  # copy to simulated output
    mrgsolve::mrgsim() %>%  # simulate using mrgsolve
    tibble::as_tibble()
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Visual Plots
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
# Prepare data
  calamt_plas <- output_tb %>%
    dplyr::filter(cmt == 6) %>%
    dplyr::pull(P)
  plas_tb <- pk_tb %>%
    dplyr::filter(cmt == 6 & !is.na(time)) %>%
    dplyr::mutate(P = calamt_plas, RES = CALAMT - calamt_plas)
  calamt_urin <- output_tb %>%
    dplyr::filter(cmt == 31) %>%
    {dplyr::filter(., evid == 1)$UCA - dplyr::filter(., evid == 0)$UCA}
  urin_tb <- pk_tb %>%
    dplyr::filter(cmt == 31 & evid == 1) %>%
    dplyr::mutate(UCA = calamt_urin, RES = CALAMT - calamt_urin)
  
# Plasma Calcium Residuals vs. Time
  p <- NULL
  p <- ggplot(aes(x = time, y = RES), data = plas_tb)
  p <- p + ggtitle("Plasma Calcium Residuals vs. Time")
  p <- p + geom_point(shape = 1, colour = "blue")
  p <- p + geom_smooth(method = "loess", size = 1, colour = "red")
  p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  p <- p + scale_x_continuous("Time (hours)")
  p1 <- p + scale_y_continuous("Residual (mmol)", limits = c(-5, 5), 
    breaks = c(-5, -3, -1, 0, 1, 3, 5))
  ggsave("output/PlasResVsTime_iCAL.png", width = 10.7, height = 8.03)
  
  p1_facet <- p1 + facet_wrap(~TRT)
  ggsave("output/PlasResVsTimeFacet_iCAL.png", width = 10.7, height = 8.03)
  
# Plasma Calcium Residuals vs. Predicted
  p <- NULL
  p <- ggplot(aes(x = P, y = RES), data = plas_tb)
  p <- p + ggtitle("Plasma Calcium Residuals vs. Predicted Values")
  p <- p + geom_point(shape = 1, colour = "blue")
  p <- p + geom_smooth(method = "loess", size = 1, colour = "red")
  p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  p <- p + scale_x_continuous("Plasma Calcium (mmol)")
  p2 <- p + scale_y_continuous("Residual (mmol)", limits = c(-5, 5), 
    breaks = c(-5, -3, -1, 0, 1, 3, 5))
  ggsave("output/PlasResVsPred_iCAL.png", width = 10.7, height = 8.03)
  
  p2_facet <- p2 + facet_wrap(~TRT)
  ggsave("output/PlasResVsPredFacet_iCAL.png", width = 10.7, height = 8.03)
  
# Urine Calcium Residuals vs. Treatment Arm
  p <- NULL
  p <- ggplot()
  p <- p + ggtitle("Urine Calcium Residuals")
  p <- p + geom_boxplot(aes(x = "Caffeine", y = RES), 
    data = dplyr::filter(urin_tb, TRT == "Caffeine"))
  p <- p + geom_boxplot(aes(x = "Placebo", y = RES), 
    data = dplyr::filter(urin_tb, TRT == "Placebo"))
  p <- p + geom_hline(yintercept = 0, linetype = "dashed")
  p <- p + scale_x_discrete("Treatment Arm")
  p3 <- p + scale_y_continuous("Residual (mmol)", limits = c(-1, 22), 
    breaks = c(0, 5, 10, 15, 20))
  ggsave("output/UrinResVsTrt_iCAL.png", width = 10.7, height = 8.03)
  