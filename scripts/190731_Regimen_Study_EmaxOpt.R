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
  cabone_mod <- mrgsolve::mread("model/cabone_caffpkpd.cpp")  # cabone model
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Number of individuals
  nsim <- 1
  trt <- with(dplyr::filter(pk_tb, !duplicated(ID)), table(TRT))
  uid <- sum(trt)  # Number of observed individuals
  nid <- uid*nsim  # Number of simulated individuals
  ID <- 1:nid  # Sequence of individual ID's
  
# Source population and add relevant info to population data
  source("scripts/190730_Population_Study.R")
  pop_tb <- dplyr::mutate_at(pop_tb, paste0("ETA", 1:9), function(x) 0)
  pop_tb <- dplyr::filter(pk_tb, !duplicated(ID)) %>%
    dplyr::select(ID, P_0, GFR_0) %>%
    merge(pop_tb, by = "ID") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(PDEMAX = 0.446)  # 0.35
  
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
  
# Optimise Emax
  # opt_emax_fn <- function(par) {
  # # Fill in EMAX
  #   popin_tb <- dplyr::mutate(pop_tb, PDEMAX = par[1])
  # # Pipe dataset to model
  #   output_tb <- cabone_mod %>%
  #     mrgsolve::data_set(input_tb) %>%  # set input data (observations/events)
  #     mrgsolve::idata_set(popin_tb) %>%  # set individual data (ETAs/ICs)
  #     mrgsolve::carry_out(amt, evid, rate, cmt, SIM) %>%  # copy to simulated output
  #     mrgsolve::mrgsim() %>%  # simulate using mrgsolve
  #     tibble::as_tibble()
  # # Extract urine data
  #   calamt_urin <- output_tb %>%
  #     dplyr::filter(cmt == 31) %>%
  #     {dplyr::filter(., evid == 1)$UCA - dplyr::filter(., evid == 0)$UCA}
  # # Calculate log-likelihood
  #   loglik <- dnorm(log(urine_tb$CALAMT), log(calamt_urin), par[2], log = T)
  #   return(-2*sum(loglik, na.rm = T))
  # }
  # init_par <- c(0.35, 1)
  # est_emax <- optim(init_par, opt_emax_fn, method = "L-BFGS-B",
  #   lower = c(0.01, 0.01), upper = c(0.99, 10), hessian = TRUE,
  #   control = list(
  #   parscale = init_par, fnscale = opt_emax_fn(init_par),
  #   factr = 1e12, pgtol = 1e-8
  # ))
  # theta = 0.446; sigma = 0.760
  
# Simulate standard dosing regimen
# Pipe dataset to model
  tictoc::tic()
  output_tb <- cabone_mod %>%
    mrgsolve::data_set(input_tb) %>%  # set input data (observations/events)
    mrgsolve::idata_set(pop_tb) %>%  # set individual data (ETAs/ICs)
    mrgsolve::carry_out(amt, evid, rate, cmt, SIM) %>%  # copy to simulated output
    mrgsolve::mrgsim() %>%  # simulate using mrgsolve
    tibble::as_tibble()
  tictoc::toc()
  
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
    dplyr::mutate(UCA = calamt_urin, RES = CALAMT - calamt_urin, logRES = log(CALAMT) - log(calamt_urin))
  
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
  ggsave("output/PlasResVsTime_EmaxOpt.png", width = 10.7, height = 8.03)
  
  p1_facet <- p1 + facet_wrap(~TRT)
  ggsave("output/PlasResVsTimeFacet_EmaxOpt.png", width = 10.7, height = 8.03)
  
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
  ggsave("output/PlasResVsPred_EmaxOpt.png", width = 10.7, height = 8.03)
  
  p2_facet <- p2 + facet_wrap(~TRT)
  ggsave("output/PlasResVsPredFacet_EmaxOpt.png", width = 10.7, height = 8.03)
  
# Urine Calcium Residuals vs. Treatment Arm
  p <- NULL
  p <- ggplot(data = urin_tb)
  p <- p + ggtitle("Urine Calcium Residuals")
  p <- p + geom_boxplot(aes(x = "Predicted", y = UCA))
  p <- p + geom_boxplot(aes(x = "Observed", y = CALAMT))
  p <- p + scale_x_discrete("Treatment Arm")
  p <- p + scale_y_log10("Residual (mmol)")
  p3 <- p + facet_wrap(~TRT)
  ggsave("output/UrinResVsTrt_EmaxOpt.png", width = 10.7, height = 8.03)
  
  
  p <- NULL
  p <- ggplot(data = urin_tb)
  p <- p + ggtitle("Urine Calcium Residuals")
  p <- p + geom_dotplot(aes(x = "Predicted", y = UCA), 
    stackdir = "center", binaxis = "y", dotsize = 0.5)
  p <- p + geom_dotplot(aes(x = "Observed", y = CALAMT), 
    stackdir = "center", binaxis = "y", dotsize = 0.5)
  p <- p + scale_x_discrete("Treatment Arm")
  p <- p + scale_y_log10("Residual (mmol)")
  p4 <- p + facet_wrap(~TRT)
  ggsave("output/UrinResVsTrtDotPlot_EmaxOpt.png", width = 10.7, height = 8.03)
