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
  library(quantreg)  # regression VPCs
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
  nsim <- 10
  trt <- with(dplyr::filter(pk_tb, !duplicated(ID)), table(TRT))
  uid <- sum(trt)  # Number of observed individuals
  nid <- uid*nsim  # Number of simulated individuals
  ID <- 1:nid  # Sequence of individual ID's
  
# Source population and add relevant info to population data
  source("scripts/190730_Population_Study.R")
  pop_tb <- dplyr::mutate(pop_tb, 
      ID_ORIG = rep(unique(pk_tb$ID_ORIG), nsim),
      SIM = rep(1:nsim, each = length(unique(pk_tb$ID_ORIG)))
    ) %>% dplyr::inner_join(trt_tb, by = "ID_ORIG")
  pop_tb <- dplyr::filter(pk_tb, !duplicated(ID)) %>%
    dplyr::select(ID_ORIG, P_0, GFR_0) %>%
    dplyr::inner_join(pop_tb, by = "ID_ORIG") %>%
    dplyr::arrange(ID) %>%
    dplyr::mutate(PDEMAX = 0.446)
# From `190731_Regimen_Study_EmaxOpt.R`:
# theta (se%): Emax = 0.446 (9.88 %); sigma = 0.760 (7.44 %)
# sigma is the variance of the proportional error (sd = 0.872)
  
# Read in simulation output
  output_tb <- readr::read_rds("output/EmaxSim.rds")
# Fix output from early sims to include ID_ORIG, SIM and TRT columns
  output_tb <- pk_tb %>%
    dplyr::filter(!duplicated(ID_ORIG)) %>%
    dplyr::select(ID_ORIG, TRT) %>%
    {tibble::tibble(
      ID = ID,
      ID_ORIG = rep(.$ID_ORIG, nsim),
      TRT = rep(.$TRT, nsim),
      SIM = rep(1:nsim, each = uid)
    )} %>%
    dplyr::inner_join(output_tb, by = "ID")
  
# Read in predicted values
  pred_tb <- readr::read_rds("output/EmaxPred.rds") %>%
    dplyr::mutate(TRT = dplyr::if_else(TRT == 1, "Caffeine", "Placebo"))
  
# Create predicted values for simulation and observed datasets
  pred_sim <- pred_tb %>% 
    dplyr::group_by(ID_ORIG) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(tb) {
    # Identify time where collection was started
      coll_tb <- dplyr::filter(tb, cmt == 31, evid == 0)
      time_val <- dplyr::pull(coll_tb, time)
      urin_val <- dplyr::pull(coll_tb, UCA)
    # For each one determine the dTIME and dUCA values
      day1_tb <- dplyr::filter(tb, time >= time_val[1] & time < time_val[2]) %>%
        dplyr::mutate(dTIME = time - time_val[1], dUCA = UCA - urin_val[1]) %>%
        dplyr::mutate(DAY = 1)
      day2_tb <- dplyr::filter(tb, time >= time_val[2]) %>%
        dplyr::mutate(dTIME = time - time_val[2], dUCA = UCA - urin_val[2]) %>%
        dplyr::mutate(DAY = 2)
      dplyr::bind_rows(day1_tb, day2_tb) %>%
        dplyr::select(dTIME, UCA_PRED = dUCA, DAY)
    })) %>% tidyr::unnest() %>% dplyr::group_by(ID_ORIG, DAY) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(tb) {
      dplyr::filter(tb, !duplicated(dTIME))
    })) %>% tidyr::unnest()
  
  pred_obs <- pred_tb %>%
    dplyr::filter(cmt == 31) %>%
    dplyr::group_by(ID_ORIG) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(tb) {
    # Identify time where collection was started
      coll_tb <- dplyr::filter(tb, cmt == 31, evid == 0)
      time_val <- dplyr::pull(coll_tb, time)
      urin_val <- dplyr::pull(coll_tb, UCA)
    # For each one determine the dTIME and dUCA values
      day1_tb <- dplyr::filter(tb, time >= time_val[1] & time < time_val[2]) %>%
        dplyr::mutate(dTIME = time - time_val[1], dUCA = UCA - urin_val[1]) %>%
        dplyr::mutate(DAY = 1)
      day2_tb <- dplyr::filter(tb, time >= time_val[2]) %>%
        dplyr::mutate(dTIME = time - time_val[2], dUCA = UCA - urin_val[2]) %>%
        dplyr::mutate(DAY = 2)
      dplyr::bind_rows(day1_tb, day2_tb) %>%
        dplyr::filter(evid == 1) %>%
        dplyr::select(dTIME, UCA_PRED = dUCA, DAY)
    })) %>% tidyr::unnest()
  
# Create simulation data
  plot_sim <- output_tb %>% 
    dplyr::select(ID, ID_ORIG, TRT, SIM, time, evid, amt, cmt, P, UCA) %>%
    dplyr::group_by(ID) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, function(tb) {
    # Identify time where collection was started
      coll_tb <- dplyr::filter(tb, cmt == 31, evid == 0)
      time_val <- dplyr::pull(coll_tb, time)
      urin_val <- dplyr::pull(coll_tb, UCA)
    # For each one determine the dTIME and dUCA values
      day1_tb <- dplyr::filter(tb, time > time_val[1] & time < time_val[2]) %>%
        dplyr::mutate(dTIME = time - time_val[1], dUCA = UCA - urin_val[1]) %>%
        dplyr::mutate(DAY = 1)
      day2_tb <- dplyr::filter(tb, time > time_val[2]) %>%
        dplyr::mutate(dTIME = time - time_val[2], dUCA = UCA - urin_val[2]) %>%
        dplyr::mutate(DAY = 2)
      dplyr::bind_rows(day1_tb, day2_tb)
    })) %>% tidyr::unnest() %>%
    dplyr::inner_join(pred_sim, by = c("ID_ORIG", "dTIME", "DAY")) %>%
    dplyr::rename(TIME = dTIME)
  caff_sim <- dplyr::filter(plot_sim, TRT == "Caffeine")
  plac_sim <- dplyr::filter(plot_sim, TRT == "Placebo")
  
# Create observed data
  plot_obs <- pk_tb %>%
    dplyr::select(ID, ID_ORIG, TRT, dTIME, evid, cmt, CALAMT) %>%
    dplyr::filter(cmt == 31 & evid == 1) %>%
    dplyr::mutate(DAY = rep(1:2, uid)) %>%
    dplyr::inner_join(pred_obs, by = c("ID_ORIG", "dTIME", "DAY")) %>%
    dplyr::rename(TIME = dTIME) %>% 
    dplyr::filter(!is.na(CALAMT))
  caff_obs <- dplyr::filter(plot_obs, TRT == "Caffeine")
  plac_obs <- dplyr::filter(plot_obs, TRT == "Placebo" & !is.na(CALAMT))
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Regression VPC (due to poor binning options)
# LOESS Smoothed PRED
# Define functions
  aicc.loess <- function(fit) { 
  # LOESS improved AIC function
  # compute AIC_C for a LOESS fit, from: 
  # Hurvich, C.M., Simonoff, J.S., and Tsai, C. L. 1998. Smoothing  
  # parameter selection in nonparametric regression using an improved  
  # Akaike Information Criterion. Journal of the Royal Statistical  
  # Society B 60: 271–293. 
  # Check that inputs are correct
    stopifnot(inherits(fit, 'loess')) 
    
  # calculate AIC_C
    n <- fit$n 
    trace <- fit$trace.hat 
    sigma2 <- sum(resid(fit)^2)/(n - 1) 
    return(log(sigma2) + 1 + (2*(trace + 1))/(n - trace - 2)) 
  } 

  autoloess <- function(fit, span = c(0.1, 0.9)) { 
  # LOESS span optimisation funhction
  # compute loess fit which has span minimizes AIC_C 
  # fit = loess fit; span parameter value doesn't matter 
  # span = a two-value vector representing the minimum and maximum span values 
  # Returns LOESS fit with span minimizing the AIC_C function 
  # check that input is correct
    stopifnot(inherits(fit, 'loess'), length(span) == 2) 
     
  # loess function in form to be used by optimize 
    f <- function(span) aicc.loess(update(fit, span = span)) 
     
  # find best loess according to loss function 
    opt.fit  <- update(fit, span = optimize(f, span)$minimum)
    opt.span <- optimize(f, span)$minimum
    return(list(fit = opt.fit, span = opt.span)) 
  } 
  
# Fit LOESS model and estimate PRED
# Caffeine
# Initially a LOESS model of PRED versus TIME is fit
  lmod <- loess(UCA_PRED ~ TIME, data = caff_obs, 
    span = 0.65, na.action = na.exclude)
  
# Then the span of the LOESS model is optimised
  lmod.opt <- autoloess(lmod)
  lmod.opt
  
# Compute estimated PRED and resulting prediction-correct values
  caff_lmod <- loess(UCA_PRED ~ TIME, data = caff_obs, 
    span = lmod.opt$span, na.action = na.exclude)
  caff_obs <- dplyr::mutate(caff_obs, UCA_ePRED = fitted(caff_lmod)) %>%
    dplyr::mutate(UCA_pcY = CALAMT*UCA_ePRED/UCA_PRED)
  
# Placebo
# Initially a LOESS model of PRED versus TIME is fit
  lmod <- loess(UCA_PRED ~ TIME, data = plac_obs, 
    span = 0.65, na.action = na.exclude)
  
# Then the span of the LOESS model is optimised
  lmod.opt <- autoloess(lmod)
  lmod.opt
  
# Compute estimated PRED and resulting prediction-correct values
  plac_lmod <- loess(UCA_PRED ~ TIME, data = plac_obs, 
    span = lmod.opt$span, na.action = na.exclude)
  plac_obs <- dplyr::mutate(plac_obs, UCA_ePRED = fitted(plac_lmod)) %>%
    dplyr::mutate(UCA_pcY = CALAMT*UCA_ePRED/UCA_PRED)
  
  # p05 <- NULL
  # p05 <- ggplot(data = plac_obs)
  # p05 <- p05 + geom_point(aes(x = TIME, y = UCA_PRED),
  #   colour = "blue", shape = 1)
  # p05 <- p05 + geom_point(aes(x = TIME, y = fitted(plac_lmod)),
  #   colour = "red", shape = 1)
  # p05 <- p05 + geom_point(aes(x = TIME, y = fitted(lmod)),
  #   colour = "green", shape = 1)
  # p05
  
# AQR characterised quantiles
# Define functions
  lambda.opt.pc <- function(llam, quant, data) {
    require(quantreg)
    a <- AIC(
      rqss(
        data$UCA_pcY ~ qss(data$TIME, lambda = exp(llam)), 
        tau = quant, na.action = na.exclude
      ), k = -1
    )
  }
  
# Caffeine
# Optimise smoothing parameter for each AQR model (one for each quantile)
  caff_llam_med <- optimise(lambda.opt.pc, quant = 0.50, data = caff_obs, 
    interval = c(0, 7))$min
  caff_llam_p05 <- optimise(lambda.opt.pc, quant = 0.05, data = caff_obs, 
    interval = c(0, 7))$min
  caff_llam_p95 <- optimise(lambda.opt.pc, quant = 0.95, data = caff_obs, 
    interval = c(0, 7))$min
  
# Fit AQR model to observed data
  caff_obs_med <- rqss(
    caff_obs$UCA_pcY ~ qss(caff_obs$TIME, lambda = exp(caff_llam_med)),
    tau = 0.50, na.action = na.exclude
  )
  caff_obs_p05 <- rqss(
    caff_obs$UCA_pcY ~ qss(caff_obs$TIME, lambda = exp(caff_llam_p05)),
    tau = 0.05, na.action = na.exclude
  )
  caff_obs_p95 <- rqss(
    caff_obs$UCA_pcY ~ qss(caff_obs$TIME, lambda = exp(caff_llam_p95)),
    tau = 0.95, na.action = na.exclude
  )
  
# Add fitted AQR data to observed dataset
  caff_obs <- dplyr::mutate(caff_obs,
    med_fit = fitted(caff_obs_med),
    p05_fit = fitted(caff_obs_p05),
    p95_fit = fitted(caff_obs_p95)
  )
  
# Placebo
# Optimise smoothing parameter for each AQR model (one for each quantile)
  plac_llam_med <- optimise(lambda.opt.pc, quant = 0.50, plac_obs, 
    interval = c(0, 7))$min
  plac_llam_p05 <- optimise(lambda.opt.pc, quant = 0.05, plac_obs, 
    interval = c(0, 7))$min
  plac_llam_p95 <- optimise(lambda.opt.pc, quant = 0.95, plac_obs, 
    interval = c(0, 7))$min
  
# Fit AQR model to observed data
  plac_obs_med <- rqss(
    plac_obs$UCA_pcY ~ qss(plac_obs$TIME, lambda = exp(plac_llam_med)),
    tau = 0.50, na.action = na.exclude
  )
  plac_obs_p05 <- rqss(
    plac_obs$UCA_pcY ~ qss(plac_obs$TIME, lambda = exp(plac_llam_p05)),
    tau = 0.05, na.action = na.exclude
  )
  plac_obs_p95 <- rqss(
    plac_obs$UCA_pcY ~ qss(plac_obs$TIME, lambda = exp(plac_llam_p95)),
    tau = 0.95, na.action = na.exclude
  )
  
# Add fitted AQR data to observed dataset
  plac_obs <- dplyr::mutate(plac_obs,
    med_fit = fitted(plac_obs_med),
    p05_fit = fitted(plac_obs_p05),
    p95_fit = fitted(plac_obs_p95)
  )
  
# Bind obs data together
  plot_obs <- rbind(caff_obs, plac_obs)
  
# Check optimised values
  # p06 <- NULL
  # p06 <- ggplot(data = plac_obs)
  # p06 <- p06 + geom_point(aes(x = TIME, y = UCA_pcY),
  #   colour = "blue", shape = 1)
  # p06 <- p06 + geom_line(aes(x = TIME, y = med_fit), size = 0.8)
  # p06 <- p06 + geom_line(aes(x = TIME, y = p05_fit),
  #   linetype = "dashed", size = 0.8)
  # p06 <- p06 + geom_line(aes(x = TIME, y = p95_fit),
  #   linetype = "dashed", size = 0.8)
  # p06
  
# Use ddply to apply LOESS and AQR smooths to simulated datasets
  caff_sim <- plyr::ddply(caff_sim, plyr::.(SIM), function(df) {
  # Calculate ePREDs from fitted LOESS models
    pred.loess <- loess(
      UCA_PRED ~ TIME, data = df, span = caff_lmod$par$span, na.action = na.exclude
    )
    df$UCA_ePRED <- fitted(pred.loess)
    
  # Calculate prediction-corrected values
    df$UCA_pcY <- df$dUCA*(df$UCA_ePRED/df$UCA_PRED)
    df <- df[df$UCA_pcY >= 0, ]
    day1_df <- df[df$DAY == 1, ]
    day2_df <- df[df$DAY == 2, ]
    
  # Fit AQR models
    day1_med <- rqss(
      day1_df$UCA_pcY ~ qss(day1_df$TIME, lambda = exp(caff_llam_med)),
      tau = 0.50, na.action = na.exclude
    )
    day2_med <- rqss(
      day2_df$UCA_pcY ~ qss(day2_df$TIME, lambda = exp(caff_llam_med)),
      tau = 0.50, na.action = na.exclude
    )
    day1_p05 <- rqss(
      day1_df$UCA_pcY ~ qss(day1_df$TIME, lambda = exp(caff_llam_p05)),
      tau = 0.05, na.action = na.exclude
    )
    day2_p05 <- rqss(
      day2_df$UCA_pcY ~ qss(day2_df$TIME, lambda = exp(caff_llam_p05)),
      tau = 0.05, na.action = na.exclude
    )
    day1_p95 <- rqss(
      day1_df$UCA_pcY ~ qss(day1_df$TIME, lambda = exp(caff_llam_p95)),
      tau = 0.95, na.action = na.exclude
    )
    day2_p95 <- rqss(
      day2_df$UCA_pcY ~ qss(day2_df$TIME, lambda = exp(caff_llam_p95)),
      tau = 0.95, na.action = na.exclude
    )
    
  # Collect fitted data from AQR models
    cbind(df,
      data.frame(
        med_fit = c(fitted(day1_med), fitted(day2_med)),
        p05_fit = c(fitted(day1_p05), fitted(day2_p05)),
        p95_fit = c(fitted(day1_p95), fitted(day2_p95))
      )
    )
  })
  plac_sim <- plyr::ddply(plac_sim, plyr::.(SIM), function(df) {
# Calculate ePREDs from fitted LOESS models
    pred.loess <- loess(
      UCA_PRED ~ TIME, data = df, span = plac_lmod$par$span, na.action = na.exclude
    )
    df$UCA_ePRED <- fitted(pred.loess)
    
  # Calculate prediction-corrected values
    df$UCA_pcY <- df$dUCA*(df$UCA_ePRED/df$UCA_PRED)
    df <- df[df$UCA_pcY >= 0, ]
    day1_df <- df[df$DAY == 1, ]
    day2_df <- df[df$DAY == 2, ]
    
  # Fit AQR models
    day1_med <- rqss(
      day1_df$UCA_pcY ~ qss(day1_df$TIME, lambda = exp(caff_llam_med)),
      tau = 0.50, na.action = na.exclude
    )
    day2_med <- rqss(
      day2_df$UCA_pcY ~ qss(day2_df$TIME, lambda = exp(caff_llam_med)),
      tau = 0.50, na.action = na.exclude
    )
    day1_p05 <- rqss(
      day1_df$UCA_pcY ~ qss(day1_df$TIME, lambda = exp(caff_llam_p05)),
      tau = 0.05, na.action = na.exclude
    )
    day2_p05 <- rqss(
      day2_df$UCA_pcY ~ qss(day2_df$TIME, lambda = exp(caff_llam_p05)),
      tau = 0.05, na.action = na.exclude
    )
    day1_p95 <- rqss(
      day1_df$UCA_pcY ~ qss(day1_df$TIME, lambda = exp(caff_llam_p95)),
      tau = 0.95, na.action = na.exclude
    )
    day2_p95 <- rqss(
      day2_df$UCA_pcY ~ qss(day2_df$TIME, lambda = exp(caff_llam_p95)),
      tau = 0.95, na.action = na.exclude
    )
    
  # Collect fitted data from AQR models
    cbind(df,
      data.frame(
        med_fit = c(fitted(day1_med), fitted(day2_med)),
        p05_fit = c(fitted(day1_p05), fitted(day2_p05)),
        p95_fit = c(fitted(day1_p95), fitted(day2_p95))
      )
    )
  })
  plot_sim <- rbind(caff_sim, plac_sim)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# VPC
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
  p07 <- NULL
  p07 <- ggplot(plot_sim)
  p07 <- p07 + geom_point(aes(x = TIME, y = UCA_pcY), data = plot_obs,
    colour = "blue", shape = 1)
  
  p07 <- p07 + geom_smooth(aes(x = TIME, y = med_fit), method = "loess", 
    colour = "black", size = 0.8)
  p07 <- p07 + geom_smooth(aes(x = TIME, y = p05_fit), method = "loess", 
    colour = "red", size = 0.8, linetype = "dashed", alpha = 0.5)
  p07 <- p07 + geom_smooth(aes(x = TIME, y = p95_fit), method = "loess", 
    colour = "red", size = 0.8, linetype = "dashed", alpha = 0.5)
  
  p07 <- p07 + facet_wrap(~TRT)
  
  p07 <- p07 + scale_x_continuous("Time (hours)", breaks = 0:6*8)
  p07 <- p07 + scale_y_log10("Calcium (mmol)")
  
  p07
  ggsave("output/urine_vpc_varPGFR_EmaxOpt.png", width = 10.8, height = 8.7)
  