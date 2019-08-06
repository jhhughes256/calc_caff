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
  nsim <- 1000
  trt <- with(dplyr::filter(pk_tb, !duplicated(ID)), table(TRT))
  uid <- sum(trt)  # Number of observed individuals
  nid <- uid*nsim  # Number of simulated individuals
  ID <- 1:nid  # Sequence of individual ID's
  
# Source population and add relevant info to population data
  source("scripts/190730_Population_Study.R")
  pop_tb <- dplyr::mutate(pop_tb, 
    ID_ORIG = rep(unique(pk_tb$ID_ORIG), nsim),
    SIM = rep(1:nsim, each = length(unique(pk_tb$ID_ORIG))))
# TODO: Incorporate TRT here
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
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Visual Plots
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
# Prepare data
# Define VPC bins
  bin_cuts <- c(0, 2, 5, 8, 12, 26, 29, 32, 35)
  
# Plasma
  calamt_plas <- pk_tb %>%
    dplyr::filter(cmt == 6 & !is.na(time)) %>%
    dplyr::pull(CALAMT)
  plas_tb <- output_tb %>%
    dplyr::filter(cmt == 6) %>%
    dplyr::mutate(CALAMT = rep(calamt_plas, nsim)) %>%
    dplyr::mutate(RES = rep(calamt_plas, nsim) - P) %>%
    dplyr::mutate(SIM = rep(1:nsim, each = length(calamt_plas))) %>%
    dplyr::mutate(TADBIN = Hmisc::cut2(time, cuts = bin_cuts, levels.mean = T)) %>%
    dplyr::mutate(TADBIN = as.numeric(paste(TADBIN)))
  
  vpc_plas_pk <- pk_tb %>%
    dplyr::filter(cmt == 6 & !is.na(time)) %>%
    dplyr::mutate(TADBIN = Hmisc::cut2(time, cuts = bin_cuts, levels.mean = T)) %>%
    dplyr::mutate(TADBIN = as.numeric(paste(TADBIN)))
  
  vpc_plas_tb <- plyr::ddply(plas_tb, plyr::.(SIM, TADBIN, TRT), function(x) {
    data.frame(
      medianS = median(x$P),
      loCI90S = CI90lo(x$P),
      hiCI90S = CI90hi(x$P)
    )
  })
  
# Urine
  calamt_urin <- pk_tb %>%
    dplyr::filter(cmt == 31 & evid == 1) %>%
    dplyr::pull(CALAMT)
  uca_urin <- output_tb %>%
    dplyr::filter(cmt == 31) %>%
    {dplyr::filter(., evid == 1)$UCA - dplyr::filter(., evid == 0)$UCA}
  uca_dtime <- output_tb %>%
    dplyr::filter(cmt == 31) %>%
    {dplyr::filter(., evid == 1)$time - dplyr::filter(., evid == 0)$time}
  urin_tb <- output_tb %>%
    dplyr::filter(cmt == 31 & evid == 1) %>%
    dplyr::mutate(UCA = uca_urin, CALAMT = rep(calamt_urin, nsim)) %>%
    dplyr::mutate(dTIME = uca_dtime, RES = CALAMT - UCA) %>%
    dplyr::mutate(SIM = rep(1:nsim, each = length(calamt_urin)))
    
  
  vpc_urin_tb <- plyr::ddply(urin_tb, plyr::.(SIM, dTIME, TRT), function(x) {
    data.frame(
      medianS = median(x$UCA),
      loCI90S = CI90lo(x$UCA),
      hiCI90S = CI90hi(x$UCA)
    )
  })
  
# Plasma Calcium VPC
  p <- NULL
	p <- ggplot(data = vpc_plas_pk)
	
  p <- p + geom_point(aes(x = TADBIN, y = CALAMT), colour = "blue", shape = 1)
  p <- p + stat_summary(aes(x = TADBIN, y = CALAMT), fun.y = median,
    geom = "line", colour = "red", size = 1)
  p <- p + stat_summary(aes(x = TADBIN, y = CALAMT), fun.y = CI90lo,
    geom = "line", colour = "red", linetype = "dashed", size = 1)
  p <- p + stat_summary(aes(x = TADBIN, y = CALAMT), fun.y = CI90hi,
    geom = "line", colour = "red", linetype = "dashed", size = 1)

	p <- p + stat_summary(aes(x = TADBIN, y = medianS), data = vpc_plas_tb,
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "red")
  p <- p + stat_summary(aes(x = TADBIN, y = medianS), data = vpc_plas_tb,
    fun.y = median, geom = "line", colour = "black", size = 1)

	p <- p + stat_summary(aes(x = TADBIN, y = loCI90S), data = vpc_plas_tb,
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
	p <- p + stat_summary(aes(x = TADBIN, y = loCI90S), data = vpc_plas_tb,
    fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

	p <- p + stat_summary(aes(x = TADBIN, y = hiCI90S), data = vpc_plas_tb,
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
	p <- p + stat_summary(aes(x = TADBIN, y = hiCI90S), data = vpc_plas_tb,
    fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)
	
	p <- p + facet_wrap(~TRT)
	p
	
	
  p <- NULL
	p <- ggplot(data = dplyr::filter(pk_tb, cmt == 31 & evid == 1))
	
  p <- p + geom_point(aes(x = dTIME, y = CALAMT), colour = "blue", shape = 1)
  p <- p + stat_summary(aes(x = dTIME, y = CALAMT), fun.y = median,
    geom = "line", colour = "red", size = 1)
  p <- p + stat_summary(aes(x = dTIME, y = CALAMT), fun.y = CI90lo,
    geom = "line", colour = "red", linetype = "dashed", size = 1)
  p <- p + stat_summary(aes(x = dTIME, y = CALAMT), fun.y = CI90hi,
    geom = "line", colour = "red", linetype = "dashed", size = 1)

	p <- p + stat_summary(aes(x = dTIME, y = medianS), data = vpc_urin_tb,
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "red")
  p <- p + stat_summary(aes(x = dTIME, y = medianS), data = vpc_urin_tb,
    fun.y = median, geom = "line", colour = "black", size = 1)

	p <- p + stat_summary(aes(x = dTIME, y = loCI90S), data = vpc_urin_tb,
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
	p <- p + stat_summary(aes(x = dTIME, y = loCI90S), data = vpc_urin_tb,
    fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)

	p <- p + stat_summary(aes(x = dTIME, y = hiCI90S), data = vpc_urin_tb,
    geom = "ribbon", fun.ymin = "CI95lo", fun.ymax = "CI95hi", alpha = 0.3, fill = "blue")
	p <- p + stat_summary(aes(x = dTIME, y = hiCI90S), data = vpc_urin_tb,
    fun.y = median, geom = "line", colour = "black", linetype = "dashed", size = 1)
	
	p <- p + coord_trans(y = "log10", limy = c(0.4, 25))
	p <- p + facet_wrap(~TRT)
	p
	
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
