# Determine equation for inhibition of calcium reabsorption by caffeine
# -----------------------------------------------------------------------------
# The caffeine concentrations and there %effect on adenosine receptors (A1 &
#   A2A) are published in the following manuscript.

# This in vitro data will provide initial estimates for EC50 and Hill parameters
#   of an sigmoidal Emax equation that will ultimately be used to describe
#   the inhibition of calcium reabsorption by caffeine.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare work environment
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../calc_caff/")

# Document packages used via namespace (e.g. )
  # pkg_info <- c(
  #   "magrittr" = "1.5",  # use of pipes %>%
  #   "readr" = "1.3.1",  # reading in digitised data
  #   "ggplot2" = "3.1.1"  # graphical presentation of data
  # )
# Load packages (if not referred to via namespace i.e `package::function()`)
  library(magrittr)
  library(ggplot2)
  
# Source functions
  # source("scripts/functions_utility.R")

# Load data
  name_str <- c("Conc", "Eff")
  a1_tb <- readr::read_csv("raw_data/Magkos_A1.csv", col_names = name_str)
  a2a_tb <- readr::read_csv("raw_data/Magkos_A2A.csv", col_names = name_str)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Estimate sigmoidal emax equation parameters
# Emax
# As effect is a percentage...
  emax <- 100
  
# EC50 (uM)
# Both datasets have concentration for effect approximately 50%
  a1_ec50 <- a1_tb %>%
    dplyr::filter(abs(Eff - 50) == min(abs(Eff - 50))) %>%
    dplyr::pull(Conc) %>% signif(3)
  
  a2a_ec50 <- a2a_tb %>%
    dplyr::filter(abs(Eff - 50) == min(abs(Eff - 50))) %>%
    dplyr::pull(Conc) %>% signif(3)
  
# Hill
# Create function to predict effect
  sig_emax_fn <- function(conc, hill, ec50) {
    emax*conc^hill/(ec50^hill + conc^hill)
  }
  
# Create function to estimate hill parameter
  opt_hill_fn <- function(par, data) {
  # Get concentration and EC50 data
    ec50 <- get(paste0(data, "_ec50"))
    tb <- get(paste0(data, "_tb")) %>%
      dplyr::filter(abs(Conc - ec50) != min(abs(Conc - ec50)))
  # Predict effect
    yhat <- sig_emax_fn(tb$Conc, par[1], ec50)
  # Determine log-likelihood and calculate objective function value
    loglik <- dnorm(tb$Eff, yhat, par[2], log = T)
    return(-2*sum(loglik))
  }
  
# Estimate hill parameter value
  init_par <- c(hill = 1, sigma = 0.3)
  a1_est <- optim(init_par, opt_hill_fn, data = "a1", method = "L-BFGS-B",
    lower = 0, upper = 9,
    control = list(
    parscale = init_par, fnscale = opt_hill_fn(init_par, "a2a"),
    factr = 1e12, pgtol = 1e-8
  ))
  a1_hill <- signif(a1_est$par["hill"], 3)
  
  a2a_est <- optim(init_par, opt_hill_fn, data = "a2a", method = "L-BFGS-B",
    lower = 0, upper = 9,
    control = list(
    parscale = init_par, fnscale = opt_hill_fn(init_par, "a2a"),
    factr = 1e12, pgtol = 1e-8
  ))
  a2a_hill <- signif(a2a_est$par["hill"], 3)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Diagnostic plots
# Predict effect using estimated parameters
  sim_conc <- 10^c(0:32/8 - 0.65)
  a1_sim <- tibble::tibble(
    Conc = sim_conc,
    Eff = sig_emax_fn(sim_conc, a1_hill, a1_ec50)
  )
  a2a_sim <- tibble::tibble(
    Conc = sim_conc,
    Eff = sig_emax_fn(sim_conc, a2a_hill, a2a_ec50)
  )
  
# Create plots
  p <- NULL
  p <- ggplot()
  p <- p + ggtitle("A1 Adenosine Receptor")
  p <- p + geom_line(aes(x = Conc*0.19419, y = Eff), data = a1_sim,  # /
    size = 1, colour = "red")
  p <- p + geom_point(aes(x = Conc*0.19419, y = Eff), data = a1_tb, 
    size = 2, colour = "blue")
  p <- p + scale_x_log10("Caffeine Concentration (mg/L)", 
    labels = comma)
  p <- p + scale_y_continuous("Percent Effect (%)")
  p1 <- p + coord_cartesian(xlim = c(0.05, 300), ylim = c(0, 100))
  
  p <- NULL
  p <- ggplot()
  p <- p + ggtitle("A2A Adenosine Receptor")
  p <- p + geom_line(aes(x = Conc*0.19419, y = Eff), data = a2a_sim, 
    size = 1, colour = "red")
  p <- p + geom_point(aes(x = Conc*0.19419, y = Eff), data = a2a_tb, 
    size = 2, colour = "blue")
  p <- p + scale_x_log10("Caffeine Concentration (mg/L)", 
    labels = comma)
  p <- p + scale_y_continuous("Percent Effect (%)")
  p2 <- p + coord_cartesian(xlim = c(0.05, 300), ylim = c(0, 100))
  
  cowplot::plot_grid(p1, p2, align = "h", nrow = 1)
  