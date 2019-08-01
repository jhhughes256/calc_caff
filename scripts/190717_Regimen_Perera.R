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
  # library(MASS)  # mvrnorm

# Source external scripts
  source("scripts/functions_utility.R")  # functions utility
  source("model/caffpk_perera_MA.R")  # PopPK model script
  # source("model/caffpk_perera_noMA.R")  # PopPK model script

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Number of individuals
  nsim <- 1000
  nid <- 30*nsim  # Number of individuals
  ID <- 1:nid  # Sequence of individual ID's
  
# Source population
  source("scripts/190628_Population_Perera.R")
  
# Create simulation input dataset
# Define time points
  conc_times <- c(0, 0.5, 1, 1.5, 2, 4, 6, 8, 10, 24, 144, 168, 192)
  dose_times <- 0  # 1 dose
  dose_amt <- 100
  nobs <- length(conc_times)
  
# Duplicate population data for number of times and alter dose times to include
#   dosage events for mrgsolve
  concdose_tb <- tibble::tibble(
    ID = rep(ID, each = length(conc_times)),
    time = rep(conc_times, nid),
    amt = dplyr::if_else(time %in% dose_times, dose_amt, 0),
    cmt = 1,
    evid = dplyr::if_else(amt != 0, 1, 0),
    rate = 0,
    SIM = rep(1:1000, each = length(conc_times)*nid/1000)
  )
  
# Residual unexplained variability
  EPS_tb <- mrgsolve::smat(capx_mod, make = TRUE) %>%
    {MASS::mvrnorm(n = nid*nobs, mu = rep(0, dim(.)[1]), Sigma = .)} %>%
    tibble::as_tibble() %>% 
    dplyr::rename_all(function(x) paste0("EPS", readr::parse_number(x)))
  
# Create final input dataset
  input_tb <- dplyr::bind_cols(concdose_tb, EPS_tb)
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulate standard dosing regimen
# Pipe dataset to model
  output_tb <- capx_mod %>%
    mrgsolve::data_set(input_tb) %>%  # set input data (observations/events)
    mrgsolve::idata_set(pop_tb) %>%  # set individual data (ETAs/ICs)
    mrgsolve::carry_out(amt, evid, rate, cmt, SIM) %>%  # copy to simulated output
    mrgsolve::mrgsim() %>%  # simulate using mrgsolve
    tibble::as_tibble()
  
# Save data as .RDS for later use
  # saveRDS(output_tb, "output/pkdata_perera.rds")
  
# Observe population parameters
  output_tb %>% dplyr::select(
      ID, time, evid, amt, cmt, SIM, KA, K12, K21, 
      V1, V2, TLAG, CLCAPX, CLCAO, CLPX
    ) %>% dplyr::filter(!duplicated(ID))
  
# Observe mean concentration after 7 days abstinence
  if (168 %in% conc_times) {
    output_tb %>% filter(time == 168) %>% pull(DVCA) %>% mean()
    output_tb %>% filter(time == 168) %>% pull(DVPX) %>% mean()
  }
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot simulated component of VPC
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))

# Calculate 5, 50 and 95 percentiles for each simulated study
  sim_vpc_prep <- function(tb) {
    with(tb, tibble::tibble(
      analyte = c("Caffeine", "Paraxanthine"),
      medianS = c(
        median(DVCA, na.rm = T),
        median(DVPX, na.rm = T)
      ),
      CI90loS = c(
        quantile(DVCA, probs = 0.05, na.rm = T),
        quantile(DVPX, probs = 0.05, na.rm = T)
      ),
      CI90hiS = c(
        quantile(DVCA, probs = 0.95, na.rm = T),
        quantile(DVPX, probs = 0.95, na.rm = T) 
      )
    ))
  }
  sim_tb <- output_tb %>%
    dplyr::group_by(SIM, time) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, sim_vpc_prep)) %>%
    tidyr::unnest()
  
  p <- NULL
  p <- ggplot(data = dplyr::filter(sim_tb, analyte == "Caffeine"))
  p <- p + stat_summary(aes(x = time, y = medianS), geom = "line", 
    fun.y = median, size = 1)
  p <- p + stat_summary(aes(x = time, y = CI90loS), geom = "line", 
    fun.y = median, size = 1, colour = "red", alpha = 0.3)
  p <- p + stat_summary(aes(x = time, y = CI90hiS), geom = "line", 
    fun.y = median, size = 1, colour = "red", alpha = 0.3)
  p <- p + labs(x = "Time (hours)", y = "Concentration (ug/mL)")
  p1 <- p + coord_cartesian(xlim = c(0, 24), ylim = c(0, 4))  # ylim = c(0, 3))
  
  p <- NULL
  p <- ggplot(data = dplyr::filter(sim_tb, analyte == "Paraxanthine"))
  p <- p + stat_summary(aes(x = time, y = medianS), geom = "line", 
    fun.y = median, size = 1)
  p <- p + stat_summary(aes(x = time, y = CI90loS), geom = "line", 
    fun.y = median, size = 1, colour = "red", alpha = 0.3)
  p <- p + stat_summary(aes(x = time, y = CI90hiS), geom = "line", 
    fun.y = median, size = 1, colour = "red", alpha = 0.3)
  p <- p + labs(x = "Time (hours)", y = "Concentration (ug/mL)")
  p2 <- p + coord_cartesian(xlim = c(0, 25), ylim = c(0, 2))  # ylim = c(0, 0.8))
  
  cowplot::plot_grid(p1, p2, nrow = 1, align = "h")