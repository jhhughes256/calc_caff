# Caffine PopPK Model Simulation - Population
# ------------------------------------------------------------------------------
# Create population of representative patients for simulations
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define demographic data
  set.seed(123456) 
  
# Population parameter variability
# Extract omega block values from model and allocate individual distributions
  ETA_tb <- mrgsolve::omat(capx_mod, make = TRUE) %>%
    {MASS::mvrnorm(n = nid, mu = rep(0, dim(.)[1]), Sigma = .)} %>%
    tibble::as_tibble() %>% 
    dplyr::rename_all(function(x) paste0("ETA", readr::parse_number(x))) %>%
    dplyr::rename(PKCENT_0 = ETA10, PKPARA_0 = ETA11) %>%
    dplyr::mutate(PKCENT_0 = init(capx_mod)[["CENT"]]*exp(PKCENT_0)) %>%
    dplyr::mutate(PKPARA_0 = init(capx_mod)[["PARA"]]*exp(PKPARA_0))
  
# Create data frame of individuals with varying demographics and ETA values
  pop_tb <- dplyr::bind_cols(ID = ID, ETA_tb)
  