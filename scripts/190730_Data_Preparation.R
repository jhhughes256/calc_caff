# Prepare calcium data for use
# -----------------------------------------------------------------------------
# Data previously published anywhere? Find out!

# Data needs to be prepared such that there are separate datasets for calcium
#   concentrations in plasma and urine. To be compared to results from QSP 
#   model, when linked with caffeine PK model by a Imax model on calcium 
#   reabsorption.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare work environment
# Load libraries - documentation purposes only, libraries don't need loading
  library(lubridate)  # required for its changes to how dates work
  # library(readxl)
  # library(dplyr)
  # library(tidyr)
  # library(purrr)

# Read in data
  file_path <- "raw_data/20160707_Caffeine_Study_Database_unfiltered.xlsx"
  raw_tb <- readxl::read_excel(file_path, sheet = "PK", na = "NR", skip = 1,
    col_names = c(
      "ID", "TRT", "MAT", "SAMPLE", "DATST", "DATEN", "VOL", 
      "SOD", "SECR", "CAL", "NACL", "CRCL", "CACL"
    ),
    col_types = c(
      "numeric", "text", "text", "numeric", "date", "date", rep("numeric", 7)
    )
  )
  
# Fix small errors in the dataset
# Incorrect month entry for blood draw (time checked against protocol)
  raw_tb <- dplyr::mutate(raw_tb, DATST = dplyr::if_else(
    ID == 1423 & SAMPLE == 14, 
    DATST - months(1), DATST)
  )
  
# Add in time after earliest sample/dose per patient
  mutate_data <- function(tb) {
  # Set reference time point to be 11pm on first day of study for that patient
    first_time <- as.Date(min(tb$DATST, na.rm = TRUE))
    lubridate::hour(first_time) <- 23
  # Calculate time after ref. time and add in other important mrgsolve columns
    dplyr::mutate(tb, time = as.numeric(DATST - first_time, units = "hours")) %>%
      dplyr::mutate(cmt = dplyr::if_else(MAT == "Plasma", 6, 31)) %>%
      dplyr::mutate(amt = 0, evid = 0, rate = 0) %>%
  # Calculate time after ref. time that urine collection ended
      dplyr::mutate(TIMEEND = as.numeric(DATEN - first_time, units = "hours")) %>%
      dplyr::mutate(dTIME = TIMEEND - time) %>%
  # Calculate patients average GFR
      dplyr::mutate(GFR_0 = mean(CRCL, na.rm = T)) %>%
      dplyr::mutate(GFR_0 = dplyr::if_else(GFR_0 > 7.2, 7.2, GFR_0)) %>%
  # Convert calcium conc to amount and determine initial calcium plasma conc
      dplyr::mutate(CALAMT = CAL*dplyr::if_else(!is.na(VOL), VOL, 14)) %>%
      dplyr::mutate(P_0 = head(CALAMT[MAT == "Plasma"], 1))
  }
  obs_tb <- dplyr::arrange(raw_tb, ID, DATST) %>%
    dplyr::group_by(ID) %>% tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, mutate_data)) %>%
    tidyr::unnest() %>%
    dplyr::rename(ID_ORIG = ID) %>%
    dplyr::mutate(ID = rep(1:length(unique(ID_ORIG)), each = 10))
# Warnings due to NAs
  
# Define additional samples for urine sample time ending
  urine_tb <- obs_tb %>%
    dplyr::filter(!is.na(TIMEEND)) %>%
    dplyr::mutate(time = TIMEEND, evid = 1)
    
  
# Define dosing data (when made available)
  dose_times <- c(2, 4, 6, 8, 26, 28, 30, 32)
  dose_uid <- length(unique(obs_tb$ID))
  dose_tb <- dplyr::filter(obs_tb, !duplicated(ID)) %>%
    dplyr::mutate_at(
      c("SAMPLE", "DATST", "DATEN", "VOL", "SOD", "SECR", 
        "CAL", "NACL", "CRCL", "CACL", "TIMEEND", "CALAMT"),
      function(x) NA
    ) %>% 
    dplyr::slice(rep(1:dose_uid, each = 8)) %>%  # rep.data.frame
    dplyr::mutate(MAT = "Dose", cmt = 25, evid = 1, rate = 0) %>%
    dplyr::mutate(time = rep(dose_times, times = dose_uid)) %>%
    dplyr::mutate(amt = dplyr::if_else(TRT == "Caffeine", 200, 0))
  
# Create final pk dataset
  pk_tb <- dplyr::bind_rows(obs_tb, dose_tb, urine_tb) %>%
    dplyr::arrange(ID, time, dplyr::desc(amt))
  
# Identify caffeine and plasma ID's
  id_caff <- dplyr::filter(pk_tb, TRT == "Caffeine" & !duplicated(ID)) %>% 
    dplyr::pull(ID)
  id_plac <- dplyr::filter(pk_tb, TRT == "Placebo" & !duplicated(ID)) %>% 
    dplyr::pull(ID)
  