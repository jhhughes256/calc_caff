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
  # library(readxl)

# Read in data
  file_path <- "raw_data/20160707 Caffeine Study Database_unfiltered.xlsx"
  pk_tb <- readxl::read_excel(file_path, sheet = "PK", na = "NR", skip = 1,
    col_names = c(
      "ID", "TRT", "MAT", "SAMPLE", "DATST", "DATEN", "VOL", 
      "SOD", "SECR", "CAL", "NACL", "CRCL", "CACL"
    ),
    col_types = c(
      "numeric", "text", "text", "numeric", "date", "date", rep("numeric", 7)
    )
  )
  
    