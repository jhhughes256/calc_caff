# Testing the use of Metrum's mrgsolve model of the calcium bone model
# -----------------------------------------------------------------------------
# Metrum provides an implementation of Peterson and Riggs Calcium Homeostasis 
#   and Bone Remodelling QSP Model in mrgsolve in their GitHub repository:
#   https://github.com/metrumresearchgroup/OpenBoneMin. Model originally from:
#
# M. C. Peterson and M. M. Riggs. A physiologically based mathematical model of 
#   integrated calcium homeostasis and bone remodeling. Bone, 46:49–63, Jan 2010
#
# Metrum model is different in that they have implemented denosumab into the
#   model alongside teriparatide. Denosumab model based orginally from:
#
# M. Peterson B. Stouch D. Chen S. Baughman D. Holloway P. Bekker and S. Martin.
#   A population pk/pd model describes the rapid profound and sustained 
#   suppression of urinary n-telopeptide following administration of amg 162 a 
#   fully human monoclonal antibody against rankl to healthy postmenopausal 
#   women. The AAPS Journal 24(6 Abstract W4340) 2004.
#
# Script seeks to load the model and explore it's use, with the aim of linking
#   it with a caffeine PKPD.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Prepare workspace
# Clear workspace
  rm(list=ls(all=TRUE))
  graphics.off()

# Set working directory
# if not working with RStudio project place working directory here
# setwd("C:/.../calc_caff/")

# # Load package libraries
  # library(dplyr)	# dplyr required for mrgsolve
  # library(mrgsolve)  # Metrum differential equation solver for pharmacometrics

# Source external scripts

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Explore calcium bone model
# Read in model
  mod <- mrgsolve::mread("model/cabone_base.cpp")
  
# Using init we can observe the compartments
  mrgsolve::init(mod)
# Compartments are named according to the README by Peterson and Riggs with the
#   following exceptions:
# A -> AOH
# SC -> ? (subcutaneous PTH compartment appears to be gone?)
# Additional compartments are:
# DENSC, DENCENT, DENPER - added for denos model, PK compartments (depot, central and peripheral)
# BMDfn, BMDls, BMDlsDEN - added for denos model, seems to calculate BMD on the fly
# EST, GFR - added for denos model
# PKGUT, PKCENT, PKPER1, PKPER2 - added as a general PK compartment (we can put caffeine in here!)
# UCA - collects calcium in urine

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Try to simulate with this model using teriparatide
# Specify dosing compartment and dose
  cmt_n <- mrgsolve::cmtn(mod, "TERISC")
  teri_amt <- c(20, 40)*1E6/4117.8  # convert to dose units (ug -> pmol)

# Set up a quick simulation template
  input_df <- mrgsolve::expand.ev(  # quick data template
    amt = teri_amt, ii = 24, dur = 9, cmt = cmt_n, addl = 10
  )
  
# Set requested data from model and simulate
  req_out <- "PTHpm,CaC"  # define what we want out of the model
  output_df <- mrgsolve::mrgsim(mod, 
    data = input_df, delta = 0.1, end = (9+1)*24, Req = req_out
  )
  
# Plot output
  mrgsolve::plot(output_df)
  
# Exploring further...
# Can also request the values in any of the compartments
  req_out <- "TERISC,UCA"  # define what we want out of the model
  output_df <- mrgsolve::mrgsim(mod, 
    data = input_df, delta = 0.1, end = (9+1)*24, Req = req_out
  )
  mrgsolve::plot(output_df)
# This can be used to capture PK predictions as well as calcium excretion
# For caffeine would compare UCA for placebo vs treatment
  