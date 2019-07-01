# Caffeine Population Pharmacokinetic Model - Perera et al.
# ------------------------------------------------------------------------------
# Parameters taken from 24 hours methylxanthine abstinence model
#
# Model sourced from:
# Perera V, Gross AS, Forrest A, et al. A Pharmacometric Approach to Investigate
#   the Impact of Methylxanthine Abstinence and Caffeine Consumption on CYP1A2 
#   Activity. Drug Metab Disposition. 2013;41(11):1957-1966.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load libraries
  # library(dplyr)
  # library(mrgsolve)

# Define model code
  code <- '
$INIT  // Initial Conditions for Compartments (mg)
  GUT = 1.37,  // Absorption Compartment // IC(4) - 1.37 FIXED
  CENT = 1.72,  // Central Compartment // IC(1) - 1.72 (229 CV%)
  PERI = 1.20,  // Peripheral Compartment // IC(2) - 1.20 FIXED
  PARA = 4.06,  // Paraxanthine Compartment // IC(3) - 4.06 (48.2 CV%)

$SET  // Set Differential Equation Solver Options			
  atol = 1e-8, rtol = 1e-8
  maxsteps = 100000

$PARAM  // Population parameters
  TVKA = 3.65,  // absorption rate constant (h-1)
  TVK12 = 1.71,  // microconstant - cent-peri (h-1)
  TVK21 = 2.02,  // microconstant - peri-cent (h-1)
  TVV1 = 31.2,  // central volume of distribution (L)
  TVV2 = 24.4,  // peripheral volume of distribution (L)
  TVTLAG = 0.15,  // lag-time (h) {says h-1 in manuscript?}
  TVCLCAPX = 2.73,  // clearance of caffeine - to paraxanthine (L/h)
  TVCLCAO = 4.86,  // clearance of caffeine - renal and other (L/h)
  TVCLPX = 3.81,  // clearance of paraxanthine (L/h)
  LLOQ = 0.05,  // assay limit of quantification (mg/L) {used as additive RUV}

  // Default Covariate Values for Simulation

  // Default ETA Values for Simulation
  // Allocated in population so set to zero
  ETA1 = 0,  // PPVKA
  ETA2 = 0,  // PPVK12
  ETA3 = 0,  // PPVK21
  ETA4 = 0,  // PPVV1
  ETA5 = 0,  // PPVV2
  ETA6 = 0,  // PPVTLAG
  ETA7 = 0,  // PPVCLCAPX
  ETA8 = 0,  // PPVCLCAO
  ETA9 = 0,  // PPVCLPX

$OMEGA  // Population parameter Variability
  name = "omega1"
  block = FALSE
  1.060900  // PPVKA
  0.105625  // PPVK12
  0.053361  // PPVK21
  0.031329  // PPVV1
  0.060516  // PPVV2
  1.144900  // PPVLAG
  0.051529  // PPVCLCAPX
  0.326041  // PPVCLCAO
  0.017956  // PPVCLPX
  5.244100  // PPVCENT
  0.232324  // PPVPARA
  
$SIGMA  // Residual Unexplained Variability	
  block = FALSE
  0.02  // RUVPROCA // SD1slope
  0.04  // RUVPROPX // SD2slope

$MAIN  // Individual Parameter Values
  double KA = TVKA*exp(ETA1);  // *exp(ETA(1))
  double K12 = TVK12*exp(ETA2);  // *exp(ETA(2))
  double K21 = TVK21*exp(ETA3);  // *exp(ETA(3))
  double V1 = TVV1*exp(ETA4);  // *exp(ETA(4))
  double V2 = TVV2*exp(ETA5);  // *exp(ETA(5))
  double TLAG = TVTLAG*exp(ETA6);  // *exp(ETA(6))
  double CLCAPX = TVCLCAPX*exp(ETA7);  // *exp(ETA(7))
  double CLCAO = TVCLCAO*exp(ETA8);  // *exp(ETA(8))
  double CLPX = TVCLPX*exp(ETA9);  // *exp(ETA(9))
  ALAG_GUT = TLAG;

$ODE  // Differential Equations
  double C1 = CENT/V1;
  double C2 = PARA/V2;

  dxdt_GUT = -KA*GUT;
  dxdt_CENT = KA*GUT - CENT*K12 + PERI*K21 - C1*CLCAPX - C1*CLCAO;
  dxdt_PERI = CENT*K12 - PERI*K21;
  dxdt_PARA = C1*CLCAPX - C2*CLPX;

$TABLE  // Determines Values and Includes in Output	
  double IPRECA = C1;  // caffeine individual prediction
  double IPREPX = C2;  // paraxanthine individual prediction

  double DVCA = IPRECA*(1 + EPS(1)) + LLOQ;  // caffeine observed conc
  double DVPX = IPREPX*(1 + EPS(2)) + LLOQ;  // paraxanthine observed conc

$CAPTURE 
  IPRECA IPREPX DVCA DVPX GUT CENT PERI PARA 
  KA K12 K21 V1 V2 TLAG CLCAPX CLCAO CLPX
  ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 ETA7 ETA8 ETA9 EPS(1) EPS(2)
'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile the model code
  mod <- mrgsolve::mcode("CaffPK", code)