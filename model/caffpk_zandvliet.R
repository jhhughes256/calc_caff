# Caffeine Population Pharmacokinetic Model - Perera et al.
# ------------------------------------------------------------------------------
# Parameters imputed assuming a bioavailability of 60% for inhaled caffeine.
#
# Model sourced from:
# Zandvliet AS, Huitema ADR, De Jonge ME, et al. Population Pharmacokinetics of 
#   Caffeine and its Metabolites Theobromine, Paraxanthine and Theophylline 
#   after Inhalation in Combination with Diacetylmorphine. Basic Clin Pharmacol 
#   Toxicol. 2005;96(1):71-79.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load libraries
  # library(dplyr)
  # library(mrgsolve)

# Define model code
  code <- '
$INIT  // Initial Conditions for Compartments (mg)
  GUT = 0,  // Absorption Compartment
  CENT = 0,  // Central Compartment
  PERI = 0,  // Peripheral Compartment
  PARA = 0,  // Paraxanthine Compartment
  THPH = 0,  // Theophylline Compartment
  THBR = 0,  // Theobromine compartment

$SET  // Set Differential Equation Solver Options			
  atol = 1e-8, rtol = 1e-8
  maxsteps = 100000

$PARAM  // Population parameters
  TVKA = 3.55,  // absorption rate constant (h-1) {taken from Perera et al.}
  TVV1 = 76.2,  // caffeine central volume of distribution (L)
  TVCL = 8.37,  // caffeine total clearance (L/h)
  TVV2 = 214,  // caffeine peripheral volume of distribution (L)
  TVQ = 9.74,  // inter-compartmental clearance (L/h)
  TVFPX = 0.0288,  // paraxanthine formation rate constant (L-1)
  TVKPX = 0.308,  // paraxanthine elimination rate constant (h-1)
  TVFTP = 0.00342,  // theophylline formation rate constant (L-1)
  TVKTP = 0.313,  // theophylline elimination rate constant (h-1)
  TVFTB = 0.00537,  // theobromine formation rate constant (L-1)
  TVKTB = 0.117,  // theobromine elimination rate constant (h-1)
  TVFCA = 0.6,  // approximate bioavailability (dimensionless)

  // Default Covariate Values for Simulation

  // Default ETA Values for Simulation
  // Allocated in population so set to zero
  ETA1 = 0,  // PPVKA
  ETA2 = 0,  // PPVV1
  ETA3 = 0,  // PPVCL
  ETA4 = 0,  // PPVFPX
  ETA5 = 0,  // PPVFTP
  ETA6 = 0,  // PPVFTB

$OMEGA  // Population parameter Variability
  name = "omega1"
  block = FALSE
  0.313600  // PPVKA

$OMEGA  // Population parameter Variability
  name = "omega2"
  block = TRUE
  0.269361  // PPVV1
  0.122723 0.316969  // PPVCL

$OMEGA  // Population parameter Variability
  name = "omega3"
  block = TRUE
  0.429025  // PPVFPX
  0.312147 0.414736  // PPVFTP
  0.446710 0.662805 1.537600  // PPVFTB

$SIGMA  // Residual Unexplained Variability	
  block = FALSE
  0.068644  // RUVEXPCA
  0.045369  // RUVEXPPX
  0.049284  // RUVEXPTP
  0.051076  // RUVEXPTB

$MAIN  // Individual Parameter Values
  double KA = TVKA*exp(ETA1);  // *exp(ETA(1))
  double V1 = TVFCA*TVV1*exp(ETA2);  // *exp(ETA(2))
  double CL = TVFCA*TVCL*exp(ETA3);  // *exp(ETA(3))
  double V2 = TVFCA*TVV2;
  double Q = TVFCA*TVQ;
  double FPX = TVFPX*exp(ETA4);  // *exp(ETA(4))
  double KPX = TVKPX;
  double FTP = TVFTP*exp(ETA5);  // *exp(ETA(5))
  double KTP = TVKTP;
  double FTB = TVFTB*exp(ETA6);  // *exp(ETA(6))
  double KTB = TVKTB;

$ODE  // Differential Equations
  double C1 = CENT/V1;
  double C2 = PERI/V2;
  double C3 = PARA/V1;
  double C4 = THPH/V1;
  double C5 = THBR/V1;

  dxdt_GUT = -KA*GUT;
  dxdt_CENT = KA*GUT - C1*Q + C2*Q - C1*CL;
  dxdt_PERI = C1*Q - C2*Q;
  dxdt_PARA = CENT*CL*FPX - PARA*KPX;
  dxdt_THPH = CENT*CL*FTP - THPH*KTP;
  dxdt_THBR = CENT*CL*FTB - THBR*KTB;

$TABLE  // Determines Values and Includes in Output	
  double IPRECA = C1;  // caffeine individual prediction
  double IPREPX = C3;  // paraxanthine individual prediction
  double IPRETP = C4;  // theophylline individual prediction
  double IPRETB = C5;  // theobromine individual prediction

  double DVCA = IPRECA*exp(EPS(1));  // caffeine observed conc
  double DVPX = IPREPX*exp(EPS(2));  // paraxanthine observed conc
  double DVTP = IPRETP*exp(EPS(3));  // theophylline observed conc
  double DVTB = IPRETB*exp(EPS(4));  // theobromine observed conc

$CAPTURE 
  IPRECA IPREPX IPRETP IPRETB DVCA DVPX DVTP DVTB 
  GUT CENT PERI PARA THPH THBR KA V1 CL V2 Q FPX KPX FTP KTP FTB KTB
  ETA1 ETA2 ETA3 ETA4 ETA5 ETA6 EPS(1) EPS(2) EPS(3) EPS(4)
'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Compile the model code
  mod <- mrgsolve::mcode("CaffPK", code)