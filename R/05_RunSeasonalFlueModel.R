#Purpose:
#Script to run Adaptive Population Monte Carlo Approximate Bayesian Computation
#Run serially & take arguments from command line as function inputs

#Fit multi-strain, non-age structured influenza transmission model
#(with immunity propagation) to empirical data

#Code Author: Ed Hill
#-------------------------------------------------------------------------------
#--------------------------------------------------------------------------
# ADD FILES TO SEARCH PATH FOR ODES/MODEL RUN FUNCTION
#--------------------------------------------------------------------------

setwd("/Users/kyueunlee/Library/CloudStorage/GoogleDrive-leex5499@umn.edu/My Drive/Z Drive/Project_all/2022_FluVOI/Influenza-VOI")
source("R/01_RunSeasonalFluModelODEs.R")
source("R/02_ExpHistUpdate.R")
source("R/03_SeasonalFluModel_Fns.R")


#-------------------------------------------------------------------------------
# LOAD REQUIRED PACKAGES
#-------------------------------------------------------------------------------
library(Matrix)
library(stats)
library(deSolve)
library(openxlsx)
library(ggplot2)
library(reshape2)
library(matlab)
library(MASS)
library(mvtnorm)


outdir = "output/"

ARGS = c("1600", # RunID
         "data/PSA_data/ILIData/EmpData_InfluenzaPositiveGPConsultRateByStrain_2009to2018.csv", # Hospitalization data
         "11", # Number of influenza seasons to simulate from 2009/2010 onwards (e.g. 7 will be up to and inc. 2015/16)
         "RunModelSimn" #APMC related functions
         )
#--------------------------------------------------------------------------
# Set RunID idx
#--------------------------------------------------------------------------
RunID = as.integer(ARGS[1])

#--------------------------------------------------------------------------
# Specify empirical data being fitted to
#--------------------------------------------------------------------------
ObvsDataFileName = ARGS[2]
ObvsDataTemp = read.table(ObvsDataFileName,sep=',')

#Specify number of influenza seasonsto be simulated
#From 2009/2010 onwards: - > 7 will be up to and inc. 2015/16
#                        - > 8 up to and inc. 2016/17 etc
SeasonsToSimulate = as.integer(ARGS[3])

#Get maximum number of seasons that could be used in inference procedure
#In case all seasons are not used, get number of seasons omitted
MaxSeasonNumToSimulate = nrow(ObvsDataTemp)
SeasonsUnused =  MaxSeasonNumToSimulate - SeasonsToSimulate

#Error check
IntCheck = (SeasonsToSimulate%%1)==0
if (IntCheck == 0 || SeasonsToSimulate <= 0 || (SeasonsToSimulate > MaxSeasonNumToSimulate)){
  stop("Incompatible SeasonsToFit value")
}

#Pick out row 4 onwards of ObvsDataTemp.
#Row 4 corresponds to 2012/13 influenza season, row 5 to 2013/14 influenza season etc
ObvsData = ObvsDataTemp[(4:(nrow(ObvsDataTemp)-SeasonsUnused)),]

#--------------------------------------------------------------------------
# SET UP APMC SCHEME RELATED PARAMETERS
#--------------------------------------------------------------------------
s_5 = ARGS[4] #Convert string to Symbol
#Make Symbols callable functions
ModelSimnFn = get(s_5) #Model simulation


#-------------------------------------------------------------------------------
# SPECIFY TYPE OF RUN THROUGH FLAG VARIABLE
#-------------------------------------------------------------------------------

#TYPE OF RUN
# (INFLUENCES VACCINE UPTAKE/EFFICACY, & USE OF ODEBurnIn, ODEH1N1OnlyTime, ODEAlleStrainTime)
#1 - exploratory; 2 - historical; 3 - alternative vacc. scheme
SimnRunType = 2 #Fixed to 2, as fitting to historical data

#Select exposure history susceptibility modification form for vaccine-related states
#0 - Fixed values every season
#1 - Relate to previous season vaccine efficacy
ExpHistVaccType = 1

if (ExpHistVaccType != 0 && ExpHistVaccType !=1){
  stop("ExpHistVaccType must take value 0 or 1")}


#--------------------------------------------------------------------------
### DISEASE DYNAMICS SIMN/FLAG PARAMETERS
#--------------------------------------------------------------------------

SimnStartDate=9 #Month of year to start simulation on

#Store population-level FOI flag option
#0 - inactive, 1 - active
StoreFlag_PopnFOI = 0

#Run time for ODE model. Take values post burn in
ODEBurnIn = 0*365 #No entity information recorded
ODEH1N1OnlyTime = 1*365 #Time spect recording entity values, but only H1N1 infection allowed
ODEAllStrainTime = (SeasonsToSimulate-1)*365 #Infection by all strains now possible

ODESampleTime = ODEH1N1OnlyTime + ODEAllStrainTime #Timeframe over which data outputs are stored
MaxTime = ODEBurnIn + ODESampleTime #Total simlation run time

#Time over which only H1N1 type infection allowed
TotalTime_H1N1Only = ODEBurnIn + ODEH1N1OnlyTime

#Get number of seasons
NumOfSeasons = as.integer(ODESampleTime/365)

#Specify timestep to use in ode45 scheme
timestep = 1

#Concatenate simulation variables
SimnParam=list(SimnStartDate,ODEBurnIn,ODESampleTime,MaxTime,timestep,TotalTime_H1N1Only,StoreFlag_PopnFOI)

#------------------------------------------------------------------------------
###  TRANSMISSION RELATED PARAMETERS
#------------------------------------------------------------------------------
NumOfStrains = 4 #Specify number of strains in the system (A/H1N1, A/H3N2, two B lineages)
beta = c(0.4390,0.5026,0.4263,0.4255) #Placeholder beta values. Will be replaced during inference process.
gamma = c(1/3.8*rep(1,NumOfStrains)) #recovery rate
sigma = c(1/1.4,1/1.4,1/0.6,1/0.6)  #latent rate
InfectionParam = matrix(c(beta,sigma,gamma), nrow = 4)

#Set proportion of population that are susceptible at end of season within
#a natural infection exposure history class to retain pre-existing immunity
MultiSeasonImmPropn = 0

#------------------------------------------------------------------------------
###  DEMOGRAHPY PARAMETERS
#------------------------------------------------------------------------------
#Number of exposure history classes
#One for naive, one for vacc with no natural infection, one per strain
#(unvacc), one per strain with vacc
ExpHistNum = (NumOfStrains*2) + 2

#Set per capita birth&death rates
life_exp = 78.0
b = 1/(life_exp*365)
d = 1/(life_exp*365)
BirthParam = rep(0,ExpHistNum)
BirthParam[1] = b
DeathParam = d


#VaxEffFactor_list = c(0.5,0.75,1,1.25,1.5)
#VaxUpFactor_list = c(0.5,0.75,1,1.25,1.5)
#for (i in VaxEffFactor_list){
  VaxEffFactor = 1
#  for (j in VaxUpFactor_list){
    VaxUpFactor = 1
#    print(i,j)
    #------------------------------------------------------------------------------
    ###  VACCINATION UPTAKE
    #--------------------------------------------------------------------------
    
    #BASED ON HISTORICAL DATA
    path = ""
    ## Import the data - Pandemic flu vacc (2009/2010 season)
    PandemicFluVaccUptake = read.xlsx(paste0(path, "data/PSA_data/VaccUptake/NonAgeStrucModel_DailyVaccUptakeCalYr_PandemicFluVacc_EMH.xlsx"),sheet="2009_2010",colNames = F, rows =5, cols = c(3:368))
      
    ## Import the data - Sesonal flu vacc
    HistoricalSeasonalFluVaccUptake = read.xlsx(paste0(path,"data/PSA_data/VaccUptake/NonAgeStrucModel_DailyVaccUptakeBySeasonCalYr_EMH.xlsx"),sheet="All popn.", colNames = F, rows = c(5:14), cols = c(3:368))
      
    #Combine uptake vectors into a single array
    VaccUptakeBySeason_DataFrame = rbind(PandemicFluVaccUptake, HistoricalSeasonalFluVaccUptake)
      
    # Collate into Array
    VaccUptakeBySeasonUpTo2020 = as.matrix(VaccUptakeBySeason_DataFrame)
      
    #Extend to cover up to and including 2049/2050 influenza season (assuming that vaccine uptake is the same as the most recent uptake)
    VaccUptakeBySeason = matrix(0,SeasonsToSimulate,365)
    VaccUptakeBySeason[c(1:nrow(VaccUptakeBySeasonUpTo2020)),] = VaccUptakeBySeasonUpTo2020
#    LastVaccUptake = VaccUptakeBySeasonUpTo2020[nrow(VaccUptakeBySeasonUpTo2020),] # last observed vacc uptake in 2020
#    VaccUptakeBySeason[c((nrow(VaccUptakeBySeasonUpTo2020)+1):SeasonsToSimulate),] = matrix(rep(LastVaccUptake,SeasonsToSimulate-nrow(VaccUptakeBySeasonUpTo2020)), ncol=365,byrow=TRUE) #nrow=nrow(LastVaccUptake)*(SeasonsToSimulate-9)) #repmat(VaccUptakeBySeasonUpTo2018[end,:]',SeasonsToSimulate-9,1)
#    VaccUptakeBySeason[c((nrow(VaccUptakeBySeasonUpTo2020)+2):SeasonsToSimulate),] = VaxUpFactor*matrix(rep(LastVaccUptake,SeasonsToSimulate-nrow(VaccUptakeBySeasonUpTo2020)-1), ncol=365,byrow=TRUE) #nrow=nrow(LastVaccUptake)*(SeasonsToSimulate-9)) #repmat(VaccUptakeBySeasonUpTo2018[end,:]',SeasonsToSimulate-9,1)
    
    
    
    #--------------------------------------------------------------------------
    ### VACCINATION - LEAKY TRANSMISSION SETTINGS
    #--------------------------------------------------------------------------
    #Set flag variable for "leaky" transmission being unactive or active
    #0 - Infectiousness of infected vacc. group unmodified.
    #1 - Infected vacc. group has reduced infectiousness
    LeakyTransFlag = 0
    
    #Validity check on LeakyTransFlag value
    if (LeakyTransFlag !=0 && LeakyTransFlag !=1){
      stop("Invalid value of $LeakyTransFlag for LeakyTransFlag, should be 0 or 1")
    }
    
    #--------------------------------------------------------------------------
    ### VACCINATION - EFFICACY
    #--------------------------------------------------------------------------
    #Set leaky vaccine efficacy parameters
    #Row i - Strain i (A/H1N1, A/H3N2, two B lineages)
    #Cell 1 - Susceptibility reduction (as propn)
    #Cell 2 - infectiousness reduction (as propn)
    
    
    #RUN ALTERNATIVE VACC. SCHEMES
    VaccEfficacyDataFrame = read.xlsx(paste0(path,"data/PSA_data/VaccEfficacy/VaccEfficacy_AllPopn.xlsx"),sheet="MidPoint", colNames = F, rows = c(3:13), cols = c(3:6))
    
    # Collate into Array
    VaccEfficacyUpTo2020 = as.matrix(VaccEfficacyDataFrame)
    
    #Extend to cover up to and including 2029/2030 influenza season
    VaccEfficacy = matrix(0,SeasonsToSimulate,NumOfStrains)
    VaccEfficacy[c(1:nrow(VaccEfficacyUpTo2020)),] = VaccEfficacyUpTo2020
    
    #Assign to LeakyVaccVar variable
    if (LeakyTransFlag == 0){
      LeakyVaccVarBySeason = t(VaccEfficacy)
    } else if (LeakyTransFlag == 1){
      alpha = t(VaccEfficacy)
      delta = t(VaccEfficacy)
      
      LeakyVaccVarBySeason = array(c(alpha,delta), dim = c(nrow(alpha),ncol(alpha),2))
    }
    
    
    #--------------------------------------------------------------------------
    ### INITIAL INF. PROPORTION
    #--------------------------------------------------------------------------
    #Column i for strain i
    #Row 1 for when only H1N1 infection allowed
    #Row 2 when infection by any strain allowed
    InfPropn_StartOfSeason = matrix(c(1e-5, 0, 0, 0, 
                                      2.5*10^-6, 2.5*10^-6, 2.5*10^-6, 2.5*10^-6), c(2,4), byrow=TRUE)
    
    
    #--------------------------------------------------------------------------
    ### GIVE DETAILS ON LOADING INITIAL SUSCEPTIBLE CLASS CONDITIONS FROM FILE IF NEEDED
    #--------------------------------------------------------------------------
    ICFromFile = c(0, NA)
    
    #--------------------------------------------------------------------------
    ### ASSIGN TO VARIABLE IF THESE SET OF SIMULATIONS ARE FORWARD SIMULATIONS
    #--------------------------------------------------------------------------
    ForwardSimnFlag = 1
    VaccEffTestFlag = 2
    
    
    #--------------------------------------------------------------------------
    ### NUMBER OF INFLUENZA SEASONS WORTH OF DATA BEING CONSIDERED
    #--------------------------------------------------------------------------
    #Running forward simulation. All seasons will be saved to output
    RetainedSeasonNum = SeasonsToSimulate
    
    #--------------------------------------------------------------------------
    # AGGREGATE FIXED PARAMETERS
    #--------------------------------------------------------------------------
    FixedModelParams = list(SimnRunType,ExpHistVaccType,StoreFlag_PopnFOI,
                            SimnParam,NumOfStrains,ExpHistNum,InfectionParam,MultiSeasonImmPropn,BirthParam,
                            DeathParam,VaccUptakeBySeason,LeakyTransFlag,LeakyVaccVarBySeason,
                            InfPropn_StartOfSeason,ICFromFile,RetainedSeasonNum)
    
    #--------------------------------------------------------------------------
    # SET END-OF-GENERATION OUTPUT TEXT FILE INFO
    #--------------------------------------------------------------------------
    OutputFileName = "output/EndOfGenAPMC_"
    
    #ObservedData = ObvsData # TODO: delete this 
    #SampleFirstGenFn = SampleFromPriorFn#TODO: delete this
    #SampleFirstGenFn_output <- SampleFirstGenFn(N)
    #FullParamSet <- SampleFirstGenFn_output[[1]]
    #Particle_Weights <- SampleFirstGenFn_output[[2]]
      
    # Read the best parameter set from calibration 
    b_param_unscaled <- read.csv("data/b_param_unscaled_8seasons.csv")
    b_param_unscaled <- b_param_unscaled[,-c(1)]
    b_param_unscaled <- cbind(b_param_unscaled[1:6], 0, b_param_unscaled[1,7:length(b_param_unscaled)])
    b_param_unscaled = as.vector(t(b_param_unscaled))
    #b_param_unscaled[8:length(b_param_unscaled)] = 1 # set no assertainment
    
    
#    VaccEfficacy2010To2020 = VaccEfficacyUpTo2020[c(2:nrow(VaccEfficacyUpTo2020)),]
    #Get mean efficacy value per strain, from the the historic flu seasons 2010/2011-2017/2018 to extrapolate vax eff
 #   MeanEffVals = apply(VaccEfficacy2010To2020,2,mean) #findmax(VaccEfficacyUpTo2018[2:nrow(VaccEfficacyUpTo2018),],1)[1]
 #   VaccEfficacy[(nrow(VaccEfficacyUpTo2020)+1):nrow(VaccEfficacy),] = matrix(rep(MeanEffVals, SeasonsToSimulate - nrow(VaccEfficacyUpTo2020)),ncol=4, byrow = TRUE) 
 #   VaccEfficacy[(nrow(VaccEfficacyUpTo2020)+2):nrow(VaccEfficacy),] = VaxEffFactor*matrix(rep(MeanEffVals, SeasonsToSimulate - nrow(VaccEfficacyUpTo2020)-1),ncol=4, byrow = TRUE) 
    
    #Pass updated vaccine efficacy array into FixedModelParams
    FixedModelParams[[13]] = t(VaccEfficacy)
    
    SimnData = ModelSimnFn(FixedModelParams,b_param_unscaled)
    ratebystrain <- SimnData[[1]]
    tot_I <- SimnData[[2]]
    immunehist <- SimnData[[4]]
    write.csv(ratebystrain, paste0(outdir, "ratebystrain_Vaxeff",VaxEffFactor,"_VaxUp",VaxUpFactor,".csv"))
    write.csv(tot_I, paste0(outdir, "tot_I_Vaxeff",VaxEffFactor,"_VaxUp",VaxUpFactor,".csv"))
    write.csv(immunehist, paste0(outdir, "immunehist_Vaxeff",VaxEffFactor,"_VaxUp",VaxUpFactor,".csv"))
#  }
#}
