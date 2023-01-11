library(Matrix)
library(stats)
library(lhs)
#Purpose:
#Contains functions to perform Adaptive Population Monte Carlo Approximate Bayesian Computation
#Run in serial

#Fit multi-strain, non-age structured influenza transmission model
#(with immunity propagation) to empirical data

#Code Author: Ed Hill, replicated in R by Kyu Lee
#-------------------------------------------------------------------------------

### GROUPING OF FUNCTIONS TO BE PASSED INTO THE APMC SCHEME ###
# (i) FUNCTION TO PERTURB SAMPLES WHEN GENERATING PROPOSED PARAMETER SETS
# (ii) FUNCTIONS TO DRAW SAMPLES FROM PRIOR
# (iii) SPECIFY PRIOR DENSITY
# (iv) SUMMARY STATISTIC CALCULATION FUNCTIONS
# (v) FUNCTION TO RUN MODEL SIMULATION & PRODUCE DESIRED OUTPUTS TO FEED INTO SUMMARY STATISTIC FUNCTION
# (vi) FUNCTIONS USED IN MODEL SIMULATION

#-------------------------------------------------------------------------------
# (i) FUNCTION TO PERTURB SAMPLES WHEN GENERATING PROPOSED PARAMETER SETS
#-------------------------------------------------------------------------------

OLCMPerturbVarFn <- function(RetainedParams,RetainedWeights,SurvivorParticleIdx,GenT){
  # Optimal local covariance matrix (OLCM)
  # Uses a multivariate normal distribution with a covariance matrix based on a subset of the particles from the
  # previous iteration, whose distances are smaller than the threshold of the current iteration
  
  #Inputs:
  #   RetainedParams,RetainedWeights - Current set of samples with associated weights
  # 	SurvivorParticleIdx - Boolean true/false states defining the subset of the particles from the previous iteration, whose distances are smaller than the threshold of the current iteratio
  #	GenT - Generation number, used to determine if should sample using a global covariance
  #Outputs:
  #   C - Variance-covariance matrix. Three dimensional, slice per particle
  
  #Define number of parameters and particles
  #RetainedParamsDims = length(RetainedParams)
  ParamNum = ncol(RetainedParams)  #Number of parameters matches number of columns in parameter array
  ParticleNum = nrow(RetainedParams) #Number of particles in generation matches number of rows in parameter array
  
  #Update covariance arrays
  if (GenT == 1){ #First generation,
    C_SingleSlice = 2.0*cov.wt(RetainedParams,wt=RetainedWeights,method="unbiased")$cov
      #cov(RetainedParams,AnalyticWeights(RetainedWeights[:]),corrected=true)
  
    #Check if C is non-singular. Modify if singular.
    while (rankMatrix(C_SingleSlice)[1] != ParamNum){ #Check
      C_SingleSlice = C_SingleSlice + 1e-12*diag(ncol(C_SingleSlice)) #C_SingleSlice + 1e-12I
    }
    
    #Repeat covariance matrix, once per parameter set
    ParticleNum = nrow(RetainedParams) #Number of particles in generation matches number of rows in parameter array
    C = replicate(ParticleNum,C_SingleSlice)
  }else{
    #Get parameter sets from the previous iteration, whose distances are smaller than the threshold of the current iteration
    # i.e. "those that survived", row ID of RetainedParams array dennoted by SurvivorParticleIdx
    SPP_thetas = RetainedParams[as.vector(SurvivorParticleIdx),] 	#SPP, subset of previous particles
  
    #Pick out weights for subset of the particles from the
    # previous iteration, whose distances are smaller than the threshold of the current iteration
    # Normalise those weights
    SPP_Weights = RetainedWeights[as.vector(SurvivorParticleIdx)]
    SPP_NormWeights = SPP_Weights/sum(SPP_Weights)
    
    #Compute mean of the survivor particle population
    m = colSums(SPP_thetas*SPP_NormWeights)
    
    #Initialise variance-covariance array
    C = array(0,c(ParamNum,ParamNum,ParticleNum))
    
    #Loop through each particle and compute covariance
    for (kk in c(1:ParticleNum)){
      Current_C = matrix(0,ParamNum, ParamNum)
      for (jj in c(1:ParamNum)){
        for (ii in c(jj:ParamNum)){
          Current_C[ii, jj] = sum(SPP_NormWeights*(SPP_thetas[,ii] - m[ii])*(SPP_thetas[, jj] - m[jj])) + (m[ii] - RetainedParams[kk,ii])*(m[jj] - RetainedParams[kk,jj])
    
          #Covariance matrix is symmetric.
          #Assign transposed off-diagonal value
          if (ii != jj){
          Current_C[jj, ii] = Current_C[ii, jj]}
        }
      }
    
    
      #Check if Current_C is non-singular. Blow it up if singular.
      while (rankMatrix(Current_C) != ParamNum){ #Check
        Current_C = Current_C + 1e-12*diag(ncol(Current_C)) #Current_C + 1e-12I
      }
    
    #Assign revised covariance to C (collection of covariance arrays)
      C[,,kk] = Current_C
    }
  
  return (C)
  }
}

#-------------------------------------------------------------------------------
# (ii) FUNCTION TO DRAW SAMPLES FROM PRIOR
#-------------------------------------------------------------------------------

SampleFirstGenFn_FromFile <- function(N){
  path = "/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Results/ParamInference/APMCOutputFiles/"
#Inputs:
#   N - Number of samples from LHC parameter set to take forward as first generation
  ParticlesFromPrior = read.table(paste0(path,"APMCsamples_JulRun#$(RunID)_RetainedParams.txt"))
  
  #Assign weights to particle
  Particle_Weights = read.table(paste0(path,"APMCsamples_JulRun#$(RunID)_RetainedWeights.txt"))
  
  #Get SurvivorParticleIdx
  SurvivorParticleIdx = read.table(paste0(path,"APMCsamples_JulRun#$(RunID)_SurvivorParticleIdx.txt"))
  
  #Get summary statistic measure for each paramter set
  Particle_SummStat = read.table(paste0(path,"APMCsamples_JulRun#$(RunID)_RetainedSummStat.txt"))

  return(list(ParticlesFromPrior, Particle_Weights, SurvivorParticleIdx, Particle_SummStat))
}

#FITTING TO 2012/2013 - 2015/2016 INFLUENZA SEASONS
#Draw samples from specified prior ranges, with ascertainment prob. per
#season
SampleFirstGenFn_FourSeasonFit <- function(N){
  #Inputs:
  #   N - Number of samples from prior to take forward as first generation
  
  #Distributions to be sampled from for each variable
  d = matrix(c(0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, #Transmissibility
       0,1, 0,1, 0,0, #Immunity propagation params
       0.01,0.01 ,0.01,0.01 ,0.01,0.01 ,0.01,0.01), ncol=2, byrow = TRUE) #Strain modifier for susceptibility
  
  ParamToFitNum = nrow(d) #Total number of parameters to be fitted
  
  #Sample and assign to array
  LHSsample = randomLHS(N, ParamToFitNum)
  ParticlesFromPrior = matrix(0,N,ParamToFitNum)
  for (ii in c(1:ParamToFitNum)){
    ParticlesFromPrior[,ii] = qunif(LHSsample[,ii], d[ii,1], d[ii,2])
  }
  
  #Assign weights to particle
  Particle_Weights = rep(1,N)
  
  return (list(ParticlesFromPrior, Particle_Weights))
}

#FITTING TO 2012/2013 - 2016/2017 INFLUENZA SEASONS
#Draw samples from specified prior ranges, with ascertainment prob. per
#season
SampleFirstGenFn_FiveSeasonFit <-function(N){
  #Inputs:
  #   N - Number of samples from prior to take forward as first generation
  
  #Distributions to be sampled from for each variable
  d = matrix(c(0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, #Transmissibility
               0,1, 0,1, 0,0, #Immunity propagation params
               0.9,1 ,0.9,1 ,0.9,1 ,0.9,1, 0.9,1), ncol=2, byrow = TRUE) #Strain modifier for susceptibility

  ParamToFitNum = nrow(d) #Total number of parameters to be fitted
  
  #Sample and assign to array
  ParticlesFromPrior = matrix(0,N,ParamToFitNum)
  for (ii in c(1:ParamToFitNum)){
    ParticlesFromPrior[,ii] = runif(N, d[ii,1], d[ii,2])
  }
  
  #Assign weights to particle
  Particle_Weights = rep(1,N)
  
  return (list(ParticlesFromPrior, Particle_Weights))
}

#FITTING TO 2012/2013 - 2017/2018 INFLUENZA SEASONS
#Draw samples from specified prior ranges, with ascertainment prob. per
#season
SampleFirstGenFn_SixSeasonFit <- function(N){
#Inputs:
#   N - Number of samples from prior to take forward as first generation

#Distributions to be sampled from for each variable
  d = matrix(c(0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, #Transmissibility
               0,1, 0,1, 0,0, #Immunity propagation params
               0.01,0.01 ,0.01,0.01 ,0.01,0.01 ,0.01,0.01, 0.01,0.01, 0.01,0.01), ncol=2, byrow = TRUE) #Strain modifier for susceptibility

  ParamToFitNum = nrow(d) #Total number of parameters to be fitted
  
  #Sample and assign to array
  ParticlesFromPrior = matrix(0,N,ParamToFitNum)
  for (ii in c(1:ParamToFitNum)){
    ParticlesFromPrior[,ii] = runif(N, d[ii,1], d[ii,2])
  }
  
  #Assign weights to particle
  Particle_Weights = rep(1,N)
  
  return (list(ParticlesFromPrior, Particle_Weights))
  }


#FITTING TO 2012/2013 - 2017/2018 INFLUENZA SEASONS
#Draw samples from specified prior ranges, with ascertainment prob. per
#season
#Includes parameter for carry-over immunity resulting from natural infection
SampleFirstGenFn_SixSeasonFitPlusMultiSeasonImm <-function(N){
  #Inputs:
  #   N - Number of samples from prior to take forward as first generation

  #Distributions to be sampled from for each variable
  d = matrix(c(0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, #Transmissibility
               0,1, 0,1, 0,0, #Immunity propagation params
               0.01,0.01 ,0.01,0.01 ,0.01,0.01 ,0.01,0.01, 0.01,0.01, 0.01,0.01), ncol=2, byrow = TRUE) #Strain modifier for susceptibility
  
  ParamToFitNum = nrow(d) #Total number of parameters to be fitted
  
  #Sample and assign to array
  ParticlesFromPrior = matrix(0,N,ParamToFitNum)
  for (ii in c(1:ParamToFitNum)){
    ParticlesFromPrior[,ii] = rand(d[ii],N)
  }
  
  #Assign weights to particle
  Particle_Weights = rep(1,N)
  
  return (list(ParticlesFromPrior, Particle_Weights))
  }

#FITTING TO 2012/2013 - 2019/2020 INFLUENZA SEASONS
#Draw samples from specified prior ranges, with ascertainment prob. per
#season
SampleFirstGenFn_EightSeasonFit <-function(N){
  #Inputs:
  #   N - Number of samples from prior to take forward as first generation
  
  #Distributions to be sampled from for each variable
  d = matrix(c(0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, #Transmissibility
               0.5,1, 0.5,1, 0,0, #Immunity propagation params
               0,0.01, 0,0.01, 0,0.01 ,0,0.01, 0,0.01, 0,0.01, 0,0.01, 0,0.01), ncol=2, byrow = TRUE) #Strain modifier for susceptibility

  ParamToFitNum = nrow(d) #Total number of parameters to be fitted
  
  #Sample and assign to array
  ParticlesFromPrior = matrix(0,N,ParamToFitNum)
  for (ii in c(1:ParamToFitNum)){
    ParticlesFromPrior[,ii] = runif(N, d[ii,1], d[ii,2])
  }
  
  #Assign weights to particle
  Particle_Weights = rep(1,N)
  
  return (list(ParticlesFromPrior, Particle_Weights))
}

#----------------------------------------------------------------------
# (iii) SPECIFY PRIOR DENSITY
#----------------------------------------------------------------------

#Per season ascertainment probability, fit to 2012/2013 to 2015/2016
#seasons inclusive
Prior_FourSeasonFit <- function(x){

  #Inputs:
  #   x - proposed particle value to be tested
  
  #Outputs:
  #   PriorProb - Prior likelihood value of current particle (parameter set)
  #   InBoundFlag - Denotes whether proposed particle values are within prior bounds
  
  #Get number of dimensions of x
  Particle_nDims = dim(x) #ndims(x)
  
  #--------------------------------------------------------------------------
  #COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
  #use)
  #Distributions each varaible were sampled from
  d = matrix(c(0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, #Transmissibility
               0,1, 0,1, 0,0, #Immunity propagation params
               0.01,0.01 ,0.01,0.01 ,0.01,0.01 ,0.01,0.01), ncol=2, byrow = TRUE) #Strain modifier for susceptibility
  
  if (is.null(Particle_nDims)){ #(Particle_nDims == 1){
    ParticleNum = 1
    ParamNum = length(x)
    
    PriorProb = dunif(x[1], d[1,1], d[1,2])
    
    for (ii in c(2:ParamNum)){
      PriorProb = PriorProb*dunif(x[ii], d[ii,1], d[ii,2]) 
    }
    
    #Alter flag value for implausible parameter set to 0
    if (PriorProb==0){
      InBoundFlag = 0
    }else{
      InBoundFlag = 1
    }
  }else{
    ParticleNum = dim(x)[1]
    ParamNum = dim(x)[2]
  
    PriorProb = dunif(x[,1], d[1,1], d[1,2])
    for (ii in c(2:ParamNum)){
      PriorProb = PriorProb*dunif(x[,ii], d[ii,1], d[ii,2])
    }
  
    #Alter flag value for implausible parameter set to 0
    InBoundFlag = rep(1,ParticleNum)
    InBoundFlag[PriorProb==0] = 0
  }
  
  return (list(PriorProb,InBoundFlag))
  
}

#Per season ascertainment probability, fit to 2012/2013 to 2016/2017
#seasons inclusive
Prior_FiveSeasonFit <- function(x){

  #Inputs:
  #   x - proposed particle value to be tested
  
  #Outputs:
  #   PriorProb - Prior likelihood value of current particle (parameter set)
  #   InBoundFlag - Denotes whether proposed particle values are within prior bounds
  
  #Get number of dimensions of x
  Particle_nDims = dim(x)
  
  #--------------------------------------------------------------------------
  #COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
  #use)
  #Distributions each varaible were sampled from
  d = matrix(c(0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, #Transmissibility
             0,1, 0,1, 0,0, #Immunity propagation params
             0.9,1 ,0.9,1 ,0.9,1 ,0.9,1 ,0.9,1), ncol=2, byrow = TRUE) #Strain modifier for susceptibility

  
  if (is.null(Particle_nDims)){
    ParticleNum = 1
    ParamNum = length(x)
    
    PriorProb = dunif(x[1], d[1,1],d[1,2])
    
    for (ii in 2:ParamNum){
      PriorProb = PriorProb*dunif(x[ii], d[ii,1], d[ii,2]) 
    }
    
    #Alter flag value for implausible parameter set to 0
    if (PriorProb==0){
      InBoundFlag = 0
    }else{
      InBoundFlag = 1}
    
  }else{
    ParticleNum = dim(x)[1]
    ParamNum = dim(x)[2]
    
    PriorProb = dunif(x[,1],d[1,1],d[1,2])
    for (ii in 2:ParamNum){
      PriorProb = PriorProb*dunif(x[,ii],d[ii,1],d[ii,2])
    }
  
  #Alter flag value for implausible parameter set to 0
  InBoundFlag = rep(1,ParticleNum)
  InBoundFlag[PriorProb==0] = 0
  }
  
  return (list(PriorProb,InBoundFlag))
  
  }

#Per season ascertainment probability, fit to 2012/2013 to 2017/2018
#seasons inclusive
Prior_SixSeasonFit <-function(x){

  #Inputs:
  #   x - proposed particle value to be tested
  
  #Outputs:
  #   PriorProb - Prior likelihood value of current particle (parameter set)
  #   InBoundFlag - Denotes whether proposed particle values are within prior bounds
  
  #Get number of dimensions of x
  Particle_nDims = dim(x)
  
  #--------------------------------------------------------------------------
  #COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
  #use)
  #Distributions each varaible were sampled from
  d = matrix(c(0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, #Transmissibility
               0,1, 0,1, 0,0, #Immunity propagation params
               0.01,0.01 ,0.01,0.01 ,0.01,0.01 ,0.01,0.01 ,0.01,0.01, 0.01,0.01), ncol=2, byrow = TRUE) #Strain modifier for susceptibility
  
  
  if (is.null(Particle_nDims)){
    ParticleNum = 1
    ParamNum = length(x)
    
    PriorProb = dunif(x[1], d[1,1],d[1,2])
    
    for (ii in 2:ParamNum){
      PriorProb = PriorProb*runif(x[ii],d[ii,1],d[ii,2])
    }
    
    #Alter flag value for implausible parameter set to 0
    if (PriorProb==0){
      InBoundFlag = 0
    }else{
      InBoundFlag = 1
    }
  }else{
    ParticleNum = dim(x)[1]
    ParamNum = dim(x)[2]
  
    PriorProb = dunif(x[,1], d[1,1],d[1,2])
    for (ii in 2:ParamNum){
      PriorProb = PriorProb*dunif(x[,ii], d[ii,1],d[ii,2])
    }
  
  #Alter flag value for implausible parameter set to 0
  InBoundFlag = ones(Int64,ParticleNum)
  InBoundFlag[PriorProb==0] = 0
  }
  
  return (list(PriorProb,InBoundFlag))
  
}

#Per season ascertainment probability, fit to 2012/2013 to 2017/2018
#seasons inclusive
#Includes parameter for carry-over immunity resulting from natural infection
Prior_SixSeasonFitPlusMultiSeasonImm <-function(x){

  #Inputs:
  #   x - proposed particle value to be tested
  
  #Outputs:
  #   PriorProb - Prior likelihood value of current particle (parameter set)
  #   InBoundFlag - Denotes whether proposed particle values are within prior bounds
  
  #Get number of dimensions of x
  Particle_nDims = dim(x)
  
  #--------------------------------------------------------------------------
  #COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
  #use)
  #Distributions each varaible were sampled from
  d = matrix(c(0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, #Transmissibility
               0,1, 0,1, 0,0, #Immunity propagation params
               0.01,0.01 ,0.01,0.01 ,0.01,0.01 ,0.01,0.01 ,0.01,0.01, 0.01,0.01), ncol=2, byrow = TRUE) #Strain modifier for susceptibility
  
  if (is.null(Particle_nDims)){
    ParticleNum = 1
    ParamNum = length(x)
    
    PriorProb = dunif(x[1], d[1,1],d[1,2])
    
    for (ii in 2:ParamNum){
      PriorProb = PriorProb*dunif(x[ii],d[ii,1],d[ii,2])
    }
  
    #Alter flag value for implausible parameter set to 0
    if (PriorProb==0){
      InBoundFlag = 0
    }else{
      InBoundFlag = 1
    }
  }else{
    ParticleNum = size(x,1)
    ParamNum = size(x,2)
  
    PriorProb = dunif(x[,1],d[1,1],d[1,2])
    for (ii in 2:ParamNum){
      PriorProb = PriorProb*dunif(x[,ii], d[ii,1],d[ii,2])
    }
    
    #Alter flag value for implausible parameter set to 0
    InBoundFlag = rep(1,ParticleNum)
    InBoundFlag[PriorProb==0] = 0
  }
  
  return (list(PriorProb,InBoundFlag))

}

#Per season ascertainment probability, fit to 2012/2013 to 2019/2020
#seasons inclusive
Prior_EightSeasonFit <- function(x){
  
  #Inputs:
  #   x - proposed particle value to be tested
  
  #Outputs:
  #   PriorProb - Prior likelihood value of current particle (parameter set)
  #   InBoundFlag - Denotes whether proposed particle values are within prior bounds
  
  #Get number of dimensions of x
  Particle_nDims = dim(x)
  
  #--------------------------------------------------------------------------
  #COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
  #use)
  #Distributions each varaible were sampled from
#  d = matrix(c(0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, #Transmissibility
#               0,1, 0,1, 0,0, #Immunity propagation params
#               0,0.01, 0,0.01 ,0,0.01 ,0,0.01 ,0,0.01, 0,0.01, 0,0.01, 0,0.01), ncol=2, byrow = TRUE) #Strain modifier for susceptibility
  d = matrix(c(0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, 0.2632,3*0.2632, #Transmissibility
               0.5,1, 0.5,1, 0,0, #Immunity propagation params
               0,0.01, 0,0.01, 0,0.01 ,0,0.01, 0,0.01, 0,0.01, 0,0.01, 0,0.01), ncol=2, byrow = TRUE) #Strain modifier for susceptibility
  
  
  if (is.null(Particle_nDims)){
    ParticleNum = 1
    ParamNum = length(x)
    
    PriorProb = dunif(x[1], d[1,1],d[1,2])
    
    for (ii in 2:ParamNum){
      PriorProb = PriorProb*dunif(x[ii], d[ii,1], d[ii,2]) 
    }
    
    #Alter flag value for implausible parameter set to 0
    if (PriorProb==0){
      InBoundFlag = 0
    }else{
      InBoundFlag = 1}
    
  }else{
    ParticleNum = dim(x)[1]
    ParamNum = dim(x)[2]
    
    PriorProb = dunif(x[,1],d[1,1],d[1,2])
    for (ii in 2:ParamNum){
      PriorProb = PriorProb*dunif(x[,ii],d[ii,1],d[ii,2])
    }
    
    #Alter flag value for implausible parameter set to 0
    InBoundFlag = rep(1,ParticleNum)
    InBoundFlag[PriorProb==0] = 0
  }
  
  return (list(PriorProb,InBoundFlag))
  
}

#-------------------------------------------------------------------------------
# (iv) SPECIFY SUMMARY STATISTIC CALCULATION
#-------------------------------------------------------------------------------
SummStatFun <- function(ObservedData,SimnData){

  #Disaggregate simulation data
  x = SimnData[[1]]
  Total_I = SimnData[[2]]
  
  #Compute poisson deviance
  SummStatIterNum_row = nrow(ObservedData)
  SummStatIterNum_col = ncol(ObservedData)
  TempSum = 0.0
  for (ii in 1:SummStatIterNum_row){
    for (jj in 1:SummStatIterNum_col){
    if (ObservedData[ii,jj] != 0){
      TempSum = TempSum + (ObservedData[ii,jj]*log(ObservedData[ii,jj]/x[ii,jj])) - (ObservedData[ii,jj] - x[ii,jj]) # calculate Poisson devinace
      }
    }
  }
  SummStatVal_Overall = 2*(TempSum + sum(x[ObservedData==0]))
  
  #Perform temporal check
  NumOfSeasons = as.integer(length(Total_I)/366) #%Number of seasons obtained by dividing number of daily records by days in yr (+1 to account for day 0 recording1)
  Total_I_Array = matrix(0,NumOfSeasons,366)
  for (jj in 1:NumOfSeasons){
    StartIdx = ((jj-1)*366) + 1
    EndIdx = jj*366
  
    Total_I_Array[jj,] = Total_I[StartIdx:EndIdx]
  }
  #Find day of seasonal year in which infection is at peak
  MaxInfValBySeason = apply(Total_I_Array[4:nrow(Total_I_Array),],1,max) #MaxInfValBySeason,inds = findmax(Total_I_Array[4:nrow(Total_I_Array),],2)
  MaxInfIdxBySeason = apply(Total_I_Array[4:nrow(Total_I_Array),],1,function(x) which(x == max(x))) #MaxInfIdxBySeason = map(x->ind2sub(Total_I_Array[4:nrow(Total_I_Array),], x)[2], inds)
  
  AmendedSummStatVal = SummStatVal_Overall
  
  #If peak outside Sep-Feb, set TemporalFlag to 0 and put amended
  #SummStatVal as infinity
   if (sum(MaxInfIdxBySeason>182) == 0){ #181 days September-February. Plus account for initial value
     peak_all = 0
   }else{
     peak_all = 1
   }
   
  # transform x (model outcome) to one-dimensional vector
  x_vector = rep(0,ncol(x)*nrow(x))
  for (kk in 1:nrow(x)){
    x_vector[(kk*4-3):(kk*4)] = as.matrix(x[kk,])
  }
  return (list(AmendedSummStatVal,x_vector,MaxInfIdxBySeason,peak_all))
  
}

#-------------------------------------------------------------------------------
# (v) FUNCTION TO RUN MODEL SIMULATION & PRODUCE DESIRED OUTPUTS TO FEED INTO
# SUMMARY STATISTIC FUNCTION
#-------------------------------------------------------------------------------
RunModelSimn <- function(FixedModelParams,x){
  #Inputs:
  #   FixedModelParams - Parameters with consistent values across all runs
  #   x - Values of parameters that are being inferred
  
  #Outputs:
  #   SimnData - Model output to be fed into summary statistic function
  
  #Disaggregate FixedModelParams inputs
  SimnRunType = FixedModelParams[[1]]
  ExpHistVaccType = FixedModelParams[[2]]
  StoreFlag_PopnFOI = FixedModelParams[[3]]
  SimnParam = FixedModelParams[[4]]
  NumOfStrains = FixedModelParams[[5]]
  ExpHistNum = FixedModelParams[[6]]
  InfectionParam = FixedModelParams[[7]]
  MultiSeasonImmPropn = FixedModelParams[[8]]
  BirthParam = FixedModelParams[[9]]
  DeathParam = FixedModelParams[[10]]
  VaccUptakeBySeason = FixedModelParams[[11]]
  LeakyTransFlag = FixedModelParams[[12]]
  LeakyVaccVarBySeason = FixedModelParams[[13]]
  InfPropn_StartOfSeason = FixedModelParams[[14]]
  ICFromFile = FixedModelParams[[15]]
  RetainedSeasonNum = FixedModelParams[[16]]
  
  #Update beta!
  InfectionParam[,1] = x[1:4]
  
  #Update ascertainment prob
  
  #Update exposure history
  #Build exposure history array. Assign to variable
  ExpHistArrayParams = x[5:7]
  ExpHistArray = BuildExpHistArray(NumOfStrains,ExpHistArrayParams)
  ExpHistArrayFnInputs = list(ExpHistArray,ExpHistArrayParams)
  
  #Run the model!
  StoreArrays = RunSeasonalFluModel(SimnParam,
                                    InfectionParam,BirthParam,DeathParam,VaccUptakeBySeason,LeakyVaccVarBySeason,LeakyTransFlag,
                                    ExpHistArrayFnInputs,ExpHistNum,ExpHistVaccType,MultiSeasonImmPropn,
                                    InfPropn_StartOfSeason,ICFromFile,SimnRunType)
  
  #Disaggregate StoreArrays
  T = StoreArrays[[1]]#::Array{Float64,1}
  C = StoreArrays[[2]]#::Array{Float64,2}
  S_NotV = StoreArrays[[3]]
  S_V = StoreArrays[[4]]
  E_NotV = StoreArrays[[5]]#::Array{Float64,2}
  E_V = StoreArrays[[6]]#::Array{Float64,2}
  I_NotV = StoreArrays[[7]]#::Array{Float64,2}
  I_V = StoreArrays[[8]]#::Array{Float64,2}
  
  #Visits in week (or month) t, denoted c_{t}, difference in the cumulative proportion
  #of consultations over the previous timestep
  #i.e. c_{t} = p(C(t)-C(t-1)),
  CumulCaseCount = rbind(C[seq(1,nrow(C),by=366),],C[nrow(C),]) #Get cumul. case count at specified intervals
  StrainSeasonRateTotalSum = CumulCaseCount[c(2:nrow(CumulCaseCount)),]-CumulCaseCount[c(1:(nrow(CumulCaseCount)-1)),]
  StrainSeasonRateTotalSum_FitCheckSeasons = StrainSeasonRateTotalSum[c(4:nrow(StrainSeasonRateTotalSum)),] #Get model simulated values for 2012/2013 onward
  
  #----------------------------------------------------------------------
  ### Get ascertainable cases, based on ascertainment prob type
  #----------------------------------------------------------------------
  AscertainProb = x[8:length(x)] #Parameter values from LHC file
  
  #Forward simulation modifications
  #Save counts unmodified by ascertainment probabilities
  if (SimnRunType == 3){
    StrainSeasonRatePerStrain = StrainSeasonRateTotalSum_FitCheckSeasons
  } else{
    StrainSeasonRatePerStrain = matrix(0,RetainedSeasonNum,4) #Initialise storage array
    
    #Update ascertainable amount per strain
    FluToGPconversion = 100000*AscertainProb #Scale to rate per 100,000 popualation
    for (jj in c(1:length(AscertainProb))){
      StrainSeasonRatePerStrain[jj,] = StrainSeasonRateTotalSum_FitCheckSeasons[jj,]*FluToGPconversion[jj]
    }
  }
  
  #Compute infected temporal profile
  Total_I = rowSums(I_NotV) + rowSums(I_V) + rowSums(E_NotV) + rowSums(E_V)
  Total_I_strain = I_NotV + I_V + E_NotV + E_V
  
  #immune history distribution in the beginning of 2021/2022 season
  immunehist = S_NotV[seq(1,nrow(S_NotV),by=366),]
  
  #Assign outputs to be used in Summary Statisitc function to array
  SimnData = list(StrainSeasonRatePerStrain, Total_I, Total_I_strain, immunehist)
  
  return (SimnData)
}

#Includes fitting parameter for carry-over immunity resulting from natural infection
RunModelSimnWithMultiSeasonImm <- function(FixedModelParams,x){
  #Inputs:
  #   FixedModelParams - Parameters with consistent values across all runs
  #   x - Values of parameters that are being inferred
  
  #Outputs:
  #   SimnData - Model output to be fed into summary statistic function
  
  #Disaggregate FixedModelParams inputs
  SimnRunType = FixedModelParams[[1]]
  ExpHistVaccType = FixedModelParams[[2]]
  StoreFlag_PopnFOI = FixedModelParams[[3]]
  SimnParam = FixedModelParams[[4]]
  NumOfStrains = FixedModelParams[[5]]
  ExpHistNum = FixedModelParams[[6]]
  InfectionParam = FixedModelParams[[7]]
  MultiSeasonImmPropn = FixedModelParams[[8]]
  BirthParam = FixedModelParams[[9]]
  DeathParam = FixedModelParams[[10]]
  VaccUptakeBySeason = FixedModelParams[[11]]
  LeakyTransFlag = FixedModelParams[[12]]
  LeakyVaccVarBySeason = FixedModelParams[[13]]
  InfPropn_StartOfSeason = FixedModelParams[[14]]
  ICFromFile = FixedModelParams[[15]]
  RetainedSeasonNum = FixedModelParams[[16]]
  
  #Update beta!
  InfectionParam[,1] = x[1:4]
  
  #Update ascertainment prob
  
  #Update exposure history
  #Build exposure history array. Assign to variable
  ExpHistArrayParams = x[5:7]
  ExpHistArray = BuildExpHistArray(NumOfStrains,ExpHistArrayParams)
  ExpHistArrayFnInputs = list(ExpHistArray,ExpHistArrayParams)
  
  #MultiSeasonImmPropn = x[8]#::Float64
  
  #Run the model!
  StoreArrays = RunSeasonalFluModel(SimnParam,
                                    InfectionParam,BirthParam,DeathParam,VaccUptakeBySeason,LeakyVaccVarBySeason,LeakyTransFlag,
                                    ExpHistArrayFnInputs,ExpHistNum,ExpHistVaccType,MultiSeasonImmPropn,
                                    InfPropn_StartOfSeason,ICFromFile,SimnRunType)
  
  
  #Disaggregate StoreArrays
  T = StoreArrays[[1]]#::Array{Float64,1}
  C = StoreArrays[[2]]#::Array{Float64,2}
  E_NotV = StoreArrays[[5]]#::Array{Float64,2}
  E_V = StoreArrays[[6]]#::Array{Float64,2}
  I_NotV = StoreArrays[[7]]#::Array{Float64,2}
  I_V = StoreArrays[[8]]#::Array{Float64,2}
  
  S_NotV = StoreArrays[[3]]
  S_V = StoreArrays[[4]]
  R_NotV = StoreArrays[[9]]
  R_V = StoreArrays[[10]]
  #Visits in week (or month) t, denoted c_{t}, difference in the cumulative proportion
  #of consultations over the previous timestep
  #i.e. c_{t} = p(C(t)-C(t-1)),
  CumulCaseCount = rbind(C[seq(1,nrow(C),by=366),],C[nrow(C),]) #Get cumul. case count at specified intervals
  StrainSeasonRateTotalSum = CumulCaseCount[c(2:nrow(CumulCaseCount)),]-CumulCaseCount[c(1:(nrow(CumulCaseCount)-1)),]
  StrainSeasonRateTotalSum_FitCheckSeasons = StrainSeasonRateTotalSum[c(4:nrow(StrainSeasonRateTotalSum)),] #Get model simulated values for 2012/2013 onward
  
  #----------------------------------------------------------------------
  ### Get ascertainable cases, based on ascertainment prob type
  #----------------------------------------------------------------------
  AscertainProb = x[9:ncol(x)] #Parameter values from LHC file
  
  #Forward simulation modifications
  #Save counts unmodified by ascertainment probabilities
  if (SimnRunType == 3){
    StrainSeasonRatePerStrain = StrainSeasonRateTotalSum_FitCheckSeasons
  } else {
    StrainSeasonRatePerStrain = matrix(0,RetainedSeasonNum,4) #Initialise storage array
    
    #Update ascertainable amount per strain
    FluToGPconversion = 100000*AscertainProb #Scale to rate per 100,000 popualation
    for (jj in c(1:nrow(AscertainProb))){
      StrainSeasonRatePerStrain[jj,] = StrainSeasonRateTotalSum_FitCheckSeasons[jj,]*FluToGPconversion[jj]
    }
  }
  
  #Compute infected temporal profile
  Total_I = rowSums(I_NotV) + rowSums(I_V) + rowSums(E_NotV) + rowSums(E_V)
  Store_Iall = I_NotV + I_V + E_NotV + E_V
  Total_V = rowSums(S_V) + rowSums(E_V) + rowSums(I_V) + rowSums(R_V)
  
  #Assign outputs to be used in Summary Statisitc function to array
  SimnData = list(StrainSeasonRatePerStrain,Total_I)
  
  return (SimnData)
}

#----------------------------------------------------------------------
# (vi) FUNCTIONS USED IN SUMMARY STATISTIC CALCULATION
#----------------------------------------------------------------------

BuildExpHistArray <- function(NumOfStrains,ExpHistArrayParams){
  #PURPOSE: Construct exposure history array based on current parameter set
  #Inputs:
  #   NumOfStrains - As titled
  #   ExpHistArrayParams - File containing parameter sets to be run
  #Outputs:
  #   ExpHistArray - interaction array between exposure history and susceptibility to
  #   the current season strain variant
  
  #Number of exposure history classes
  #One for naive, one for vacc with no natural infection, one per strain
  #(unvacc), one per strain with vacc
  ExpHistNum = (NumOfStrains*2) + 2
  
  #Column per exposure history. Row per strain.
  ExpHistArray = matrix(1,NumOfStrains,ExpHistNum)
  HalfExpHistNum = as.integer(ExpHistNum/2)
  
  #Column 1 - no exposure in previous season
  #Columns 2-5 - natural infection by one of the strains
  #Column 6  - vacc. in previous season, no natural infection
  #Columns 7-10 - Natural infection by ones of the sratins AND vaccinated
  
  #Right hand columns, vaccinated previous season
  ExpHistArray[,((HalfExpHistNum+1):ncol(ExpHistArray))] = ExpHistArrayParams[3]
  
  for (ii in c(2:HalfExpHistNum)){
    #Update entries: natural infection
    ExpHistArray[ii-1,ii] = ExpHistArrayParams[1]
    ExpHistArray[ii-1,(ii+HalfExpHistNum)] = ExpHistArrayParams[1]
  }
  
  #Update entries: influenza B cross-reactivity
  ExpHistArray[(NumOfStrains-1),c(HalfExpHistNum, ExpHistNum)] = ExpHistArrayParams[2]
  ExpHistArray[NumOfStrains,c(HalfExpHistNum-1, ExpHistNum-1)] = ExpHistArrayParams[2]
  
  return (ExpHistArray)
}


RunInferenceAPMC <- function(RunID,ObvsData,SeasonsToSimulate,
                          N_alpha,alpha,MinAcceptRate,MaxGen,PerturbVarScale,
                          PriorFn,SummStatFn,SampleFromPriorFn,ModelSimnFn,
                          FirstGenFromFileFlag){

  #-------------------------------------------------------------------------------
  # SET UP APMC SCHEME RELATED PARAMETERS
  #-------------------------------------------------------------------------------
  
  #Calculate number of samples before retention phase
  ScalingFactor = 1/alpha
  N = as.integer(round(N_alpha*ScalingFactor))
  
  #-------------------------------------------------------------------------------
  # DEFINE AND GROUP MODEL SIMULATION FIXED PARAMETERS
  #-------------------------------------------------------------------------------
  
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
  life_exp = 81.0
  b = 1/(life_exp*365)
  d = 1/(life_exp*365)
  BirthParam = rep(0,ExpHistNum)
  BirthParam[1] = b
  DeathParam = d
  
  #------------------------------------------------------------------------------
  ###  VACCINATION UPTAKE
  #--------------------------------------------------------------------------
  
  #BASED ON HISTORICAL DATA
  path = "/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/SeasonalFluImmunityPropagation-master/"
  ## Import the data - Pandemic flu vacc (2009/2010 season)
  PandemicFluVaccUptake = read.xlsx(paste0(path, "Data/VaccUptake/NonAgeStrucModel_DailyVaccUptakeCalYr_PandemicFluVacc_EMH.xlsx"),sheet="2009_2010",colNames = F, rows =5, cols = c(3:368))
  
  ## Import the data - Sesonal flu vacc
  HistoricalSeasonalFluVaccUptake = read.xlsx(paste0(path,"Data/VaccUptake/NonAgeStrucModel_DailyVaccUptakeBySeasonCalYr_EMH.xlsx"),sheet="All popn.", colNames = F, rows = c(5:12), cols = c(3:368))
  
  #Combine uptake vectors into a single array
  VaccUptakeBySeason_DataFrame = rbind(PandemicFluVaccUptake, HistoricalSeasonalFluVaccUptake)
  
  # Collate into Array
  VaccUptakeBySeason = as.matrix(VaccUptakeBySeason_DataFrame) #Array{Float64, 2}(VaccUptakeBySeason_DataFrame)
  
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
  
  #BASED ON HISTORICAL DATA
  VaccEfficacyDataFrame = read.xlsx(paste0(path,"Data/VaccEfficacy/VaccEfficacy_AllPopn.xlsx"),sheet="MidPoint", colNames = F, rows = c(3:11), cols = c(3:6))
  
  # Collate into Array
  VaccEfficacy = as.matrix(VaccEfficacyDataFrame) #Array{Float64, 2}(VaccEfficacyDataFrame)
  
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
  InfPropn_StartOfSeason = matrix(c(1e-5, 0, 0, 0, 2.5*10^-6, 2.5*10^-6, 2.5*10^-6, 2.5*10^-6), c(2,3), byrow=TRUE)
  
  
  #--------------------------------------------------------------------------
  ### GIVE DETAILS ON LOADING INITIAL SUSCEPTIBLE CLASS CONDITIONS FROM FILE IF NEEDED
  #--------------------------------------------------------------------------
  ICFromFile = c(0, NA)
  
  #--------------------------------------------------------------------------
  ### NUMBER OF INFLUENZA SEASONS WORTH OF DATA BEING CONSIDERED
  #--------------------------------------------------------------------------
  RetainedSeasonNum = nrow(ObvsData)
  
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
  OutputFileName = "/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Results/ParamInference/APMCOutputFiles/EndOfGenAPMC_"
  
  #--------------------------------------------------------------------------
  # RUN APMC SCHEME
  #--------------------------------------------------------------------------
  APMC_outcome <- APMC_loop(ObvsData,FirstGenFromFileFlag,
  	N,alpha,N_alpha,MinAcceptRate,MaxGen,PerturbVarScale,
      PriorFn,SummStatFn,SampleFromPriorFn,FixedModelParams,ModelSimnFn,RunID,OutputFileName)
  RetainedParams <- APMC_outcome[[1]]
  RetainedWeights <- APMC_outcome[[2]]
  RetainedSummStat <- APMC_outcome[[3]]
  SurvivorParticleIdx <- APMC_outcome[[4]]
  GenT <- APMC_outcome[[5]]
  #-------------------------------------------------------------------------------
  # SAVE TO FILE
  #-------------------------------------------------------------------------------
  
  #Specify filename to write sample values
  FName_GenT = "/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Results/ParamInference/APMCOutputFiles/APMCsamples_JulRun#$(RunID)_GenT.txt"
  FName_RetainedParams = "/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Results/ParamInference/APMCOutputFiles/APMCsamples_JulRun#$(RunID)_RetainedParams.txt"
  FName_RetainedSummStat = "/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Results/ParamInference/APMCOutputFiles/APMCsamples_JulRun#$(RunID)_RetainedSummStat.txt"
  FName_RetainedWeights = "/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Results/ParamInference/APMCOutputFiles/APMCsamples_JulRun#$(RunID)_RetainedWeights.txt"
  FName_SurvivorParticleIdx = "/Users/kyulee/Google Drive/Z Drive/Postdoc_UPitt/FRED/CDC_Influenza/Results/ParamInference/APMCOutputFiles/APMCsamples_JulRun#$(RunID)_SurvivorParticleIdx.txt"
  
  #--------------------------------------------------------------------------
  ### SAVE SUMMARY STATISTICS TO FILE
  #--------------------------------------------------------------------------
  write.table(GenT,FName_GenT)
  write.table(RetainedParams,FName_RetainedParams)
  write.table(RetainedSummStat,FName_RetainedSummStat)
  write.table(RetainedWeights,FName_RetainedWeights)
  write.table(SurvivorParticleIdx,FName_SurvivorParticleIdx)
  
  }
