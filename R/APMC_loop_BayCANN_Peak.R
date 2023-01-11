library(matlab)
library(MASS)
library(mvtnorm)

APMC_loop <-function(ObservedData,FirstGenFromFileFlag,N,alpha,N_alpha,
                   MinAcceptRate,MaxGen,PerturbVarFn,
                   PriorFn,SummStatFn,SampleFirstGenFn,
                   FixedModelParams,ModelSimnFn,RunID,OutputFileName){

  #Inputs:
  #   ObservedData - (float array) The data!
  #   FirstGenFromFileFlag - (binary) Indicator for method of obtaining initial sample (1 from file with particles from previous APMC run, 0 to draw from prior)
  #   N - (Int scalar) Number of particles in generation pre retention phase
  #   alpha - (Float scalar) Proportion of samples kept in retention phase
  #   N_alpha -  (Int scalar) Number of retained samples
  #   MinAcceptRate - (Float scalar) Stopping criterion: Threshold for propn of newly genertated samples that
  #                   must attain superior summary statistic to threshold summary statistic
  #                   value epsilon for simulation to continue
  #   MaxGen - (Int scalar) Forcibly stop inference procedure if this amount of geneations
  #                   is reached
  #   PerturbVarFn - (fn handle) Used as covariance matrix in Gaussian distribution when perturbing samples
  #   PriorFun - (fn handle) Prior calculation & check if perturbed sample is in prior bounds
  #   SummStatFun - (fn handle) summary statistic calculation
  #   SampleFirstGenFn - (fn handle) specifying generation of first generation
  #                           of particles & associated weights
  #                      -> generated from prior distribution, equal weights
  #                      -> import samples/weights from file
  #   FixedModelParams - (Tuple) Batch of quantities to be passed into model simulation function
  #   ModelSimnFn - (fn handle) Perform model simulation
  #   RunID - (Int scalar) Identifier for APMC run
  #   OutputFileName - (string) Directory location for end-of-generation data to be
  #                     saved to
  
  #Outputs:
  # RetainedParams,RetainedWeights,RetainedSummStat - Result: N_alpha
  #           particles with assoicated weights and summary statistic values
  # SurvivorParticleIdx - Logical (true/false) vector, gives those N_alpha particles
  #                       NOT newly generated in final generation
  # GenT - Number of generations until stopping criterion satisfied.
  
  #-------------------------------------------------------------------------------
  #RUN INITIAL LOOP
  #-------------------------------------------------------------------------------
  if (FirstGenFromFileFlag==1){ #Get first generation particles from previous inference run
    #Initialisation
    SampleFirstGenFn_output <- SampleFirstGenFn(N_alpha)
    RetainedParams <- SampleFirstGenFn_output[[1]]
    Particle_Weights_Temp <- SampleFirstGenFn_output[[2]]
    SurvivorParticleIdxTemp <- SampleFirstGenFn_output[[3]]
    SummStat_Temp <- SampleFirstGenFn_output[[4]]
    #Convert SurvivorParticleIdx to BitArray
    SurvivorParticleIdx = SurvivorParticleIdxTemp #SurvivorParticleIdx = convert(BitArray,SurvivorParticleIdxTemp) #TODO: Bitarray?
    
    #Get number of parameters being inferred
    ParamInferredNum = ncol(RetainedParams)
    FullParamSet = matrix(0,N,ParamInferredNum) #Initialise full parameter set array
    
    #Assign weights to particle
    Particle_Weights = rep(1,N)
    Particle_Weights[1:N_alpha] = as.vector(Particle_Weights_Temp$x)
    
    #Initialise summary statistic storage vector
    SummStat = rep(0,N)
    SummStat[1:N_alpha] = as.vector(SummStat_Temp$x)
    
    #Return the first alpha-quartile of SummStat
    epsilon = SummStat[N_alpha]
    
    #Obtain first-generation of accepted particles
    RetainedWeights = Particle_Weights[1:N_alpha]
    RetainedSummStat = SummStat[1:N_alpha]
    
    #Set covariance matrix from accepted samples (based on information from file)
    GenNum = 1
    C = PerturbVarFn(RetainedParams,RetainedWeights,SurvivorParticleIdx,GenNum)
  }else {#Get first generation particles from prior
    #Initialisation
    ObservedData = ObvsData # TODO: delete this 
    SampleFirstGenFn = SampleFirstGenFn_FourSeasonFit#TODO: delete this
    SampleFirstGenFn_output <- SampleFirstGenFn(N)
    FullParamSet <- SampleFirstGenFn_output[[1]]
    Particle_Weights <- SampleFirstGenFn_output[[2]]
  
    #Initialise summary statistic storage vector
    SummStat = rep(0,N)
    SimnRateAll = matrix(0,N,ncol(ObservedData)*nrow(ObservedData))
    PeakDateAll = matrix(0,N,nrow(ObservedData))
    PeakAll = rep(0,N)
    # First tolerance, sample from prior
    for (ii in 1:N){
      print(paste0("Prior param:",ii))
      #Run model with designated parameter set
      #println("FullParamSet: $(FullParamSet[ii,:])")
      SimnData = ModelSimnFn(FixedModelParams,FullParamSet[ii,])
      
      #Generate summary statistic for current prior set
      SoS_out = SummStatFn(ObservedData,SimnData)
      SoS = SoS_out[[1]]
      SimnRate = SoS_out[[2]]
      PeakDate = SoS_out[[3]]
      PeakBi = SoS_out[[4]]
      SummStat[ii] = SoS
      SimnRateAll[ii,] = as.matrix(SimnRate)
      PeakDateAll[ii,] = as.matrix(PeakDate)
      PeakAll[ii] = PeakBi
      #println("SoS: $SoS")
      
#      if (is.na(SoS)){    #Convert implausible (NaN) outputs to inf
#        SummStat[ii] = Inf
#      }else{
#        SummStat[ii] = SoS
#      }
    }
    #Find particle idxs that are to be retained (below epsilon)
    # EpsilonIdx = order(SummStat)
    # ParamRetained = EpsilonIdx[1:N_alpha]
    # 
    # #Obtain first-generation of accepted particles
    # RetainedParams = FullParamSet[ParamRetained,]
    # RetainedWeights = Particle_Weights[ParamRetained]
    # RetainedSummStat = SummStat[ParamRetained]
    # 
    # #Return the first alpha-quartile of SummStat
    # epsilon = RetainedSummStat[length(RetainedSummStat)]
    # #Set covariance matrix from accepted samples (in first generation set)
    # GenNum = 1
    # SurvivorParticleIdx = 0 #Placeholder value. As first generation SurvivorParticleIdx not used within PerturbVarFn
    # C = PerturbVarFn(RetainedParams,RetainedWeights,SurvivorParticleIdx,GenNum)
  }
  
  # print("sample from prior is completed")
  # #Reorder FullParamSet, Weights, SummStat. Place retained values before unretained
  # FullParamSet[1:N_alpha,] = as.matrix(RetainedParams)
  # Particle_Weights[1:N_alpha] = RetainedWeights
  # SummStat[1:N_alpha] = RetainedSummStat
  # 
  # #Set acceptance rate to 1
  # AcceptRate = 1.0
  # 
  # #Set up generation index
  # GenT = 1
  #-------------------------------------------------------------------------------
  # Write param values, weights, etc to output files
  #-------------------------------------------------------------------------------
  ParamSetFileName = paste0(OutputFileName,paste0("_AllParamSets_Run#",paste0(RunID,".csv")))
#  WeightsFileName = paste0(OutputFileName,paste0("_Weights_Run#",paste0(RunID,".txt")))
  SummStatFileName =  paste0(OutputFileName,paste0("_AllSummStat_Run#",paste0(RunID,".csv")))
#  ThresholdValFileName = paste0(OutputFileName,paste0("_ThresholdVal_Run#",paste0(RunID,".txt")))
#  AccRateFileName = paste0(OutputFileName,paste0("_AccRate_Run#",paste0(RunID,".txt")))
#  SurvivorParticleIdxFileName = paste0(OutputFileName,paste0("_SurvivorParticleIdx_Run#",paste0(RunID,".txt")))
  SimnOutcomeFileName =  paste0(OutputFileName,paste0("_AllSimnOutcome_Run#",paste0(RunID,".csv")))
  SimnPeakDateFileName =  paste0(OutputFileName,paste0("_AllSimnPeakDate_Run#",paste0(RunID,".csv")))
  SimnPeakAllFileName =  paste0(OutputFileName,paste0("_AllSimnPeakAll_Run#",paste0(RunID,".csv")))
  
  write.csv(FullParamSet,ParamSetFileName, row.names=FALSE)
#  write.table(RetainedWeights,WeightsFileName, row.names=FALSE)
  write.csv(SummStat,SummStatFileName, row.names=FALSE)
  write.csv(SimnRateAll,SimnOutcomeFileName, row.names=FALSE)
  write.csv(PeakDateAll,SimnPeakDateFileName, row.names=FALSE)
  write.csv(PeakAll,SimnPeakAllFileName, row.names=FALSE)
  #  write.table(epsilon,ThresholdValFileName, row.names=FALSE)
#  write.table(AcceptRate,AccRateFileName, row.names=FALSE)
#  write.table(SurvivorParticleIdx,SurvivorParticleIdxFileName, row.names=FALSE)
  #-------------------------------------------------------------------------------
  print("writing prior is completed")
  #-------------------------------------------------------------------------------
  # RUN APMC LOOP
  #-------------------------------------------------------------------------------
  ParamNum = size(RetainedParams,2) #Number of parameters
  while ((AcceptRate > MinAcceptRate) && (GenT < MaxGen)){
    #tic()
    print(paste0("GenT is",GenT))
    #Update generation index
    GenT = GenT + 1
  
    #----------------------------------------------------------------------
    # Sample N-N_alpha particles from the previous weighted sample
    #----------------------------------------------------------------------
    #Generate random numbers to be used for drawing from previous parameter set
    #generation
    PickParticleRand = runif(N - N_alpha)
    
    #Construct vector of particle weights. Append 0 as first entry.
    CumSumWeights = cumsum(RetainedWeights)/sum(RetainedWeights)
  
    #Draw from previous parameter set generation
    #If using local covariance, store selected covariance (associated with sample chsoen to be perturbed)
    if (ndims(C) == 3){
      CovarSetsDrawn = rep(list(NA),N-N_alpha) #Array{Array{Float64,2}}(N-N_alpha)
    }
    ParamSetsDrawn = zeros(N-N_alpha,ParamNum)
    for (ii in 1:length(PickParticleRand)){
      ParamDrawIdx = min(which(CumSumWeights >= PickParticleRand[ii])) #searchsortedfirst(CumSumWeights, PickParticleRand[ii])
      ParamSetsDrawn[ii,] = as.matrix(RetainedParams[ParamDrawIdx,])
    
      if (ndims(C) == 3){
        CovarSetsDrawn[[ii]] = C[,,ParamDrawIdx]}
    }
  
  
    #----------------------------------------------------------------------
    
    #----------------------------------------------------------------------
    # Generate new particles with the proposal distribution
    #----------------------------------------------------------------------
    PerturbedParamSets = matrix(0,N-N_alpha,ParamNum)
    for (jj in 1:(N-N_alpha)){ #Resample out of range param. sets
    
      if (ndims(C) == 1){ #C is vector of component-wise variances
        d2_mean =ParamSetsDrawn[jj,]
        d2_sd = sqrt(C)
      } else if (ndims(C) == 2){ #C is a global covariance
        d2_mean=ParamSetsDrawn[jj,]
        d2_sd = C
      } else if (ndims(C) == 3){ #C is a local covariance
        d2_mean=as.vector(ParamSetsDrawn[jj,])
        d2_sd = as.matrix(CovarSetsDrawn[[jj]])
      }else{
        stop("Covariance array has a dimension above 3. Not compatible. Program ended.")
      }
      
      PerturbedParamSets[jj,] = mvrnorm(mu = d2_mean, Sigma = d2_sd)
    }
  
  
    #Check if perturbed parameter set is in prior bounds.
    PriorFnOutputTuple = PriorFn(PerturbedParamSets)
    InBoundFlag_NewParamSet = PriorFnOutputTuple[[2]]
    
    #Indicate if particle is contained within prior
    InBoundFlag_AllParamSet = c(rep(1,N_alpha),InBoundFlag_NewParamSet)
    #----------------------------------------------------------------------
  
    #Assign perturbed parameter sets to concatenated paramter set array
    #PerturbedParamSets
    FullParamSet[((N_alpha + 1):N),] = PerturbedParamSets
    for (ii in ((N_alpha + 1):N)){
      #tic()
      #Assign current particle/param set to variable
      CurrentParamSet = FullParamSet[ii,]
      
      #Assign validity of current particle to variable (does it lie in prior
      #dist?)
      InBoundFlag_CurrentParamSet = InBoundFlag_AllParamSet[ii]
      
      #If particle not in prior bounds, set weight to zero and SummStat value to
      #Inf. Otherwise, run model.
      if (InBoundFlag_CurrentParamSet == 0){
        SummStat[ii] = Inf
        Particle_Weights[ii] = 0
      }else{
    
        #Run model with designated parameter set
        SimnData = ModelSimnFn(FixedModelParams,CurrentParamSet)
      
        #Generate summary statistic for current prior set
        SummStatVal_NoPrior = SummStatFn(ObservedData,SimnData)
      
        #elapsed_time_1 = toc()
        #println("elapsed_time_1: $elapsed_time_1")
  
        if (is.na(SummStatVal_NoPrior) ){    #Convert implausible (NaN) outputs to inf
          SummStat[ii] = Inf
        }else{
          SummStat[ii] = SummStatVal_NoPrior
        }
  
        #-------------------------------------------------------------------
        #Assign weights to particle
        #-------------------------------------------------------------------
        
        #Assign weight to particle - initialise count variables
        WeightsDenom = 0.0
        RetainedWeightsSum = sum(RetainedWeights)
        NormWeights = RetainedWeights/RetainedWeightsSum
        
        #Assign weights to particle - variance variables
        if (ndims(C) == 1){ #C is vector of component-wise variances
          sigma = sqrt(C)
        
          #Assign weights to particle - evaluation loop
          for (jj in (1:N_alpha)){
            d1 = Normal(RetainedParams[jj,1], sigma[1]) #Declare distribution
            KernelVal = pdf(d1,CurrentParamSet[1]) #Get pdf for current parameter set
            for (kk in (2:ParamNum)){
              d1 = Normal(RetainedParams[jj,kk], sigma[kk])
              KernelVal = KernelVal*pdf(d1,CurrentParamSet[kk])
            }
  
            WeightsDenom = WeightsDenom + (NormWeights[jj]*KernelVal)
  
          }
        }else if (ndims(C) == 2){ # C is a global variance-covariance array
          #Assign weights to particle - evaluation loop
          for (jj in (1:N_alpha)){
            d1_mean <- RetainedParams[jj,]
            d1_sd <- C #Declare distribution
            KernelVal = dmvnorm(CurrentParamSet, mean=d1_mean, sigma=d1_sd) #Get pdf for current parameter set
          
            WeightsDenom = WeightsDenom + (NormWeights[jj]*KernelVal)
          }
  
        }else if (ndims(C) == 3){  # C is a local variance-covariance array
          #Assign weights to particle - evaluation loop
          for (jj in (1:N_alpha)){
            d1_mean <- as.matrix(t(RetainedParams[jj,]))
            d1_sd <- as.matrix(C[,,jj])             #d1 = MvNormal(RetainedParams[jj,], C[,,jj]) #Declare distribution
            KernelVal = dmvnorm(as.matrix(t(CurrentParamSet)), mean=d1_mean, sigma=d1_sd) #pdf(d1,CurrentParamSet) #Get pdf for current parameter set
          
            WeightsDenom = WeightsDenom + (NormWeights[jj]*KernelVal)
          }
  
        }else{
          stop("Covariance array has a dimension above 3. Not compatible. Program ended.")
        }
  
        #For particles lying outside prior range, amend weights to be zero
        PriorFn_output <- PriorFn(CurrentParamSet)
        PriorVal <- PriorFn_output[[1]]
        InBoundFlagCheck <- PriorFn_output[[2]]
        
        Particle_Weights[ii] = PriorVal/WeightsDenom
        #-------------------------------------------------------------------
      }
  
    }
  
    #Update acceptance rate
    AcceptRate = ((sum(SummStat[(N_alpha+1):length(SummStat)]<epsilon))/(N-N_alpha))
  
    #Find particle idxs that are to be retained
    EpsilonIdx = order(SummStat)
  
    ParamRetained = EpsilonIdx[1:N_alpha]
  
    #Obtain those particles that were present at beginnng of generation,
    #and have "survived" (been retained) at end of generation
    SurvivorParticleIdx = ParamRetained<=N_alpha
    
    #Obtain next-generation of accepted particles
    RetainedParams = FullParamSet[ParamRetained,]
    RetainedWeights = Particle_Weights[ParamRetained]
    RetainedSummStat = SummStat[ParamRetained]
    
    #Return the first alpha-quartile of SummStat
    epsilon = RetainedSummStat[length(RetainedSummStat)]
    
    #Set covariance matrix from accepted samples
    C = PerturbVarFn(RetainedParams,RetainedWeights,SurvivorParticleIdx,GenT)
    
    #Reorder FullParamSet, Weights, SummStat. Place retained values before unretained
    FullParamSet[1:N_alpha,] = RetainedParams
    Particle_Weights[1:N_alpha] = RetainedWeights
    SummStat[1:N_alpha] = RetainedSummStat
    
    print("GenT number $GenT complete")
  
    #-------------------------------------------------------------------------------
    # Write param values, weights, etc to output files
    #-------------------------------------------------------------------------------
    write.table(RetainedParams,ParamSetFileName,append=TRUE, row.names=FALSE, col.names=FALSE)
    write.table(t(RetainedWeights),WeightsFileName,append=TRUE, row.names=FALSE, col.names=FALSE)
    write.table(RetainedSummStat,SummStatFileName,append=TRUE, row.names=FALSE, col.names=FALSE)
    write.table(epsilon,ThresholdValFileName,append=TRUE, row.names=FALSE, col.names=FALSE)
    write.table(AcceptRate,AccRateFileName,append=TRUE, row.names=FALSE, col.names=FALSE)
    write.table(t(SurvivorParticleIdx),SurvivorParticleIdxFileName,append=TRUE, row.names=FALSE, col.names=FALSE)
    
#    io1 = open(ParamSetFileName, "a")
#    io2 = open(WeightsFileName, "a")
#    io3 = open(SummStatFileName, "a")
    #    io4 = open(ThresholdValFileName, "a")
    #io5 = open(AccRateFileName, "a")
    #io6 = open(SurvivorParticleIdxFileName, "a")
    
    #writedlm(io1,RetainedParams)
    #writedlm(io2,RetainedWeights)
    #writedlm(io3,RetainedSummStat)
    #writedlm(io4,epsilon)
    #writedlm(io5,AcceptRate)
    #writedlm(io6,SurvivorParticleIdx)
    
    #close(io1); close(io2); close(io3); close(io4); close(io5); close(io6);
    
    #-------------------------------------------------------------------------------
    
  }
  
  return (list(RetainedParams,RetainedWeights,RetainedSummStat,SurvivorParticleIdx,GenT))

}

