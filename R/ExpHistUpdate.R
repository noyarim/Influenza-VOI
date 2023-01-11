ExpHistUpdate <- function(ExpHistVaccType,ExpHistArray,ExpHistArrayParams,VaccEfficacy,ExpHistNum,NumOfStrains){
  
  #Inputs
  # ExpHistVaccType - Flag variable to specify how susceptibility will be modified for
  #                   vaccine-related exposure history classes
  #                  vacc. groups is reduced (0 - Not reduced, 1 - Reduced)
  # ExpHistArray -
  # ExpHistArrayParams - Contains initial ExpHistArray values, with
  #                       additional paramter inputs as extra cell entries
  # VaccEfficacy - For last season, effectiveness of vaccine
  # ExpHistNum - Total number of exposure histories
  # NumOfStrains
  
  #Outputs
  # ExpHistArray - Updated version of exposure history array
  if (ExpHistVaccType == 1){
    
    #Values to access particular array indexes
    HalfExpHistNum = as.integer(ExpHistNum/2)
    
    #Assign previous season efficacy scaling value to new variable name
    ExpHistVaccScaling = ExpHistArrayParams[3]
    
    #Update entry in ExpHistArray corresponding to vaccianted in previous season
    #ExpHistArray(:,HalfExpHistNum+1:end) = ones(NumOfStrains,HalfExpHistNum) - repmat(ExpHistVaccScaling*VaccEfficacy,1,HalfExpHistNum);
    for (ii in c(1:HalfExpHistNum)){
      ExpHistArray[,(HalfExpHistNum+ii)] = rep(1,NumOfStrains) - (ExpHistVaccScaling*VaccEfficacy)
    }
    
    for (ii in c(2:HalfExpHistNum)){
      #Update entries: natural infection
      ExpHistArray[(ii-1),(ii+HalfExpHistNum)] = min(ExpHistArrayParams[1],ExpHistArray[(ii-1),(ii+HalfExpHistNum)])
    }
    
    #Update entries: influenza B cross-reactivity
    ExpHistArray[(NumOfStrains-1),ncol(ExpHistArray)] = min(ExpHistArrayParams[2],ExpHistArray[(NumOfStrains-1),ncol(ExpHistArray)])
    ExpHistArray[NumOfStrains,(ncol(ExpHistArray)-1)] = min(ExpHistArrayParams[2],ExpHistArray[NumOfStrains,(ncol(ExpHistArray)-1)])
    
  }
  return (ExpHistArray)
  
}