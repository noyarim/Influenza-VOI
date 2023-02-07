%Purpose:
%Using RCGP data, compute cumulative cases (sum of weekly ILI rate per 100,000)
%for 2009/10 - 2017/2018 seasons, subdivided by class (A/B), subtype, lineage

%Values are scaled using weekly influenza positivity data

%Author: Ed Hill
%--------------------------------------------------------------------------

clear variables

%% APPLY STRAIN COMPOSITION DATA TO ILI DATA (Undetermined Flu B split evenly)

%Load strain composition data
%Variable name: InfAB_Mix_IgnoreUndetSamples 
load('../WHOFluNet/StrainMixResultsUK/UK_InfAB_Mix_IgnoreUndetSamples.mat')

%Load ILI data
%Variable name: NonAgeModel_SeasonRateTotalSum_FluPositiveScaled 
load('NonAgeModel_InfluenzaPositiveGPConsultRateAllStrains_2009to2018.mat')
%%
%Amend InfAB_Mix_IgnoreUndetSamples so in 2009/10 season, Inf B is split
%evenly between both lineages. Similar for 2011/12 season
InfAB_Mix_IgnoreUndetSamples_Amended = InfAB_Mix_IgnoreUndetSamples(:,1:4);
InfAB_Mix_IgnoreUndetSamples_Amended(1,3:4) = InfAB_Mix_IgnoreUndetSamples(1,5)/2;
InfAB_Mix_IgnoreUndetSamples_Amended(3,3:4) = InfAB_Mix_IgnoreUndetSamples(3,5)/2;
%%
%Multiply strain proportions by ILI data from same season
% StrainSeasonRateTotalSum: Array 
%   -> rows - Season 2009/10, 2010/11, ... , 2015/16. 2016/17 
%   -> cols - A/H1N1, A/H3, B/Yam, B/Vic, B/All undet.  
StrainSeasonRateTotalSum = NonAgeModel_SeasonRateTotalSum_FluPositiveScaled.*InfAB_Mix_IgnoreUndetSamples_Amended;

%Update 2009/2010 weights (row 1) to be solely assigned to A/H1N1
StrainSeasonRateTotalSum(1,:) = 0;
StrainSeasonRateTotalSum(1,1) = NonAgeModel_SeasonRateTotalSum_FluPositiveScaled(1);

%% Save ILI rate, strain stratified data to file
save('InfluenzaPositiveGPConsultRateByStrain_2009to2018.mat','StrainSeasonRateTotalSum')
