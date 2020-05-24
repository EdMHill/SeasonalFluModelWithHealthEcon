#Purpose:
#Script to call MultiParamSetEconEvalTest function
#Run health economic evaluation multiple times, with distinct
#epidemiological data and health economic parameter sets

#Set up to run in parallel as multiple jobs
#Each job runs altered vaccine scheme

#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# LOAD REQUIRED ENVIRONMENT (FILE PATH BASED ON CWD BEING THE LOCATION THIS FILE RESIDES)
#--------------------------------------------------------------------------
@everywhere using Pkg
@everywhere Pkg.activate("../")


#--------------------------------------------------------------------------
#Load required packages
#--------------------------------------------------------------------------
@everywhere using MAT
@everywhere using XLSX
@everywhere using Distributions
@everywhere using Random
@everywhere using DelimitedFiles

#--------------------------------------------------------------------------
#Load to path directories containing required function files
#--------------------------------------------------------------------------
@everywhere include("FatalQALYLossFnsJulV1.jl")
@everywhere include("DoHealthEconEvalJulV1.jl")
@everywhere include("HealthEconParamByAgeConstructionJulV1.jl")

#--------------------------------------------------------------------------
### Functions to carry out probabilistic sensitivity analysis
#--------------------------------------------------------------------------

###################################
### Vaccine schedule data      ####
###################################

@everywhere function VaccEoSUptake_ExpandSchoolAndElderlyNoUnder2s_SpanAgeSpace(JobID,AgeGrpNum,RiskGroupFlag)
#Inputs:
#   JobID - (integer) Used to determine how VaccUptake array should be populated
#   AgeGrpNum - (integer) Number of single year age-brackets in use
#   RiskGroupFlag - (Flag variable) 0: Low risk; 1: At risk.
# Outputs:
#   VaccUptake - (tuple) Uptake propn at end of season per age. Tuple entry per influenza season

    #State relevant data file names
    AtRiskOnlyVaccUptake_SpreadsheetName = "../../Data/VaccUptake/WholePopnVaccUptakeCalc/AlternativeVaccSchedules/VaccUptakePerSeason_AtRiskSchedule_EMH.xlsx"
    ExpandedProgVaccUptake_SpreadsheetName = "../../Data/VaccUptake/WholePopnVaccUptakeCalc/AlternativeVaccSchedules/VaccUptakePerSeason_AllPopnHighVaccUptake_EMH.xlsx"

    #Declare sheetnames within files of interest
    SheetNames = ["2012_2013","2013_2014","2014_2015","2015_2016","2016_2017","2017_2018"]

    #Specify relevant ranges
    if RiskGroupFlag == 0
        FieldRange = "D2:D102"
    elseif RiskGroupFlag == 1
        FieldRange = "C2:C102"
    end

    #Assign rows of uptake array based on JobIdx (read from text file)
    VaccAgeBoundsData = readdlm("VaccUptakeSweepFiles/VaccUptakeAgeSweep_BoundsData_NoUnder2sVers.txt",',')
    YoungClusterUB = VaccAgeBoundsData[JobID,1] #First entry the upper bound for young age cluster
    ElderClusterLB = VaccAgeBoundsData[JobID,2]   #Second entry for lower bound of elder age cluster

    println("JobID: $JobID")
    println("YoungClusterUB: $YoungClusterUB")
    println("ElderClusterLB: $ElderClusterLB")

    #Initialise storage array. Column per influenza season. Row per single year age bracket.
    VaccUptake = zeros(AgeGrpNum,length(SheetNames))

    #Iterate through each entry
    for ii = 1:length(SheetNames)
        #Import vaccine uptake data
        AtRiskOnlyProg_VaccUptakePercentages = XLSX.readdata(AtRiskOnlyVaccUptake_SpreadsheetName,"$(SheetNames[ii])",FieldRange)
        ExpandedProg_VaccUptakePercentages = XLSX.readdata(ExpandedProgVaccUptake_SpreadsheetName,"$(SheetNames[ii])",FieldRange)

        # Collate input data into Arrays & convert from percentages to proportions
        AtRiskOnlyProg_VaccUptakePropns = Array{Float64, 2}(AtRiskOnlyProg_VaccUptakePercentages/100)
        ExpandedProg_VaccUptakePropns = Array{Float64, 2}(ExpandedProg_VaccUptakePercentages/100)

        #Populate vaccine uptake array for current season iteration
        VaccUptakeBySeasonArray = AtRiskOnlyProg_VaccUptakePropns #Initial data. Match at risk only programme

        if isnan(YoungClusterUB)  #Elder cluster only
            ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
            VaccUptakeBySeasonArray[(ElderClusterLB+1):end] = ExpandedProg_VaccUptakePropns[(ElderClusterLB+1):end] #Amend elder age cluster rows
        elseif isnan(ElderClusterLB) #Young ages only
            YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
            VaccUptakeBySeasonArray[3:(YoungClusterUB+1)] = ExpandedProg_VaccUptakePropns[3:(YoungClusterUB+1)] #Amend young age cluster rows
        else #Both ages
            YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
            ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
            VaccUptakeBySeasonArray[3:(YoungClusterUB+1)] = ExpandedProg_VaccUptakePropns[3:(YoungClusterUB+1)] #Amend young age cluster rows
            VaccUptakeBySeasonArray[(ElderClusterLB+1):end] = ExpandedProg_VaccUptakePropns[(ElderClusterLB+1):end] #Amend elder age cluster rows
        end

        #Assign to output array
        VaccUptake[:,ii] = VaccUptakeBySeasonArray
    end

    #println("VaccUptake: $(VaccUptake[51,:])")

return VaccUptake
end

@everywhere function VaccEoSUptake_ExpandSchoolAndElderlyNoUnder4s_SpanAgeSpace(JobID,AgeGrpNum,RiskGroupFlag)
#Inputs:
#   JobID - (integer) Used to determine how VaccUptake array should be populated
#   AgeGrpNum - (integer) Number of single year age-brackets in use
#   RiskGroupFlag - (Flag variable) 0: Low risk; 1: At risk.
# Outputs:
#   VaccUptake - (tuple) Uptake propn at end of season per age. Tuple entry per influenza season

    #State relevant data file names
    AtRiskOnlyVaccUptake_SpreadsheetName = "../../Data/VaccUptake/WholePopnVaccUptakeCalc/AlternativeVaccSchedules/VaccUptakePerSeason_AtRiskSchedule_EMH.xlsx"
    ExpandedProgVaccUptake_SpreadsheetName = "../../Data/VaccUptake/WholePopnVaccUptakeCalc/AlternativeVaccSchedules/VaccUptakePerSeason_AllPopnHighVaccUptake_EMH.xlsx"

    #Declare sheetnames within files of interest
    SheetNames = ["2012_2013","2013_2014","2014_2015","2015_2016","2016_2017","2017_2018"]

    #Specify relevant ranges
    if RiskGroupFlag == 0
        FieldRange = "D2:D102"
    elseif RiskGroupFlag == 1
        FieldRange = "C2:C102"
    end

    #Assign rows of uptake array based on JobIdx (read from text file)
    VaccAgeBoundsData = readdlm("VaccUptakeSweepFiles/VaccUptakeAgeSweep_BoundsData_NoUnder4sVers.txt",',')
    YoungClusterUB = VaccAgeBoundsData[JobID,1] #First entry the upper bound for young age cluster
    ElderClusterLB = VaccAgeBoundsData[JobID,2]   #Second entry for lower bound of elder age cluster

    println("JobID: $JobID")
    println("YoungClusterUB: $YoungClusterUB")
    println("ElderClusterLB: $ElderClusterLB")

    #Initialise storage array. Column per influenza season. Row per single year age bracket.
    VaccUptake = zeros(AgeGrpNum,length(SheetNames))

    #Iterate through each entry
    for ii = 1:length(SheetNames)
        #Import vaccine uptake data
        AtRiskOnlyProg_VaccUptakePercentages = XLSX.readdata(AtRiskOnlyVaccUptake_SpreadsheetName,"$(SheetNames[ii])",FieldRange)
        ExpandedProg_VaccUptakePercentages = XLSX.readdata(ExpandedProgVaccUptake_SpreadsheetName,"$(SheetNames[ii])",FieldRange)

        # Collate input data into Arrays & convert from percentages to proportions
        AtRiskOnlyProg_VaccUptakePropns = Array{Float64, 2}(AtRiskOnlyProg_VaccUptakePercentages/100)
        ExpandedProg_VaccUptakePropns = Array{Float64, 2}(ExpandedProg_VaccUptakePercentages/100)

        #Populate vaccine uptake array for current season iteration
        VaccUptakeBySeasonArray = AtRiskOnlyProg_VaccUptakePropns #Initial data. Match at risk only programme

        if isnan(YoungClusterUB)  #Elder cluster only
            ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
            VaccUptakeBySeasonArray[(ElderClusterLB+1):end] = ExpandedProg_VaccUptakePropns[(ElderClusterLB+1):end] #Amend elder age cluster rows
        elseif isnan(ElderClusterLB) #Young ages only
            YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
            VaccUptakeBySeasonArray[5:(YoungClusterUB+1)] = ExpandedProg_VaccUptakePropns[5:(YoungClusterUB+1)] #Amend young age cluster rows
        else #Both ages
            YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
            ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
            VaccUptakeBySeasonArray[5:(YoungClusterUB+1)] = ExpandedProg_VaccUptakePropns[5:(YoungClusterUB+1)] #Amend young age cluster rows
            VaccUptakeBySeasonArray[(ElderClusterLB+1):end] = ExpandedProg_VaccUptakePropns[(ElderClusterLB+1):end] #Amend elder age cluster rows
        end

        #Assign to output array
        VaccUptake[:,ii] = VaccUptakeBySeasonArray
    end

    #println("VaccUptake: $(VaccUptake[51,:])")

return VaccUptake
end

@everywhere function VaccEoSUptake_ExpandSchoolAndElderly_SpanAgeSpace(JobID,AgeGrpNum,RiskGroupFlag)
#Inputs:
#   JobID - (integer) Used to determine how VaccUptake array should be populated
#   AgeGrpNum - (integer) Number of single year age-brackets in use
#   RiskGroupFlag - (Flag variable) 0: Low risk; 1: At risk.
# Outputs:
#   VaccUptake - (tuple) Uptake propn at end of season per age. Tuple entry per influenza season

    #State relevant data file names
    AtRiskOnlyVaccUptake_SpreadsheetName = "../../Data/VaccUptake/WholePopnVaccUptakeCalc/AlternativeVaccSchedules/VaccUptakePerSeason_AtRiskSchedule_EMH.xlsx"
    ExpandedProgVaccUptake_SpreadsheetName = "../../Data/VaccUptake/WholePopnVaccUptakeCalc/AlternativeVaccSchedules/VaccUptakePerSeason_AllPopnHighVaccUptake_EMH.xlsx"

    #Declare sheetnames within files of interest
    SheetNames = ["2012_2013","2013_2014","2014_2015","2015_2016","2016_2017","2017_2018"]

    #Specify relevant ranges
    if RiskGroupFlag == 0
        FieldRange = "D2:D102"
    elseif RiskGroupFlag == 1
        FieldRange = "C2:C102"
    end

    #Assign rows of uptake array based on JobIdx (read from text file)
    VaccAgeBoundsData = readdlm("VaccUptakeSweepFiles/VaccUptakeAgeSweep_BoundsData_EntireAgeSpaceVers.txt",',')
    YoungClusterUB = VaccAgeBoundsData[JobID,1] #First entry the upper bound for young age cluster
    ElderClusterLB = VaccAgeBoundsData[JobID,2]   #Second entry for lower bound of elder age cluster

    println("JobID: $JobID")
    println("YoungClusterUB: $YoungClusterUB")
    println("ElderClusterLB: $ElderClusterLB")

    #Initialise storage array. Column per influenza season. Row per single year age bracket.
    VaccUptake = zeros(AgeGrpNum,length(SheetNames))

    #Iterate through each entry
    for ii = 1:length(SheetNames)
        #Import vaccine uptake data
        AtRiskOnlyProg_VaccUptakePercentages = XLSX.readdata(AtRiskOnlyVaccUptake_SpreadsheetName,"$(SheetNames[ii])",FieldRange)
        ExpandedProg_VaccUptakePercentages = XLSX.readdata(ExpandedProgVaccUptake_SpreadsheetName,"$(SheetNames[ii])",FieldRange)

        # Collate input data into Arrays & convert from percentages to proportions
        AtRiskOnlyProg_VaccUptakePropns = Array{Float64, 2}(AtRiskOnlyProg_VaccUptakePercentages/100)
        ExpandedProg_VaccUptakePropns = Array{Float64, 2}(ExpandedProg_VaccUptakePercentages/100)

        #Populate vaccine uptake array for current season iteration
        VaccUptakeBySeasonArray = AtRiskOnlyProg_VaccUptakePropns #Initial data. Match at risk only programme

        if isnan(YoungClusterUB)  #Elder cluster only
            ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
            VaccUptakeBySeasonArray[(ElderClusterLB+1):end] = ExpandedProg_VaccUptakePropns[(ElderClusterLB+1):end] #Amend elder age cluster rows
        elseif isnan(ElderClusterLB) #Young ages only
            YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
            VaccUptakeBySeasonArray[1:(YoungClusterUB+1)] = ExpandedProg_VaccUptakePropns[1:(YoungClusterUB+1)] #Amend young age cluster rows
        else #Both ages
            YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
            ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
            VaccUptakeBySeasonArray[1:(YoungClusterUB+1)] = ExpandedProg_VaccUptakePropns[1:(YoungClusterUB+1)] #Amend young age cluster rows
            VaccUptakeBySeasonArray[(ElderClusterLB+1):end] = ExpandedProg_VaccUptakePropns[(ElderClusterLB+1):end] #Amend elder age cluster rows
        end

        #Assign to output array
        VaccUptake[:,ii] = VaccUptakeBySeasonArray
    end

    #println("VaccUptake: $(VaccUptake[51,:])")

return VaccUptake
end

@everywhere function VaccEoSUptake_ExpandSchoolAndElderly(JobID,AgeGrpNum,RiskGroupFlag)
#Inputs:
#   JobID - (integer) Used to determine how VaccUptake array should be populated
#   AgeGrpNum - (integer) Number of single year age-brackets in use
#   RiskGroupFlag - (Flag variable) 0: Low risk; 1: At risk.
# Outputs:
#   VaccUptake - (tuple) Uptake propn at end of season per age. Tuple entry per influenza season

    #State relevant data file names
    AtRiskOnlyVaccUptake_SpreadsheetName = "../../Data/VaccUptake/WholePopnVaccUptakeCalc/AlternativeVaccSchedules/VaccUptakePerSeason_AtRiskSchedule_EMH.xlsx"
    ExpandedProgVaccUptake_SpreadsheetName = "../../Data/VaccUptake/WholePopnVaccUptakeCalc/AlternativeVaccSchedules/VaccUptakePerSeason_AtRiskPlusElderly&0-20Schedule_EMH.xlsx"

    #Declare sheetnames within files of interest
    SheetNames = ["2012_2013","2013_2014","2014_2015","2015_2016","2016_2017","2017_2018"]

    #Specify relevant ranges
    if RiskGroupFlag == 0
        FieldRange = "D2:D102"
    elseif RiskGroupFlag == 1
        FieldRange = "C2:C102"
    end

    #Assign rows of uptake array based on JobIdx (read from text file)
    VaccAgeBoundsData = readdlm("VaccUptakeSweepFiles/VaccUptakeAgeSweep_BoundsData.txt",',')
    YoungClusterUB = VaccAgeBoundsData[JobID,1] #First entry the upper bound for young age cluster
    ElderClusterLB = VaccAgeBoundsData[JobID,2]   #Second entry for lower bound of elder age cluster

    println("JobID: $JobID")
    println("YoungClusterUB: $YoungClusterUB")
    println("ElderClusterLB: $ElderClusterLB")

    #Initialise storage array. Column per influenza season. Row per single year age bracket.
    VaccUptake = zeros(AgeGrpNum,length(SheetNames))

    #Iterate through each entry
    for ii = 1:length(SheetNames)
        #Import vaccine uptake data
        AtRiskOnlyProg_VaccUptakePercentages = XLSX.readdata(AtRiskOnlyVaccUptake_SpreadsheetName,"$(SheetNames[ii])",FieldRange)
        ExpandedProg_VaccUptakePercentages = XLSX.readdata(ExpandedProgVaccUptake_SpreadsheetName,"$(SheetNames[ii])",FieldRange)

        # Collate input data into Arrays & convert from percentages to proportions
        AtRiskOnlyProg_VaccUptakePropns = Array{Float64, 2}(AtRiskOnlyProg_VaccUptakePercentages/100)
        ExpandedProg_VaccUptakePropns = Array{Float64, 2}(ExpandedProg_VaccUptakePercentages/100)

        #Populate vaccine uptake array for current season iteration
        VaccUptakeBySeasonArray = AtRiskOnlyProg_VaccUptakePropns #Initial data. Match at risk only programme

        if isnan(YoungClusterUB)  #Elder cluster only
            ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
            VaccUptakeBySeasonArray[(ElderClusterLB+1):end] = ExpandedProg_VaccUptakePropns[(ElderClusterLB+1):end] #Amend elder age cluster rows
        elseif isnan(ElderClusterLB) #Young ages only
            YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
            VaccUptakeBySeasonArray[1:(YoungClusterUB+1)] = ExpandedProg_VaccUptakePropns[1:(YoungClusterUB+1)] #Amend young age cluster rows
        else #Both ages
            YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
            ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
            VaccUptakeBySeasonArray[1:(YoungClusterUB+1)] = ExpandedProg_VaccUptakePropns[1:(YoungClusterUB+1)] #Amend young age cluster rows
            VaccUptakeBySeasonArray[(ElderClusterLB+1):end] = ExpandedProg_VaccUptakePropns[(ElderClusterLB+1):end] #Amend elder age cluster rows
        end

        #Assign to output array
        VaccUptake[:,ii] = VaccUptakeBySeasonArray
    end

    #println("VaccUptake: $(VaccUptake[51,:])")

return VaccUptake
end

###################################
### Adjust risk group infected ####
###################################

@everywhere function CorrectRiskGrpInfected_MultiJobVers(LowRiskInfDataFile,
                                AtRiskInfDataFile,
                                AllRiskInfDataFile,
                                RiskGrpPropnData)
                                #JobID)
#Inputs:
# LowRiskInfDataFile,AtRiskInfDataFile,AllRiskInfDataFile - (strings) Datafiles containing epidemiological data
# RiskGrpPropnData - (tuple) Gives proportion of age group that is (i) at risk, (ii) at-risk. Row by influenza season, column by year of age.
# JobID - (integer) ID for vaccine scheme under consideration. Different vaccine uptake & therefore different epi outcomes for each job.

#Outputs:
#PropnLowRiskInfectedCorrected,PropnAtRiskInfectedCorrected - (ctuple) Infected counts. Row by flu season, column by age
#PropnLowRiskAscertainCorrected,PropnAtRiskAscertainCorrected - (tuple) Flu attributed GP visit counts. Row by flu season, column by age.



#Load epi data from files
LowRiskEpiDataFile = matopen(LowRiskInfDataFile)
AtRiskEpiDataFile = matopen(AtRiskInfDataFile)
AllRiskEpiDataFile = matopen(AllRiskInfDataFile)

#Get variable names in each epi file
LowRiskVarNames = names(LowRiskEpiDataFile)
AtRiskVarNames = names(AtRiskEpiDataFile)
AllRiskVarNames = names(AllRiskEpiDataFile)

#Load data. Using separate MAT file per JobID. Thus, only one varaible per MAT file
LowRiskEpiData = read(LowRiskEpiDataFile, LowRiskVarNames[1])
AtRiskEpiData = read(AtRiskEpiDataFile, AtRiskVarNames[1])
AllRiskEpiData = read(AllRiskEpiDataFile, AllRiskVarNames[1])

close(LowRiskEpiDataFile)
close(AtRiskEpiDataFile)
close(AllRiskEpiDataFile)

#Extract relevant tuple entries for infection counts
AllRiskInfectedData_ByStrain = AllRiskEpiData[:,2]
LowRiskInfectedData_ByStrain = LowRiskEpiData[:,1]
AtRiskInfectedData_ByStrain = AtRiskEpiData[:,1]

#Extract relevant tuple entries for ascertained counts
AllRiskAscertainedData_ByStrain = AllRiskEpiData[:,1]
LowRiskAscertainedData_ByStrain = LowRiskEpiData[:,2]
AtRiskAscertainedData_ByStrain = AtRiskEpiData[:,2]

#Get number of evaluations to be performed
NumberOfEvals = length(AtRiskInfectedData_ByStrain)

#Initialise storage cells
PropnLowRiskInfectedCorrected = Array{Array{Float64,3},1}(undef,NumberOfEvals)
PropnAtRiskInfectedCorrected = Array{Array{Float64,3},1}(undef,NumberOfEvals)
PropnLowRiskAscertainCorrected = Array{Array{Float64,3},1}(undef,NumberOfEvals)
PropnAtRiskAscertainCorrected = Array{Array{Float64,3},1}(undef,NumberOfEvals)

#Disaggregate RiskGrpPropnData
LowRiskGrpPropn = RiskGrpPropnData[1]
AtRiskGrpPropn = RiskGrpPropnData[2]

#For each age, get low risk ascertainment rates
for EvalIdx = 1:NumberOfEvals

    #Get indexing values for epi data
    #CurrentEval_EpiParamIdx = EpiParamIdx(EvalIdx);
    CurrentEval_EpiParamIdx = EvalIdx

    #Assign relevant infected and ascertain counts to variables
    AtRiskInfectedData_ByStrain_CurrentEval = AtRiskInfectedData_ByStrain[CurrentEval_EpiParamIdx]
    LowRiskInfectedData_ByStrain_CurrentEval = LowRiskInfectedData_ByStrain[CurrentEval_EpiParamIdx]

    AtRiskAscertainedData_ByStrain_CurrentEval = AtRiskAscertainedData_ByStrain[CurrentEval_EpiParamIdx]
    LowRiskAscertainedData_ByStrain_CurrentEval = LowRiskAscertainedData_ByStrain[CurrentEval_EpiParamIdx]

    #Get infected counts by season, strain and age
    #Find discrepency between population measure and weighted sum of risk
    #group counts
    AtRiskWeightedContribution_InfCount = AtRiskGrpPropn[4:end,:].*AtRiskInfectedData_ByStrain_CurrentEval
    LowRiskWeightedContribution_InfCount = LowRiskGrpPropn[4:end,:].*LowRiskInfectedData_ByStrain_CurrentEval

    GrpSumTotal_InfCount = AtRiskWeightedContribution_InfCount .+ LowRiskWeightedContribution_InfCount
    ContributionRatio_InfCount_AtRiskVsLowRisk = AtRiskWeightedContribution_InfCount./GrpSumTotal_InfCount

    OverallTotal_InfCount = AllRiskInfectedData_ByStrain[EvalIdx]

    Discrep_InfCount = OverallTotal_InfCount .- GrpSumTotal_InfCount #Get difference in values

    #Get ascertained counts by season, strain and age
    #Find discrepency between population measure and weighted sum of risk
    #group counts
    AtRiskWeightedContribution_Ascertained = AtRiskGrpPropn[4:end,:].*AtRiskAscertainedData_ByStrain_CurrentEval
    LowRiskWeightedContribution_Ascertained = LowRiskGrpPropn[4:end,:].*LowRiskAscertainedData_ByStrain_CurrentEval

    GrpSumTotal_Ascertained = AtRiskWeightedContribution_Ascertained .+ LowRiskWeightedContribution_Ascertained
    ContributionRatio_Ascertained_AtRiskVsLowRisk = AtRiskWeightedContribution_Ascertained./GrpSumTotal_Ascertained

    OverallTotal_Ascertained = AllRiskInfectedData_ByStrain[EvalIdx]

    Discrep_Ascertained = OverallTotal_Ascertained .- GrpSumTotal_Ascertained  #Get difference in values

    #Compute adjustments for each risk group
    AtRiskAdjust_InfCount = Discrep_InfCount.*ContributionRatio_InfCount_AtRiskVsLowRisk
    LowRiskAdjust_InfCount = Discrep_InfCount .- AtRiskAdjust_InfCount

    AtRiskAdjust_Ascertained = Discrep_Ascertained.*ContributionRatio_Ascertained_AtRiskVsLowRisk
    LowRiskAdjust_Ascertained = Discrep_Ascertained .- AtRiskAdjust_Ascertained

    #Apply adjustments. Divide adjustment by risk group proportion to attain risk group specific adjustment
    #Due to AtRiskAdjust_InfCount = AbsoluteAdjustAtRisk.*AtRiskGrpPropn[4:end,:]
    #Therefore, AbsoluteAdjustAtRisk = AtRiskAdjust_InfCount./AtRiskGrpPropn[4:end,:]
    PropnAtRiskInfectedCorrected[EvalIdx] = AtRiskInfectedData_ByStrain_CurrentEval .+ (AtRiskAdjust_InfCount./AtRiskGrpPropn[4:end,:])
    PropnLowRiskInfectedCorrected[EvalIdx] = LowRiskInfectedData_ByStrain_CurrentEval .+ (LowRiskAdjust_InfCount./LowRiskGrpPropn[4:end,:])

    PropnAtRiskAscertainCorrected[EvalIdx] = AtRiskAscertainedData_ByStrain_CurrentEval .+ (AtRiskAdjust_Ascertained./AtRiskGrpPropn[4:end,:])
    PropnLowRiskAscertainCorrected[EvalIdx] = LowRiskAscertainedData_ByStrain_CurrentEval .+ (LowRiskAdjust_Ascertained./LowRiskGrpPropn[4:end,:])

    # PropnAtRiskInfectedCorrected[EvalIdx] = OverallTotal_InfCount.*ContributionRatio_InfCount_AtRiskVsLowRisk
    # PropnLowRiskInfectedCorrected[EvalIdx] = OverallTotal_InfCount.*(1 .- ContributionRatio_InfCount_AtRiskVsLowRisk)
    #
    # PropnAtRiskAscertainCorrected[EvalIdx] = OverallTotal_Ascertained.*ContributionRatio_Ascertained_AtRiskVsLowRisk
    # PropnLowRiskAscertainCorrected[EvalIdx] = OverallTotal_Ascertained.*(1 .- ContributionRatio_Ascertained_AtRiskVsLowRisk)

    #Sanity check
    RevisedGrpSumTotal_InfCount = AtRiskInfectedData_ByStrain_CurrentEval .+ (AtRiskAdjust_InfCount./AtRiskGrpPropn[4:end,:]) .+
        LowRiskInfectedData_ByStrain_CurrentEval .+ (LowRiskAdjust_InfCount./LowRiskGrpPropn[4:end,:])

    RevisedDiscrep_InfCount = OverallTotal_InfCount .- RevisedGrpSumTotal_InfCount #Get difference in values
    #println("RevisedDiscrep_InfCount: $RevisedDiscrep_InfCount")
end

return PropnLowRiskInfectedCorrected::Array{Array{Float64,3},1},
    PropnAtRiskInfectedCorrected::Array{Array{Float64,3},1},
    PropnLowRiskAscertainCorrected::Array{Array{Float64,3},1},
    PropnAtRiskAscertainCorrected::Array{Array{Float64,3},1}

end

###################################
### Risk group GP visit rate ######
###################################

@everywhere function RiskGrpAscertainCalc_MultiJobVers(PropnLowRiskInfectedCorrected,
                                PropnAtRiskInfectedCorrected,
                                AllRiskInfDataFile,
                                PopnDataTuple,
                                EpiParamIdx,
                                AtRiskScaleFactor)
                                #JobID)

#Inputs:
# PropnLowRiskInfectedCorrected,PropnAtRiskInfectedCorrected - (tuples) Proportion of risk group infected (age specific)
# AllRiskInfDataFile - (strings) Datafiles containing epidemiological data
# PopnDataTuple - (tuple) Three entries. (i) Low risk; (ii) At risk; (iii) All risk.
# EpiParamIdx - (vector, integers) Run IDs from epidemiological runs that should be accessed
# AtRiskScaleFactor - (scalar,float) Linear scaling for at risk group
#                               ascertainent relative to low risk group
# JobID - (integer) ID for vaccine scheme under consideration. Different vaccine uptake & therefore different epi outcomes for each job.

#Outputs: LowRiskAscertainRates,AtRiskAscertainRates - (vector, floats) By age, ascertainment rate for specified risk group

#Disaggregate PopnDataCell
PopnData_LowRisk = PopnDataTuple[1]
PopnData_AtRisk = PopnDataTuple[2]
PopnData_AllRisk = PopnDataTuple[3]

#Load epi data from files
AllRiskEpiDataFile = matopen(AllRiskInfDataFile)
AllRiskVarNames = names(AllRiskEpiDataFile)
AllRiskEpiData = read(AllRiskEpiDataFile, AllRiskVarNames[1])
close(AllRiskEpiDataFile)
AllRiskAscertainedData_ByStrain = AllRiskEpiData[:,1]
AllRiskInfectedData_ByStrain = AllRiskEpiData[:,2]

#Get number of evaluations to be performed
NumberOfEvals = length(EpiParamIdx)

#Initialise storage cells
LowRiskAscertainRates = Array{Array{Float64,2},1}(undef,NumberOfEvals)
AtRiskAscertainRates = Array{Array{Float64,2},1}(undef,NumberOfEvals)

#For each age, get low risk ascertainment rates
for EvalIdx = 1:NumberOfEvals

    #Get indexing values for epi data and at risk group ascertainment
    #scaling
    CurrentEval_EpiParamIdx = EpiParamIdx[EvalIdx]
    CurrentEval_AtRiskScaleFactor = AtRiskScaleFactor[EvalIdx]

    #Get total case counts by
    #(i) summing over all strains
    #(ii) dropping singleton dimensions
    #(iii) multiplying by population counts
    #To be used in ascertainment calculation
    AllRiskTotalAscertainSum = sum(AllRiskAscertainedData_ByStrain[CurrentEval_EpiParamIdx],dims=3)
    AllRiskTotalAscertainDropDims = dropdims(AllRiskTotalAscertainSum; dims=3)
    AllRiskTotalAscertain = AllRiskTotalAscertainDropDims.*PopnData_AllRisk

    LowRiskTotalInfectedSum = sum(PropnLowRiskInfectedCorrected[CurrentEval_EpiParamIdx],dims=3)
    LowRiskTotalInfectedDropDims = dropdims(LowRiskTotalInfectedSum; dims=3)
    LowRiskTotalInfected = LowRiskTotalInfectedDropDims.*PopnData_LowRisk

    AtRiskTotalInfectedSum = sum(PropnAtRiskInfectedCorrected[CurrentEval_EpiParamIdx],dims=3)
    AtRiskTotalInfectedDropDims = dropdims(AtRiskTotalInfectedSum; dims=3)
    AtRiskTotalInfected = AtRiskTotalInfectedDropDims.*PopnData_AtRisk

    #Fix number of seasons variable and age groups. Matches RowsxCols of AllRiskTotalAscertain
    NumOfSeasons = size(AllRiskTotalAscertain,1)
    AgeGrpNum = size(AllRiskTotalAscertain,2)

    #Work through each influenza season one by one
    #Get season-specific ascertainment rates using
    #Ascertainment^LR = (Overall popn ascertained cases) / (Total LR cases + (Total AtRisk cases*AtRiskScaleFactor))
    CurrentSeason_LowRiskAscertainRates = zeros(NumOfSeasons,AgeGrpNum)
    for SeasonIdx = 1:NumOfSeasons
        CurrentSeason_LowRiskAscertainRates[SeasonIdx,:] =
            AllRiskTotalAscertain[SeasonIdx,:] ./ (LowRiskTotalInfected[SeasonIdx,:] .+
                        (CurrentEval_AtRiskScaleFactor.*AtRiskTotalInfected[SeasonIdx,:]))
    end

    #Compute low risk ascertainment
    LowRiskAscertainRates[EvalIdx] = CurrentSeason_LowRiskAscertainRates

    #Compute at risk group ascertainent rates
    AtRiskAscertainRates[EvalIdx] = CurrentEval_AtRiskScaleFactor.*LowRiskAscertainRates[EvalIdx]

    ### Sanity check ###
    #Check for difference in risk group specific estimates compared to population level estimates
    AllRiskTotalInfectedSum = sum(AllRiskInfectedData_ByStrain[CurrentEval_EpiParamIdx],dims=3)
    AllRiskTotalInfectedDropDims = dropdims(AllRiskTotalInfectedSum; dims=3)
    AllRiskTotalInfected = AllRiskTotalInfectedDropDims.*PopnData_AllRisk

    Difference_Inf = AllRiskTotalInfected .- (LowRiskTotalInfected .+ AtRiskTotalInfected)
    println("Difference_Inf sum: $(sum(Difference_Inf))")

    GrpRiskAscertainment = (CurrentSeason_LowRiskAscertainRates.*LowRiskTotalInfected) .+ (AtRiskAscertainRates[EvalIdx].*AtRiskTotalInfected)
    Difference_Ascertain = AllRiskTotalAscertain - GrpRiskAscertainment
    println("Difference_Ascertain sum: $(sum(Difference_Ascertain))")

    ######################
end

return LowRiskAscertainRates::Array{Array{Float64,2},1},AtRiskAscertainRates::Array{Array{Float64,2},1}

end

#################
### GP visits ###
#################

#Return point estimate
@everywhere function GPAscertainParamGen_ReturnPtEstimate(RndNumGenInput,TotalParamSetsEval,AgeGrpNum)
#Inputs:
# RndNumGenInput - Parameters used for generating the random numbers
# TotalParamSetsEval - (scalar, integer) Number of runs carried out in PSA
# AgeGrpNum - (integer)

#Outputs:
# FebrileTriangularDist_r - (vector) Febrile case parameter value per parameter set evaluation (age and strain independent)
# GPFluAscertain - (tuple) Within each tuple, entry per age group of GP visit
#                                        propensity (given infected & symptomatic)

    #Create variable for febrile cases
    modeFebrile = RndNumGenInput[1]::Float64

    #Create variable for GP consultation propensity
    modeConsultPropens = RndNumGenInput[2]::Float64

    #Generate random numbers
    Febrile_r = modeFebrile.*ones(1,TotalParamSetsEval)
    ConsultPropens_r = modeConsultPropens.*ones(1,TotalParamSetsEval)

    #Expand to cover stated number of age groups (AgeGrpNum)
    GPFluAscertain= Array{Array{Float64}}(undef,TotalParamSetsEval) #Initialise tuple
    for ii = 1:TotalParamSetsEval
        GPFluAscertain[ii] = Febrile_r[ii].*ConsultPropens_r[ii].*ones(AgeGrpNum,1) #Assign to output variable
    end

return  Febrile_r::Array{Float64,2}, GPFluAscertain::Array{Array{Float64}}
end


#Triangular distribution
@everywhere function GPAscertainParamGen_TriangularDist(RndNumGenInput,TotalParamSetsEval,AgeGrpNum)
#Inputs:
# RndNumGenInput - Parameters used for generating the random numbers
# TotalParamSetsEval - (scalar, integer) Number of runs carried out in PSA
# AgeGrpNum - (integer)

#Outputs:
# FebrileTriangularDist_r - (vector) Febrile case parameter value per parameter set evaluation (age and strain independent)
# GPFluAscertain - (tuple) Within each tuple, entry per age group of GP visit
#                                        propensity (given infected & symptomatic)

    #Create triangular probability distribution objects for febrile cases
    lowerFebrile = RndNumGenInput[1,1]::Float64
    modeFebrile = RndNumGenInput[1,2]::Float64
    upperFebrile = RndNumGenInput[1,3]::Float64
    pdFebrile = TriangularDist(lowerFebrile,upperFebrile,modeFebrile)
    #TriangularDist(a,b,c)
    #The *triangular distribution* with lower limit `a`, upper limit `b` and mode `c` has probability density function

    #Create triangular probability distribution objects for GP consultation propensity
    lowerConsultPropens = RndNumGenInput[2,1]::Float64
    modeConsultPropens = RndNumGenInput[2,2]::Float64
    upperConsultPropens = RndNumGenInput[2,3]::Float64
    pdConsultPropens = TriangularDist(lowerConsultPropens,upperConsultPropens,modeConsultPropens)

    #Generate random numbers
    FebrileTriangularDist_r = rand(pdFebrile,1,TotalParamSetsEval)
    ConsultPropensTriangularDist_r  = rand(pdConsultPropens,1,TotalParamSetsEval)

    #println("FebrileTriangularDist_r: $FebrileTriangularDist_r")
    #println("ConsultPropensTriangularDist_r: $ConsultPropensTriangularDist_r")

    #Expand to cover stated number of age groups (AgeGrpNum)
    GPFluAscertain= Array{Array{Float64}}(undef,TotalParamSetsEval) #Initialise tuple
    for ii = 1:TotalParamSetsEval
        GPFluAscertain[ii] = FebrileTriangularDist_r[ii].*ConsultPropensTriangularDist_r[ii].*ones(AgeGrpNum,1) #Assign to output variable
    end

return  FebrileTriangularDist_r::Array{Float64,2}, GPFluAscertain::Array{Array{Float64}}
end


@everywhere function  GPAscertainParamGen_TriangularDistFebrile_RiskSpecificAscertain(RndNumGenInput,TotalParamSetsEval,AgeGrpNum)

#Inputs:
# RndNumGenInput - Parameters used for generating the random numbers
# TotalParamSetsEval - (scalar, integer) Number of runs carried out in PSA
# AgeGrpNum - (integer)

#Outputs:
# FebrileTriangularDist_r - (vector) Febrile case parameter value per parameter set evaluation (age and strain independent)
# GPFluAscertain - (tuple) Within each tuple, entry per age group of GP visit
#                                        propensity (given infected & symptomatic)

    #Create triangular probability distribution objects for febrile cases
    lowerFebrile = RndNumGenInput[1,1]::Float64
    modeFebrile = RndNumGenInput[1,2]::Float64
    upperFebrile = RndNumGenInput[1,3]::Float64
    pdFebrile = TriangularDist(lowerFebrile,upperFebrile,modeFebrile)
    #TriangularDist(a,b,c)
    #The *triangular distribution* with lower limit `a`, upper limit `b` and mode `c` has probability density function

    #Generate random numbers
    FebrileTriangularDist_r = rand(pdFebrile,1,TotalParamSetsEval)

    #Set GPFluAscertain as a cell variable
    GPFluAscertain = Array{Array{Float64,2},1}(undef,TotalParamSetsEval)

    return FebrileTriangularDist_r::Array{Float64,2}, GPFluAscertain::Array{Array{Float64,2},1}
end

#Return point estimate, with risk groups having distinct ascertainment rates
@everywhere function GPAscertainParamGen_ReturnPtEstimate_RiskSpecificAscertain(RndNumGenInput,TotalParamSetsEval,AgeGrpNum)
#Inputs:
# RndNumGenInput - Parameters used for generating the random numbers
# TotalParamSetsEval - (scalar, integer) Number of runs carried out in PSA
# AgeGrpNum - (integer)

#Outputs:
# FebrileTriangularDist_r - (vector) Febrile case parameter value per parameter set evaluation (age and strain independent)
# GPFluAscertain - (tuple) Within each tuple, entry per age group of GP visit
#                                        propensity (given infected & symptomatic)

    #Create variable for febrile cases
    modeFebrile = RndNumGenInput[1]::Float64

    #Generate random numbers
    Febrile_r = modeFebrile.*ones(1,TotalParamSetsEval)

    #Set GPFluAscertain as a cell variable
    GPFluAscertain = Array{Array{Float64,2},1}(undef,TotalParamSetsEval)

return  Febrile_r::Array{Float64,2},  GPFluAscertain::Array{Array{Float64,2},1}
end

##############################################################
### Health episode event likelihood (relative to GP visit) ###
##############################################################

#Return point estimate
@everywhere function HealthEpsRelLhoodParamGen_ReturnPtEstimate(ValsHealthEpsRelLhood_RndNumGenInput)
#%Inputs:
# ValsHealthEpsRelLhood - (Tuple) Parameters used for generating the returned values
#                              -> Entry 1: Mean values;
#                               -> Per entry, 3D array (row by age band,
#                                                       column by severity of health outcome,
#                                                       slice by strain)
#Outputs: ValsHealthEpsRelLhood - (3D array)

    ValsHealthEpsRelLhood_means = ValsHealthEpsRelLhood_RndNumGenInput

    ValsHealthEpsRelLhood = ValsHealthEpsRelLhood_means;  #3D array: Slice by strain; Row by age band; column by clinical event

    ValsHealthEpsRelLhood[:,:,2] = ValsHealthEpsRelLhood[:,:,1]  #Make type A subtypes have same costs
    ValsHealthEpsRelLhood[:,:,4] = ValsHealthEpsRelLhood[:,:,3]  #Make type B lineages have same costs

    return ValsHealthEpsRelLhood::Array{Float64,3}

end

@everywhere function  HealthEpsRelLhoodParamGen_Normal(ValsHealthEpsRelLhood_RndNumGenInput)
#Inputs:
# ValsHealthEpsRelLhood - (Tuple) Parameters used for generating the normally distributed random numbers
#                              -> Entry 1: Mean values; Entry 2; Standard deviations
#                               -> Per entry, 3D array (row by age band,
#                                                       column by severity of health outcome,
#                                                       slice by strain)
#Outputs: ValsHealthEpsRelLhood - (3D array)

    ValsHealthEpsRelLhood_means = ValsHealthEpsRelLhood_RndNumGenInput[1]
    ValsHealthEpsRelLhood_sd = ValsHealthEpsRelLhood_RndNumGenInput[2]

    #Build distribution
    d_EpsRelLhood = Normal.(ValsHealthEpsRelLhood_means,ValsHealthEpsRelLhood_sd)

    #Generate random numbers
    #2D array: Row by age band; two columns (non-hospitalised, hospitalised)
    ValsHealthEpsRelLhood = rand.(d_EpsRelLhood)

    ValsHealthEpsRelLhood[:,:,2] = ValsHealthEpsRelLhood[:,:,1]  #Make type A subtypes have same costs
    ValsHealthEpsRelLhood[:,:,4] = ValsHealthEpsRelLhood[:,:,3]  #Make type B lineages have same costs

return ValsHealthEpsRelLhood::Array{Float64,3}

end

@everywhere function  HealthEpsRelLhoodParamGen_Bootstrap(ValsHealthEpsRelLhood_RndNumGenInput)
#Inputs:
# ValsHealthEpsRelLhood - (Tuple) Parameters used for generating the normally distributed random numbers
#                              -> Entry 1: Mean values; Entry 2; Standard deviations
#                               -> Per entry, 3D array (row by age band,
#                                                       column by severity of health outcome,
#                                                       slice by strain)
#Outputs: ValsHealthEpsRelLhood - (3D array)

    #Set reference probability and reference population size
    RefProb = 0.1
    RefPopnSize = 100000

    #generate bootstrap replicate
    r_RefProb = rand.(Binomial.(RefPopnSize,RefProb))/RefPopnSize
    r_RefProb = 0.1

    #From event ratio point estimates get relative probability of each
    #health event
    ValsHealthEpsRelLhood_means = ValsHealthEpsRelLhood_RndNumGenInput
    ValsHealthEpsRelLhood_RelProb = RefProb.*ValsHealthEpsRelLhood_means;

    #Check relative probabilities do not exceed 1
    if sum(ValsHealthEpsRelLhood_RelProb .> 1) > 0
       error("Health episode relative probability exceeds 1")
    end

    #Obtain bootstrap ratio replicates
    r_HealthEpsRelProb = rand.(Binomial.(RefPopnSize,ValsHealthEpsRelLhood_RelProb))./RefPopnSize
    ValsHealthEpsRelLhood = r_HealthEpsRelProb./r_RefProb

    ValsHealthEpsRelLhood[:,:,2] = ValsHealthEpsRelLhood[:,:,1]  #Make type A subtypes have same costs
    ValsHealthEpsRelLhood[:,:,4] = ValsHealthEpsRelLhood[:,:,3]  #Make type B lineages have same costs

return ValsHealthEpsRelLhood::Array{Float64,3}

end

#############
### Costs ###
#############

#Return point estimate
@everywhere function  CostsParamGen_ReturnPtEstimate(ValsCosts_RndNumGenInput)
#Inputs:
#ValsCosts -  Parameters used for returning point estimates
#                              -> Entry 1: Mean values; Entry 2; Standard deviations
#                               -> Per entry, 3D array (row by age band,
#                                                       column by severity of health outcome,
#                                                       slice by strain)

#Outputs: ValsCosts - (3D array)

    ValsCosts_means = ValsCosts_RndNumGenInput
    ValsCosts = ValsCosts_means

    ValsCosts[:,:,2] = ValsCosts[:,:,1]  #Make type A subtypes have same costs
    ValsCosts[:,:,4] = ValsCosts[:,:,3]  #Make type B lineages have same costs

    #Events with point estimate of no cost,
    #Ensure all replicates will have zero cost.
    ValsCosts[ValsCosts_means.==0] .= 0

return ValsCosts::Array{Float64,3}

end

#LogNormal
@everywhere function  CostsParamGen_LogNormal(ValsCosts_RndNumGenInput)
#Inputs:
#ValsCosts - (Tuple) Parameters used for generating the normally distributed random numbers
#                              -> Entry 1: Mean values; Entry 2; Standard deviations
#                               -> Per entry, 3D array (row by age band,
#                                                       column by severity of health outcome,
#                                                       slice by strain)

#Outputs: ValsCosts - (3D array)

    ValsCosts_means = ValsCosts_RndNumGenInput[1]
    ValsCosts_var = ValsCosts_RndNumGenInput[2]

    mu = log.((ValsCosts_means.^2)./sqrt.(ValsCosts_var .+ ValsCosts_means.^2))
    sigma = sqrt.(log.((ValsCosts_var./(ValsCosts_means.^2)).+1))

    #Events with point estimate of no cost,
    #Ensure all replicates will have zero cost.
    mu[ValsCosts_means.==0] .= 0
    sigma[ValsCosts_means.==0] .= 1e-14 #Amend sigma value for those with zero mean cost so LogNormal function runs!

    #Build distribution
    d_Costs = LogNormal.(mu,sigma)

    #Generate random numbers
    #2D array: Row by age band; two columns (non-hospitalised, hospitalised)
    ValsCosts = rand.(d_Costs)

    ValsCosts[:,:,2] = ValsCosts[:,:,1]  #Make type A subtypes have same costs
    ValsCosts[:,:,4] = ValsCosts[:,:,3]  #Make type B lineages have same costs

    #Events with point estimate of no cost,
    #Ensure all replicates will have zero cost.
    ValsCosts[ValsCosts_means.==0] .= 0

return ValsCosts::Array{Float64,3}

end

#Gamma
@everywhere function  CostsParamGen_Gamma(ValsCosts_RndNumGenInput)
#Inputs:
#ValsCosts - (Tuple) Parameters used for generating the normally distributed random numbers
#                              -> Entry 1: Mean values; Entry 2; Standard deviations
#                               -> Per entry, 3D array (row by age band,
#                                                       column by severity of health outcome,
#                                                       slice by strain)

#Outputs: ValsCosts - (3D array)

    ValsCosts_means = ValsCosts_RndNumGenInput

    #Build distribution
    d_Costs = Gamma.(ValsCosts_means,1)

    #Generate random numbers
    #2D array: Row by age band; two columns (non-hospitalised, hospitalised)
    ValsCosts = rand.(d_Costs)

    ValsCosts[:,:,2] = ValsCosts[:,:,1]  #Make type A subtypes have same costs
    ValsCosts[:,:,4] = ValsCosts[:,:,3]  #Make type B lineages have same costs

    #Events with point estimate of no cost,
    #Ensure all replicates will have zero cost.
    ValsCosts[ValsCosts_means.==0] .= 0


return ValsCosts::Array{Float64,3}

end

###########################
### Non-fatal QALY loss ###
###########################

#Return point estimate
@everywhere function QALYlossParamGen_ReturnPtEstimate(ValsQALYloss_RndNumGenInput,StrainNum)
#Inputs:
# ValsQALYloss_RndNumGenInput - (Tuple) Parameters used for generating values
#                              -> Entry 1: Mean values;
#                               -> Per entry, 2D array (row by age band, column by severity of health outcome)
# StrainNum - (integer)

#Outputs: ValsQALYloss - (3D array) Row by age band, column by severity of health outcome,
#                                       slice by strain

    ValsQALYloss_means = ValsQALYloss_RndNumGenInput[1]

    ValsQALYloss_SingleStrain = ValsQALYloss_means  #2D array: Row by age band; two columns (non-hospitalised, hospitalised)

    ValsQALYloss = repeat(ValsQALYloss_SingleStrain,outer=(1,1,StrainNum)) #Expand to cater for multiple strains

    return ValsQALYloss::Array{Float64,3}

end

#Normal distributed uncertainty
@everywhere function  QALYlossParamGen_Normal(ValsQALYloss_RndNumGenInput,StrainNum)
#Inputs:
# ValsQALYloss_RndNumGenInput - (Tuple) Parameters used for generating the normally distributed random numbers
#                              -> Entry 1: Mean values; Entry 2; Standard deviations
#                               -> Per entry, 2D array (row by age band, column by severity of health outcome)
# StrainNum - (integer)

#Outputs: ValsQALYloss - (3D array) Row by age band, column by severity of health outcome,
#                                       slice by strain


    ValsQALYloss_means = ValsQALYloss_RndNumGenInput[1]
    ValsQALYloss_sd = ValsQALYloss_RndNumGenInput[2]

    #Build distribution
    d_QALYloss = Normal.(ValsQALYloss_means,ValsQALYloss_sd)

    #Generate random numbers
    #2D array: Row by age band; two columns (non-hospitalised, hospitalised)
    ValsQALYloss_SingleStrain = rand.(d_QALYloss)

    ValsQALYloss = repeat(ValsQALYloss_SingleStrain,outer=(1,1,StrainNum)) #Expand to cater for multiple strains

return ValsQALYloss::Array{Float64,3}
end

####################
### Vaccine cost ###
####################

#Return point estimate
@everywhere function VaccineCostParamGen_ReturnPtEstimate(VaccCost_RndNumGenInput,TotalParamSetsEval)
#Inputs:
# RndNumGenInput - Parameters used for generating the random numbers
# TotalParamSetsEval - (scalar, integer) Number of runs carried out in PSA

#Outputs: CostPerVaccineDose - (scalar) Generated vaccine cost
    CostPerVaccineDose = VaccCost_RndNumGenInput.*ones(TotalParamSetsEval,1)

    return CostPerVaccineDose

end

#Uniform
@everywhere function VaccineCostParamGen_UniformInt(RndNumGenInput,TotalParamSetsEval)
#Inputs:
# RndNumGenInput - Parameters used for generating the random numbers
# TotalParamSetsEval - (scalar, integer) Number of runs carried out in PSA

#Outputs: CostPerVaccineDose - (scalar) Generated vaccine cost
    CostPerVaccineDose = rand(1:RndNumGenInput,(TotalParamSetsEval))

return CostPerVaccineDose
end

#Triangular distribution
@everywhere function VaccineCostParamGen_TriangularDist(RndNumGenInput,TotalParamSetsEval)
#Inputs:
# RndNumGenInput - Parameters used for generating the random numbers
# TotalParamSetsEval - (scalar, integer) Number of runs carried out in PSA

#Outputs: CostPerVaccineDose - (scalar) Generated vaccine cost

    #Create a triangular probability distribution object
    lower = RndNumGenInput[1]::Float64
    mode = RndNumGenInput[2]::Float64
    upper = RndNumGenInput[3]::Float64

    #Build distribution
    d_VaccCosts = TriangularDist(lower,upper,mode)
    #TriangularDist(a,b,c)
    #The *triangular distribution* with lower limit `a`, upper limit `b` and mode `c` has probability density function

    #Generate random numbers
    #2D array: Row by age band; two columns (non-hospitalised, hospitalised)
    CostPerVaccineDose = rand(d_VaccCosts,TotalParamSetsEval)

return CostPerVaccineDose
end

@everywhere function  MultiParamSetEconEval(AgeGrpNum,StrainNum,
                        PropnInfData,PropnAscertainedData,GPVisitFromAscertainProbFlag,
                        PopnData,FebrileCaseParam,GPFluAscertain,
                           HealthEpsRelLhoodTuple,HealthEpsCostsTuple,HealthEpsNonFatalQALYLossTuple,
                         WeightsQoL,RemainLifeExpAllSeasons,
                          DiscountRate,VaccCostsData)

#Inputs:
#    AgeGrpNum - (scalar, interger) Number of age classes in use. If set to 1, using non-age structured model
#    StrainNum - (scalar, interger)
#    PropnInfData - (3D array) Row per season, column per age, slice per strain.
#    PropnAscertainedData - (3D array) Row per season, column per age, slice per strain.
#    GPVisitFromAscertainProbFlag - (flag variable) Determines GP visits should be obtained from ascertainment cases or ALL infected cases
#    PopnData - (2D array) Population size. Row per season. Column per age.
#    FebrileCaseParam - (scalar) Febrile case parameter value (age and strain independent)
#    GPFluAscertain - (row vector) GP visit proportion. Entry per age.
#    HealthEpsRelLhoodTuple - (Tuple) Entry 1: Relative episode likelihood (relative to GP visit) by age band
#                       -> Row per single year age group,
#                       -> Column per episode type (non-fatal inpatient, outpatient, fatal case in-hospital, fatal case out of hospital)
#                       -> Slice per strain type (Slice one&two - Influenza A; Slice three&four - Influenza B)
#                                   - Entry 2: Upper boundary of each age band
#    HealthEpsCostsTuple - (Tuple) Entry 1: Cost of health episode by age band
#                       -> Row per single year age group,
#                       -> Column per episode type (GP consultation, non-fatal inpatient, outpatient, fatal case in-hospital, fatal case out of hospital)
#                       -> Slice per strain type (Slice one&two - Influenza A; Slice three&four - Influenza B)
#                                   - Entry 2: Upper boundary of each age band
#    HealthEpsNonFatalQALYLossTuple - (Tuple) Entry 1: QALY Loss data by age band
#                       -> Row per single year age group,
#                       -> Column per non-fatal episode type (non-hospitalised, hospitalised)
#                       -> Slice per strain type (Slice one&two - Influenza A; Slice three&four - Influenza B)
#                                   - Entry 2: Upper boundary of each age band
#     WeightsQoL - (vector) Adjusted quality of life. Entry per year of age (0-100+)
#     RemainLifeExpAllSeasons - (2D array) Row per age. Column per reference year
#     DiscountRate - (two element vector,float) [CostDiscountRate,QALYDiscountRate]
#     VaccCostsData - (Tuple) Entry 1: Vaccine uptake ->Row per age, ->Column per season.
#                              Entry 2: Administration charge per vaccine dose
#                              Entry 3: Cost per vaccine dose

#Outputs:
#   In following tuples, arrays of 3 dimensions. Age x Season x Strain
#   ByYrOfAgeHealthEpsCounts, AllAgeHealthEpsCounts (7 event type entries)
#       - Counts by year of age AND overall: (i) GP visit; (ii) non-fatal hospital admissions, (iii) outpatient
#           visits, (iv) fatal cases (in hospital); (v) fatal cases (out of hosptial);
#           (vi) all infected cases (including asymptomatic); (vii) symptomatic cases only.
#   ByYrOfAgeHealthEpsCosts, AllAgeHealthEpsCounts(6 event type entries)
#       -   Monetary costs by year of age AND overall: (i) GP visit; (ii) non-fatal hospital admissions, (iii) outpatient
#           visits, (iv) fatal cases (in hospital only) (v) fatal cases (out of hospital)
#           (vi) Total cost by age (all treatments) and over all ages (all treatments)
#    ByYrOfAgeHealthEpsQALYloss ,AllAgeHealthEpsQALYloss (5 event type entries)
#       -  QALY losses by year of age AND overall: (i) Non-hospitalised; (ii) non-fatal hospital admissions,
#           (iii) fatal cases (in hospital only), (iv) fatal cases (out of hospital)
#           (v) Total cost by age (all treatments) and over all ages (all treatments)
#    VaccCosts (tuple, 3 entries)
#       - in each tuple entry, 2D array, Row by year of age. Column by season.
#       - (i) Admin charges; (ii) Non-admin charges; (iii) combined cost per dose
#    VaccDeployedTuple (tuple, 2 entries)
#       - Within each entry, Row by year of age. Column by season.
#       - (i) Aboslute count
#       - (ii) Discounted count

#--------------------------------------------------------------------------
#    HealthEpsRelLhood - (3D array)
#          -> Row per single year age group,
#          -> Column per episode type (non-fatal inpatient, outpatient, fatal case in-hospital, fatal case out of hospital)
#          -> Slice per strain type (Slice one&two - Influenza A; Slice three&four - Influenza B)
#--------------------------------------------------------------------------
ValsEpsLhood = HealthEpsRelLhoodTuple[1]::Array{Float64,3}
HealthEpsLhoodParams_UpperAgeBound = HealthEpsRelLhoodTuple[2]::Vector{Int64}
EpsLhoodEntityNum = size(ValsEpsLhood,2) #Number of healh episodes matches number of columns in ValsEpsLhood

#Call function to construct health episode likelihood array for single yrs of age
HealthEpsRelLhood = HealthEconParamByAgeConstruction(ValsEpsLhood,HealthEpsLhoodParams_UpperAgeBound,EpsLhoodEntityNum,
                                AgeGrpNum,StrainNum)

#--------------------------------------------------------------------------
#    HealthEpsCosts - (3D array)
#          -> Row per single year age group,
#          -> Column per episode type (GP consultation, non-fatal inpatient, outpatient, fatal case in-hospital, fatal case out of hospital)
#          -> Slice per strain type (Slice one&two - Influenza A; Slice three&four - Influenza B)
#--------------------------------------------------------------------------
ValsCosts = HealthEpsCostsTuple[1]::Array{Float64,3}
CostsParams_UpperAgeBound = HealthEpsCostsTuple[2]::Vector{Int64}
CostsEntityNum = size(ValsCosts,2); #Number of healh episodes matches number of columns in ValsCost

#Call function to construct costs array for single yrs of age
HealthEpsCosts = HealthEconParamByAgeConstruction(ValsCosts,CostsParams_UpperAgeBound,CostsEntityNum,
                                AgeGrpNum,StrainNum)

#--------------------------------------------------------------------------
#    HealthEpsNonFatalQALYLoss - (3D array)
#          -> Row per single year age group,
#          -> Column per non-fatal episode type (non-hospitalised, hospitalised)
#          -> Slice per strain type (Slice one&two - Influenza A; Slice three&four - Influenza B)
#--------------------------------------------------------------------------
ValsQALYloss = HealthEpsNonFatalQALYLossTuple[1]::Array{Float64,3}
QALYlossParams_UpperAgeBound = HealthEpsNonFatalQALYLossTuple[2]::Vector{Int64}
QALYEntityNum = size(ValsQALYloss,2) #Number of healh episodes matches number of columns in ValsQALYloss

#Call function to construct QALY losses array for single yrs of age
HealthEpsNonFatalQALYLoss = HealthEconParamByAgeConstruction(ValsQALYloss,QALYlossParams_UpperAgeBound,QALYEntityNum,
                                AgeGrpNum,StrainNum)

#--------------------------------------------------------------------------
### Call function
#--------------------------------------------------------------------------
ByYrOfAgeHealthEpsCounts,ByYrOfAgeHealthEpsCosts,ByYrOfAgeHealthEpsQALYloss,
    AllAgeHealthEpsCounts,AllAgeHealthEpsCosts,AllAgeHealthEpsQALYloss,
    VaccCosts,VaccDeployedTuple =
    DoHealthEconEval(PropnInfData,PropnAscertainedData,GPVisitFromAscertainProbFlag,
                    PopnData,FebrileCaseParam,GPFluAscertain,
                           HealthEpsRelLhood,HealthEpsCosts,HealthEpsNonFatalQALYLoss,
                         WeightsQoL,RemainLifeExpAllSeasons,
                          DiscountRate,VaccCostsData)

return ByYrOfAgeHealthEpsCounts,ByYrOfAgeHealthEpsCosts,ByYrOfAgeHealthEpsQALYloss,
    AllAgeHealthEpsCounts,AllAgeHealthEpsCosts,AllAgeHealthEpsQALYloss,
    VaccCosts,VaccDeployedTuple

end


#--------------------------------------------------------------------------
### Construct uncertainty distributions based on input flag options
#--------------------------------------------------------------------------
@everywhere function RunHealthEconEval(RiskGroupFlag,AgeStrucModelFlag,AgeGrpNum,StrainNum,TotalParamSetsEval,
                                DiscountRate,VaccUptake,AdminChargePerVaccDose,
                                PropnInfData,PropnAscertainedData,GPVisitFromAscertainProbFlag,
                                GrpSpecificAscertainRates,
                                UncertaintyDistFlag_GPAscertain,UncertaintyDistFlag_HealthEpsRelLhood,
                                UncertaintyDistFlag_HealthEpsCosts,UncertaintyDistFlag_QALYloss,
                                UncertaintyDistFlag_VaccCosts)

#--------------------------------------------------------------------------
### Import datasets that are consistent across all health economic
### evaluation runs
#--------------------------------------------------------------------------

#Load population level data
if RiskGroupFlag == 0
    PopnDistEngland20102018LR = readdlm("../../Data/HealthEconParamData/LowRiskPopn_CountPerAge.txt",',')
    if AgeStrucModelFlag == 1
        PopnData = PopnDistEngland20102018LR[4:end,:]
    else
        PopnData = sum(PopnDistEngland20102018LR[4:end,:],2) #For non-age structured model, get overall population
    end
elseif RiskGroupFlag == 1
     PopnDistEngland20102018HR = readdlm("../../Data/HealthEconParamData/AtRiskPopn_CountPerAge.txt",',')
    if AgeStrucModelFlag == 1
        PopnData = PopnDistEngland20102018HR[4:end,:]
    else
        PopnData = sum(PopnDistEngland20102018HR[4:end,:],2) #For non-age structured model, get overall population
    end
end

#Load quality of life weights data
WeightsQoL = readdlm("WeightsQoL_ByYrOfAge_JanssenStudies.txt",',')

#Imported RemainLifeExp data. Column by year. Row by age
RemainLifeExpTemp = readdlm("../../Data/DemographicData/ExpectationOfLifeEstimates/ONSExpectedLifeEng_PopnAvgCohort_SingleYrAgeGenderDist_Ages0to100.txt",',')


#Assign RemainLifeExp data based on risk group
if RiskGroupFlag == 0
    RemainLifeExpAllSeasons = RemainLifeExpTemp[:,13:18]
elseif RiskGroupFlag == 1
    AtRiskRemainingLifeExpScale = 1. #Scaling factor for remaining life expectancy of at risk individual (relative to low risk individual)
    RemainLifeExpAllSeasons = AtRiskRemainingLifeExpScale.*RemainLifeExpTemp[:,13:18]
end


#--------------------------------------------------------------------------
### GPFluAscertain (row vector) GP visit proportion. Entry per age.
#--------------------------------------------------------------------------

#Specify function to be used to generate GP ascertainment param values
if UncertaintyDistFlag_GPAscertain == 0 #Return point estimate
    GPAscertainParamGen = GPAscertainParamGen_ReturnPtEstimate

    #Peak value inputs
    RndNumGenInput = [0.406;0.10] #Row 1 - febrile cases; row 2 - GP visit propensity

elseif UncertaintyDistFlag_GPAscertain == 1 #Triangular distribution
    GPAscertainParamGen = GPAscertainParamGen_TriangularDist

    #Inputs for triangular probability distributed random number generator
    #Lower, Peak, Upper values
    RndNumGenInput = [0.309 0.406 0.513;0.05 0.10 0.12] #Row 1 - febrile cases; row 2 - GP visit propensity
elseif UncertaintyDistFlag_GPAscertain == 2 #Febrile symptomns from triangular distribution. Load GP ascertainment rates from file.
    GPAscertainParamGen = GPAscertainParamGen_TriangularDistFebrile_RiskSpecificAscertain

    #Inputs for triangular probability distributed random number generator
    #Lower, Peak, Upper values
    RndNumGenInput = [0.309 0.406 0.513;0.05 0.10 0.12] #Row 1 - febrile cases; row 2 - GP visit propensity

elseif UncertaintyDistFlag_GPAscertain == 3 #Febrile symptomns from point estimate. Load GP ascertainment rates from file.
    GPAscertainParamGen = GPAscertainParamGen_ReturnPtEstimate_RiskSpecificAscertain

    #Peak value inputs
    RndNumGenInput = [0.406;0.10] #Row 1 - febrile cases; row 2 - GP visit propensity
else
    error("Uncompatible UncertaintyDistFlag_GPAscertain value of $UncertaintyDistFlag_GPAscertain used. Must take value 0, 1, 2 or 3.")
end

FebrileCaseParam, GPFluAscertain = GPAscertainParamGen(RndNumGenInput,TotalParamSetsEval,AgeGrpNum)

println("FebrileCaseParam: $FebrileCaseParam")

#Overwrite GP ascertainment rate if using group specific values!
if UncertaintyDistFlag_GPAscertain == 2 || UncertaintyDistFlag_GPAscertain == 3

    #Based on risk group, get all relevant cell entries
    if RiskGroupFlag == 0 #Low risk
        GPFluAscertain = GrpSpecificAscertainRates[1]
    elseif RiskGroupFlag == 1 #At risk
        GPFluAscertain = GrpSpecificAscertainRates[2]
    end
end

#--------------------------------------------------------------------------
### Construct VaccCostsData tuple
#--------------------------------------------------------------------------
#Set vaccine price
if UncertaintyDistFlag_VaccCosts == 0 #Return point estimate
    VaccineCostParamGen = VaccineCostParamGen_ReturnPtEstimate

    #Set vaccine price distribution inputs
    RndNumGenInput = 0.
elseif UncertaintyDistFlag_VaccCosts == 1 #Uniform integer

    VaccineCostParamGen = VaccineCostParamGen_UniformInt

    #Set vaccine price distribution inputs
    RndNumGenInput = 1.
elseif UncertaintyDistFlag_VaccCosts == 2 #Triangular distribution
    VaccineCostParamGen = VaccineCostParamGen_TriangularDist

    #Set vaccine price distribution inputs
    RndNumGenInput = [10.,15.,20.]
else
    error("Uncompatible UncertaintyDistFlag_VaccCosts value of $UncertaintyDistFlag_VaccCosts used. Must take value 0, 1 or 2.")
end

CostPerVaccineDose = VaccineCostParamGen(RndNumGenInput,TotalParamSetsEval)
VaccineDoseCosts = [ones(TotalParamSetsEval)*AdminChargePerVaccDose CostPerVaccineDose]

#println("VaccineDoseCosts: $VaccineDoseCosts")

#Aggregate uptake data and cost per vaccine into tuple
#VaccCostsData = Array{Tuple{Array{Float64,2},Float64},1}(TotalParamSetsEval) #Initialise tuple
VaccCostsData = Array{Array{Any,1}}(undef,TotalParamSetsEval) #Initialise tuple
for ii = 1:TotalParamSetsEval
    VaccCostsData[ii] = [VaccUptake,AdminChargePerVaccDose,CostPerVaccineDose[ii]]
end

#--------------------------------------------------------------------------
### HealthEpsNonFatalQALYLossTuple
#--------------------------------------------------------------------------

#Declare ValsQALYloss_RndNumGenInput
ValsQALYloss_means = [0.00749 0.016; 0.0082 0.018]
ValsQALYloss_sd = [8.5e-4 1.8e-3; 1.8e-3 1.8e-3]

ValsQALYloss_RndNumGenInput = [ValsQALYloss_means,ValsQALYloss_sd]

#Specify age bounds
if AgeStrucModelFlag == 1
    QALYlossParams_UpperAgeBound = vec([16 100])
else
    QALYlossParams_UpperAgeBound = vec([0])
end

#Specify function to be used to generate QALY loss param values
if UncertaintyDistFlag_QALYloss == 0 #Return point estimate
    QALYlossParamGen = QALYlossParamGen_ReturnPtEstimate;
elseif UncertaintyDistFlag_QALYloss == 1 #Normal distribution
    QALYlossParamGen = QALYlossParamGen_Normal
else
    error("Uncompatible UncertaintyDistFlag_HealthEpsCosts value of $UncertaintyDistFlag_QALYloss used. Must take value 1.")
end

#Initialise Tuple and generate TotalParamSetsEval sets of parameters
HealthEpsNonFatalQALYLossTuple = Array{Any,2}(undef,TotalParamSetsEval,2) #Initialise tuple
for ii = 1:TotalParamSetsEval
    #Generate parameter values for non-fatal QALY loss for simn ii
    ValsQALYloss = QALYlossParamGen(ValsQALYloss_RndNumGenInput,StrainNum)

    #Print outputs. Can check correspondance of values between runs.
    println("ItrIdx: $(ii)")
    println("ValsQALYloss: $(ValsQALYloss)")

    #Assign to tuple
    HealthEpsNonFatalQALYLossTuple[ii] = [ValsQALYloss,QALYlossParams_UpperAgeBound]
end

#--------------------------------------------------------------------------
### Construct HealthEpsRelLhoodTuple
#--------------------------------------------------------------------------

#Specify age bounds
if AgeStrucModelFlag == 1
    HealthEpsLhoodParams_UpperAgeBound = vec([1 5 9 19 29 39 49 59 64 74 84 100])
else
    HealthEpsLhoodParams_UpperAgeBound = vec([0]) #Non-age structrued version
end

HealthEpsLhoodAgeBandNum = length(HealthEpsLhoodParams_UpperAgeBound)

### ValsHealthEpsRelLhood_RndNumGenInput - BEGIN ###

#Files containing inpatient, outpatient and in-hospital mortality event info
NonFatalInpatient_SpreadsheetName = "../../Data/HealthEconParamData/HES_NonFatalEmergencyAdmisProbAndCosts_MultiYrAvg_J10&J11_Prim&Sec.xlsx"
Outpatient_SpreadsheetName = "../../Data/HealthEconParamData/HES_OutpatientProbAndCosts_MultiYrAvg_J10&J11_Prim&Sec.xlsx"
FatalHosp_SpreadsheetName = "../../Data/HealthEconParamData/HES_DeathProbAndCosts_MultiYrAvg_J10&J11_Prim&Sec.xlsx"

#Ratio of mortality events occurring out of hospital compared to in-hospital
#Ages 50-64, 75% in hospital, 25% out of hospital
#Ages 65-74, 65% in hospital, 35% out of hospital
#Ages 75+, 50% in hospital, 50% out of hospital.
if AgeStrucModelFlag == 1
    OutOfHospFatalCaseRatio = [0;0;0;0;0;0;0;25/75;25/75;35/65;1;1]
else
    OutOfHospFatalCaseRatio = 0
end

#Specify fields in speadsheet containing relevant data
if RiskGroupFlag == 0
     if AgeStrucModelFlag == 1
         EpsLhoodFieldRange = "B2:B13"
     else
         EpsLhoodFieldRange = "B14:B14"
     end
 elseif RiskGroupFlag == 1
      if AgeStrucModelFlag == 1
         EpsLhoodFieldRange = "E2:E13"
     else
         EpsLhoodFieldRange = "E14:E14"
     end
 end

#Set up ValsEpsLhood_means array
HealthEpsEntityNum = 4
ValsHealthEpsRelLhood_means  = zeros(HealthEpsLhoodAgeBandNum,HealthEpsEntityNum,StrainNum) #Slice per strain type (Type A & Type B)
SheetNames = ["Type A","Type B"]
for ii = 1:length(SheetNames)
    NonFatalInpatientCostTemp = XLSX.readdata(NonFatalInpatient_SpreadsheetName,"$(SheetNames[ii])",EpsLhoodFieldRange)
    OutpatientCostTemp = XLSX.readdata(Outpatient_SpreadsheetName,"$(SheetNames[ii])",EpsLhoodFieldRange)
    FatalHospCostTemp = XLSX.readdata(FatalHosp_SpreadsheetName,"$(SheetNames[ii])",EpsLhoodFieldRange)

    if ii == 1 #Type A
        # Collate into Array, Assign to storage cell
        ValsHealthEpsRelLhood_means[:,1,1] = Array{Float64, 2}(NonFatalInpatientCostTemp)
        ValsHealthEpsRelLhood_means[:,2,1] = Array{Float64, 2}(OutpatientCostTemp)
        ValsHealthEpsRelLhood_means[:,3,1] = Array{Float64, 2}(FatalHospCostTemp)
        ValsHealthEpsRelLhood_means[:,4,1] = OutOfHospFatalCaseRatio

        # Collate into Array, Assign to storage cell
        ValsHealthEpsRelLhood_means[:,1,2] = Array{Float64, 2}(NonFatalInpatientCostTemp)
        ValsHealthEpsRelLhood_means[:,2,2] = Array{Float64, 2}(OutpatientCostTemp)
        ValsHealthEpsRelLhood_means[:,3,2] = Array{Float64, 2}(FatalHospCostTemp)
        ValsHealthEpsRelLhood_means[:,4,2] = OutOfHospFatalCaseRatio
    else #Type B
        # Collate into Array, Assign to storage cell
        ValsHealthEpsRelLhood_means[:,1,3] = Array{Float64, 2}(NonFatalInpatientCostTemp)
        ValsHealthEpsRelLhood_means[:,2,3] = Array{Float64, 2}(OutpatientCostTemp)
        ValsHealthEpsRelLhood_means[:,3,3] = Array{Float64, 2}(FatalHospCostTemp)
        ValsHealthEpsRelLhood_means[:,4,3] = OutOfHospFatalCaseRatio

        # Collate into Array, Assign to storage cell
        ValsHealthEpsRelLhood_means[:,1,4] = Array{Float64, 2}(NonFatalInpatientCostTemp)
        ValsHealthEpsRelLhood_means[:,2,4] = Array{Float64, 2}(OutpatientCostTemp)
        ValsHealthEpsRelLhood_means[:,3,4] = Array{Float64, 2}(FatalHospCostTemp)
        ValsHealthEpsRelLhood_means[:,4,4] = OutOfHospFatalCaseRatio
    end
end

#########################

### ValsHealthEpsRelLhood_RndNumGenInput - END ###

#Specify function to be used to generate health episode likelihood from
#uncertainty distribution
if UncertaintyDistFlag_HealthEpsRelLhood == 0 #Return point estimates
    #Specify inputs to economic parameter generator for sensitivity
    #analysis
    ValsHealthEpsRelLhood_RndNumGenInput = ValsHealthEpsRelLhood_means

    #Define function handle
    HealthEpsRelLhoodParamGen = HealthEpsRelLhoodParamGen_ReturnPtEstimate
elseif UncertaintyDistFlag_HealthEpsRelLhood == 1 #Normal distribution

    #Set up ValsCost standard deviation array
    ValsHealthEpsRelLhood_sd = 0.01*ones(HealthEpsLhoodAgeBandNum,HealthEpsEntityNum,StrainNum)

    #Aggregate means and sd params into a tuple
    ValsHealthEpsRelLhood_RndNumGenInput = [ValsHealthEpsRelLhood_means,ValsHealthEpsRelLhood_sd]

    #Define function handle
    HealthEpsRelLhoodParamGen = HealthEpsRelLhoodParamGen_Normal
elseif UncertaintyDistFlag_HealthEpsRelLhood == 2 #Bootstrap

    #Specify inputs to economic parameter generator for sensitivity analysis
    ValsHealthEpsRelLhood_RndNumGenInput = ValsHealthEpsRelLhood_means

    #Define function handle
    HealthEpsRelLhoodParamGen = HealthEpsRelLhoodParamGen_Bootstrap
else
    error("UncompatibleUncertaintyDistFlag_HealthEpsRelLhood value of $UncertaintyDistFlag_HealthEpsRelLhood used. Must take value 0, 1 or 2.")
end


#Initialise Tuple and generate TotalParamSetsEval sets of parameters
HealthEpsRelLhoodTuple = Array{Any,2}(undef,TotalParamSetsEval,2) #Initialise tuple
for ii = 1:TotalParamSetsEval
    #Generate parameter values for health episode likelihood for simn ii
    ValsHealthEpsRelLhood = HealthEpsRelLhoodParamGen(ValsHealthEpsRelLhood_RndNumGenInput)

    #Print outputs. Can check correspondance of values between runs.
    println("ItrIdx: $(ii)")
    println("ValsHealthEpsRelLhood: $(ValsHealthEpsRelLhood)")

    #Assign to tuple
    HealthEpsRelLhoodTuple[ii] = [ValsHealthEpsRelLhood,HealthEpsLhoodParams_UpperAgeBound]

end

#--------------------------------------------------------------------------
### Construct HealthEpsCostsTuple
#--------------------------------------------------------------------------

#Specify age bounds
if AgeStrucModelFlag == 1
    CostParams_UpperAgeBound = vec([1 5 9 19 29 39 49 59 64 74 84 100])
else
    CostParams_UpperAgeBound = vec([0]) #Non-age structured version
end

CostAgeBandNum = length(CostParams_UpperAgeBound)

### Declare ValsCost_RndNumGenInput - BEGIN ###

#GP visit costs
GPConsultCost = 37.40*ones(CostAgeBandNum)

#Out-of-hospital mortality event cost
FatalOutOfHospCost = zeros(CostAgeBandNum)

#File Names
NonFatalInpatient_SpreadsheetName = "../../Data/HealthEconParamData/HES_NonFatalEmergencyAdmisProbAndCosts_MultiYrAvg_J10&J11_Prim&Sec.xlsx"
Outpatient_SpreadsheetName = "../../Data/HealthEconParamData/HES_OutpatientProbAndCosts_MultiYrAvg_J10&J11_Prim&Sec.xlsx"
FatalHosp_SpreadsheetName = "../../Data/HealthEconParamData/HES_DeathProbAndCosts_MultiYrAvg_J10&J11_Prim&Sec.xlsx"

#Specify fields in speadsheet containing relevant data
if RiskGroupFlag == 0
    if AgeStrucModelFlag == 1
        CostFieldRange = "D2:D13"
    else
        CostFieldRange = "D14:D14"
    end
elseif RiskGroupFlag == 1
    if AgeStrucModelFlag == 1
        CostFieldRange = "G2:G13"
    else
        CostFieldRange = "G14:G14"
    end
end

#Set up ValsCost array
CostEntityNum = 5
ValsCost_means = zeros(CostAgeBandNum,CostEntityNum,StrainNum) #Slice per strain type (Type A & Type B)
SheetNames = ["Type A","Type B"]
for ii = 1:length(SheetNames)
    NonFatalInpatientCostTemp = XLSX.readdata(NonFatalInpatient_SpreadsheetName,"$(SheetNames[ii])",CostFieldRange)
    OutpatientCostTemp = XLSX.readdata(Outpatient_SpreadsheetName,"$(SheetNames[ii])",CostFieldRange)
    FatalHospCostTemp = XLSX.readdata(FatalHosp_SpreadsheetName,"$(SheetNames[ii])",CostFieldRange)

    if ii == 1 #Type A
        # Collate into Array, Assign to storage cell
        ValsCost_means[:,1,1] = GPConsultCost
        ValsCost_means[:,2,1] = Array{Float64, 2}(NonFatalInpatientCostTemp)
        ValsCost_means[:,3,1] = Array{Float64, 2}(OutpatientCostTemp)
        ValsCost_means[:,4,1] = Array{Float64, 2}(FatalHospCostTemp)
        ValsCost_means[:,5,1] = FatalOutOfHospCost

        # Collate into Array, Assign to storage cell
        ValsCost_means[:,1,2] = GPConsultCost
        ValsCost_means[:,2,2] = Array{Float64, 2}(NonFatalInpatientCostTemp)
        ValsCost_means[:,3,2] = Array{Float64, 2}(OutpatientCostTemp)
        ValsCost_means[:,4,2] = Array{Float64, 2}(FatalHospCostTemp)
        ValsCost_means[:,5,2] = FatalOutOfHospCost
    else #Type B
        # Collate into Array, Assign to storage cell
        ValsCost_means[:,1,3] = GPConsultCost
        ValsCost_means[:,2,3] = Array{Float64, 2}(NonFatalInpatientCostTemp)
        ValsCost_means[:,3,3] = Array{Float64, 2}(OutpatientCostTemp)
        ValsCost_means[:,4,3] = Array{Float64, 2}(FatalHospCostTemp)
        ValsCost_means[:,5,3] = FatalOutOfHospCost

        ValsCost_means[:,1,4] = GPConsultCost
        ValsCost_means[:,2,4] = Array{Float64, 2}(NonFatalInpatientCostTemp)
        ValsCost_means[:,3,4] = Array{Float64, 2}(OutpatientCostTemp)
        ValsCost_means[:,4,4] = Array{Float64, 2}(FatalHospCostTemp)
        ValsCost_means[:,5,4] = FatalOutOfHospCost
    end
end

### Declare ValsCost_RndNumGenInput - END ###

#Specify function to be used to generate health episode costs from
#uncertainty distribution
if UncertaintyDistFlag_HealthEpsCosts == 0 #Return point estimate
    #Assign means to input variable
    ValsCost_RndNumGenInput = ValsCost_means;

    #Define function handle
    CostsParamGen = CostsParamGen_ReturnPtEstimate;
elseif UncertaintyDistFlag_HealthEpsCosts == 1 #Lognormal distribution

    #Set up ValsCost standard deviation array
    #ValsCost_var = (10*ones(CostAgeBandNum,CostEntityNum,StrainNum)).^2
    ValsCost_var = zeros(CostAgeBandNum,CostEntityNum,StrainNum)

    ValsCost_var[:,1,:] .= 8.4^2 #Gp consultation
    ValsCost_var[:,5,:] .= 0 #Out of hosptital

    #Set in-patient and outpatient standard deviation to be specified fraction of mean
    sdFrac = 0.2246
    ValsCost_var[:,2:4,:] .= (ValsCost_means[:,2:4,:].*sdFrac).^2

    #Aggregate means and var params into a tuple
    ValsCost_RndNumGenInput = [ValsCost_means,ValsCost_var]

    #Define function handle
    CostsParamGen = CostsParamGen_LogNormal
elseif UncertaintyDistFlag_HealthEpsCosts == 2 #Gamma distribution
    #Aggregate means params into a tuple
    ValsCost_RndNumGenInput = ValsCost_means

    #Define function handle
    CostsParamGen = CostsParamGen_Gamma
else
    error("Uncompatible UncertaintyDistFlag_HealthEpsCosts value of $UncertaintyDistFlag_HealthEpsCosts used. Must take value 1.")
end

#Initialise Tuple and generate TotalParamSetsEval sets of parameters
HealthEpsCostsTuple = Array{Any,2}(undef,TotalParamSetsEval,2) #Initialise tuple
for ii = 1:TotalParamSetsEval
    #Generate parameter values for health episode costs for simn ii
    ValsCost = CostsParamGen(ValsCost_RndNumGenInput)

    #Print outputs. Can check correspondance of values between runs.
    println("ItrIdx: $(ii)")
    println("ValsCost: $(ValsCost)")

    #Assign to tuple
    HealthEpsCostsTuple[ii] = [ValsCost,CostParams_UpperAgeBound]
end

#--------------------------------------------------------------------------
### Perform health economic evaluation for multiple parameter sets
#--------------------------------------------------------------------------

#Initialise storage tuples
ByYrOfAgeHealthEpsCounts = Array{Array{Array{Float64}},1}(undef,TotalParamSetsEval)
ByYrOfAgeHealthEpsCosts = Array{Array{Array{Float64}},1}(undef,TotalParamSetsEval)
ByYrOfAgeHealthEpsQALYloss = Array{Array{Array{Float64}},1}(undef,TotalParamSetsEval)
AllAgeHealthEpsCounts = Array{Array{Array{Float64}},1}(undef,TotalParamSetsEval)
AllAgeHealthEpsCosts = Array{Array{Array{Float64}},1}(undef,TotalParamSetsEval)
AllAgeHealthEpsQALYloss = Array{Array{Array{Float64}},1}(undef,TotalParamSetsEval)
VaccCosts = Array{Array{Array{Float64}},1}(undef,TotalParamSetsEval)
VaccDeployedTuple = Array{Array{Array{Float64}},1}(undef,TotalParamSetsEval)

#Call health economic evaluation function
for ParamSetIdx = 1:TotalParamSetsEval
    ByYrOfAgeHealthEpsCounts[ParamSetIdx],ByYrOfAgeHealthEpsCosts[ParamSetIdx],ByYrOfAgeHealthEpsQALYloss[ParamSetIdx],
        AllAgeHealthEpsCounts[ParamSetIdx],AllAgeHealthEpsCosts[ParamSetIdx],
        AllAgeHealthEpsQALYloss[ParamSetIdx],VaccCosts[ParamSetIdx],VaccDeployedTuple[ParamSetIdx] =
        MultiParamSetEconEval(AgeGrpNum,StrainNum,
        PropnInfData[ParamSetIdx],PropnAscertainedData[ParamSetIdx],GPVisitFromAscertainProbFlag,
        PopnData,FebrileCaseParam[ParamSetIdx],GPFluAscertain[ParamSetIdx],
        HealthEpsRelLhoodTuple[ParamSetIdx],HealthEpsCostsTuple[ParamSetIdx],HealthEpsNonFatalQALYLossTuple[ParamSetIdx],
        WeightsQoL,RemainLifeExpAllSeasons,DiscountRate,VaccCostsData[ParamSetIdx])
    println("ParamSetIdx: $ParamSetIdx")
end

return ByYrOfAgeHealthEpsCounts,ByYrOfAgeHealthEpsCosts,ByYrOfAgeHealthEpsQALYloss,
            AllAgeHealthEpsCounts,AllAgeHealthEpsCosts,AllAgeHealthEpsQALYloss,
            VaccCosts,VaccDeployedTuple,VaccineDoseCosts

end


### MAIN SCRIPT BEGINS ###

function RunHealthEcon_VaccAgeSweep(ARGS)

#Take command line arguments, ARGS, assign to variable names
#To convert strings to numbers, use parse

#--------------------------------------------------------------------------
### Seed RNG.
#--------------------------------------------------------------------------
RNG_SeedNum = parse(Int64, ARGS[1])
Random.seed!(RNG_SeedNum)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
### Set ID for economic evaluation
#--------------------------------------------------------------------------
HealthEconID_RiskGrpSpecific = ARGS[2]
HealthEconID_RiskGrpIndep = ARGS[3]

#--------------------------------------------------------------------------
### Import group proportion data
#--------------------------------------------------------------------------

ByYrOfAgeRiskGrpPropnPerSeasonModifiedAges2to4EMH  = XLSX.readdata("../../Data/RiskGrpPropnsData/AtRiskPropnSpreadsheets/ByYrOfAge_RiskGrpPropnPerSeason_ModifiedAges2to4_EMH.xlsx","Sheet1","B11:CX19")

# Collate into Array
# Row per influenza season; Column per age class
AtRiskGrpPropn = Array{Float64, 2}(ByYrOfAgeRiskGrpPropnPerSeasonModifiedAges2to4EMH )
LowRiskGrpPropn = 1 .- AtRiskGrpPropn

#Aggregate risk group data into tuple
RiskGrpPropnData = [LowRiskGrpPropn,AtRiskGrpPropn]

#--------------------------------------------------------------------------
### Declare fixed constants
#--------------------------------------------------------------------------
StrainNum = 4
AgeStrucModelFlag = 1 #Indicator variable. 1 for age structured model. 0 for non-age structured model.
RiskGroupFlag = parse(Int64, ARGS[4]) #0: Low risk; 1: At-risk

if AgeStrucModelFlag == 1
    AgeGrpNum = 101
else
    AgeGrpNum = 1
end

#--------------------------------------------------------------------------
### Error checks for fixed constants
#--------------------------------------------------------------------------
if AgeStrucModelFlag !=0 && AgeStrucModelFlag !=1
   error("AgeStrucModelFlag has value $AgeStrucModelFlag. Must take value 0 or 1.")
end

if RiskGroupFlag !=0 && RiskGroupFlag !=1
   error("RiskGroupFlag has value $RiskGroupFlag. Must take value 0 or 1.")
end

#--------------------------------------------------------------------------
### Set number of runs to be carried out & number of strategies (jobs) to be analysed
#--------------------------------------------------------------------------
TotalParamSetsEval = parse(Int64, ARGS[5])
JobNum = parse(Int64, ARGS[6])

#--------------------------------------------------------------------------
### Generate epi params to be used
#--------------------------------------------------------------------------
#The first argument to the general rand function gives a "thing to sample from"
UncertaintyDistFlag = parse(Int64, ARGS[7])
if UncertaintyDistFlag == 1 #Use multiple particles
    EpiParamIdx = rand(1:10,TotalParamSetsEval,1)
elseif UncertaintyDistFlag == 2 #Use particle retaining lowest error
    EpiParamIdx = rand(1:1,TotalParamSetsEval,1)
else #Throw error
    error("UncertaintyDistFlag mus take value 1 or 2. UncertaintyDistFlag currently has value $UncertaintyDistFlag.")
end
println("EpiParamIdx: $EpiParamIdx")

#--------------------------------------------------------------------------
### Set discount rate
#--------------------------------------------------------------------------
CostDiscountRate = parse(Float64, ARGS[8])
QALYDiscountRate = parse(Float64, ARGS[9])

DiscountRate = [CostDiscountRate QALYDiscountRate]

#--------------------------------------------------------------------------
### Correct risk group infection counts to match population level!
#--------------------------------------------------------------------------

#When summing over risk groups, the infected total does not match the
#expected overall total (per age group)!

#Get population data for each risk group
PopnDistLowRiskFile = readdlm("../../Data/HealthEconParamData/LowRiskPopn_CountPerAge.txt",',')
PopnData_LowRisk = PopnDistLowRiskFile[4:end,:] #Get data for 2012/13 influenza season onwards

PopnDistAtRiskFile = readdlm("../../Data/HealthEconParamData/AtRiskPopn_CountPerAge.txt",',')
PopnData_AtRisk = PopnDistAtRiskFile[4:end,:]

PopnData_AllRisk = PopnData_LowRisk + PopnData_AtRisk #Recover all risk group population counts

PopnDataCell = [PopnData_LowRisk,PopnData_AtRisk,PopnData_AllRisk] #Concantenate into single cell


#--------------------------------------------------------------------------
### Generate group specific ascertainment rates (from population level
### data)
#--------------------------------------------------------------------------

#Generate scaling factor values, set mean and sd for normal distribution
if UncertaintyDistFlag == 1 #Use multiple particles
  GPAtRiskScaleFactor_mean = 1.51
  GPAtRiskScaleFactor_sd = 0.18

  #Build distribution
  d_GPAtRiskScaleFactor = Normal(GPAtRiskScaleFactor_mean,GPAtRiskScaleFactor_sd)

  #Generate random numbers
  #Vector: Row for each evaluation run
  AtRiskScaleFactor = rand(d_GPAtRiskScaleFactor,TotalParamSetsEval,1)

elseif UncertaintyDistFlag == 2 #Using a single particle. No variation in paramter distribution required.
  GPAtRiskScaleFactor_mean = 1.51
  GPAtRiskScaleFactor_sd = 0

  #Build distribution
  d_GPAtRiskScaleFactor = Normal(GPAtRiskScaleFactor_mean,GPAtRiskScaleFactor_sd)

  #Generate random numbers
  #Vector: Row for each evaluation run
  AtRiskScaleFactor = rand(d_GPAtRiskScaleFactor,TotalParamSetsEval,1)
else #Throw error
  error("UncertaintyDistFlag must take value 1 or 2. UncertaintyDistFlag currently has value $UncertaintyDistFlag.")
end



#Set flag variable decalring whether GP visits should be obtained from
#ascertainment probability
GPVisitFromAscertainProbFlag = 2
#1 for yes (ascertained cases correspond to GP cases); 0 for no (instead computed from febrile prob and propensity to consult GP params)
#2 for Ascertained cases correspond to GP cases with risk group specific probability

#--------------------------------------------------------------------------
### Vaccine uptake data
#--------------------------------------------------------------------------
VaccScheduleSymb = Symbol(ARGS[10]) #Convert string to Symbol
VaccScheduleFn  = getfield(Main, VaccScheduleSymb) #Make Symbol callable functions

#--------------------------------------------------------------------------
### Vaccine administration cost (fixed reference price)
#--------------------------------------------------------------------------
AdminChargePerVaccDose = 10.

#--------------------------------------------------------------------------
### Establish uncertainty distribution options to generate health economic parameters from
### desired uncertainty distributions
#--------------------------------------------------------------------------

#Note, value of 0 means use point estimate only!
if UncertaintyDistFlag == 1 #Use uncertainty distributions
    UncertaintyDistFlag_GPAscertain = 2 #1: Triangular distribution;
                                      #2: Febrile symptomns from triangular distribution. Load GP ascertainment rates from file.
                                      #3: Febrile symptomns pointwise. Load GP ascertainment rates from file.
    UncertaintyDistFlag_HealthEpsRelLhood = 2 #1: Normal; #2 Bootstrap
    UncertaintyDistFlag_HealthEpsCosts = 1 #1: Lognormal; #2: Gamma
    UncertaintyDistFlag_QALYloss = 1 #1:Normal;
    UncertaintyDistFlag_VaccCosts = 0 #1:Uniform integer; #2: Triangular distribution
elseif UncertaintyDistFlag == 2 #Use point estimates
    UncertaintyDistFlag_GPAscertain = 3 #1: Triangular distribution;
                                    #2: Febrile symptomns from triangular distribution. Load GP ascertainment rates from file.
                                    #3: Febrile symptomns pointwise. Load GP ascertainment rates from file.
    UncertaintyDistFlag_HealthEpsRelLhood = 0 #1: Normal; #2 Bootstrap
    UncertaintyDistFlag_HealthEpsCosts = 0 #1: Lognormal; #2: Gamma
    UncertaintyDistFlag_QALYloss = 0 #1:Normal;
    UncertaintyDistFlag_VaccCosts = 0 #1:Uniform integer; #2: Triangular distribution
else
    error("UncertaintyDistFlag has value $UncertaintyDistFlag. Must take value 1 or 2.")
end

#--------------------------------------------------------------------------
### Carry out health economic evaluation for each scheme
#--------------------------------------------------------------------------

#Loop over each job
#Save outputs to distinct MAT file
@sync @distributed for JobID = 1:JobNum

    #--------------------------------------------------------------------------
    ### Reseed RNG.
    #--------------------------------------------------------------------------
    Random.seed!(RNG_SeedNum)
    #--------------------------------------------------------------------------

    #Call vaccine uptake function. Uses JobID to populate age groups distinctly for each job.
    VaccUptake= VaccScheduleFn(JobID,AgeGrpNum,RiskGroupFlag)

    #Use JOBID to access correct part of file

    #--------------------------------------------------------------------------
    ### Load relevant epi data files (for risk group specific and population
    ### level runs)
    #--------------------------------------------------------------------------
    LowRiskInfDataFile = "LowRiskFileName.mat"
    AtRiskInfDataFile = "AtRiskFileName.mat"
    AllRiskInfDataFile = "AllRiskFileName.mat"

    #--------------------------------------------------------------------------
    #Call function to correct discrpency by finding difference between risk group sum and overall estimates,
    #Aportion difference weighted by the relative contribution to total from each risk group
    #--------------------------------------------------------------------------
    PropnLowRiskInfectedCorrected,PropnAtRiskInfectedCorrected,
        PropnLowRiskAscertainCorrected,PropnAtRiskAscertainCorrected =
        CorrectRiskGrpInfected_MultiJobVers(LowRiskInfDataFile,
                                 AtRiskInfDataFile,
                                 AllRiskInfDataFile,
                                 RiskGrpPropnData)
                                 #JobID)

    #--------------------------------------------------------------------------
    #Call function to comptue risk group specific ascertainment rates
    #LowRiskAscertainRates & AtRiskAscertainRates are tuples, with entry per
    #evaluation run. Each entry contains a SeasonNum*AgeGrpNum array.
    #--------------------------------------------------------------------------
    LowRiskAscertainRates,AtRiskAscertainRates = RiskGrpAscertainCalc_MultiJobVers(PropnLowRiskInfectedCorrected,
                                                                        PropnAtRiskInfectedCorrected,
                                                                        AllRiskInfDataFile,
                                                                        PopnDataCell,
                                                                        EpiParamIdx,AtRiskScaleFactor)
                                                                        #JobID)

    #Aggregate group specific ascertainment rates into a cell array
    GrpSpecificAscertainRates = [LowRiskAscertainRates,AtRiskAscertainRates]

    #--------------------------------------------------------------------------
    ### Construct PropnInfData input, each entry a single model simulation output
    #--------------------------------------------------------------------------

    #Initialise tuples
    PropnLowRiskInfectedCorrected_EvalRunList = Array{Array{Float64,3},1}(undef,TotalParamSetsEval)
    PropnLowRiskAscertainCorrected_EvalRunList = Array{Array{Float64,3},1}(undef,TotalParamSetsEval)
    PropnAtRiskInfectedCorrected_EvalRunList = Array{Array{Float64,3},1}(undef,TotalParamSetsEval)
    PropnAtRiskAscertainCorrected_EvalRunList = Array{Array{Float64,3},1}(undef,TotalParamSetsEval)

    #For each evaluation, assign infection & ascertained cases based on EpiParamIdx
    for EvalIdx = 1:TotalParamSetsEval

        #Get entry to be accessed from infection&ascertainment datasets
        CurrentEval_EpiParamIdx = EpiParamIdx[EvalIdx]

        #Assign selected data to EvalRunList tuples
        PropnLowRiskInfectedCorrected_EvalRunList[EvalIdx] = PropnLowRiskInfectedCorrected[CurrentEval_EpiParamIdx]
        PropnLowRiskAscertainCorrected_EvalRunList[EvalIdx] = PropnLowRiskAscertainCorrected[CurrentEval_EpiParamIdx]

        PropnAtRiskInfectedCorrected_EvalRunList[EvalIdx] = PropnAtRiskInfectedCorrected[CurrentEval_EpiParamIdx]
        PropnAtRiskAscertainCorrected_EvalRunList[EvalIdx] = PropnAtRiskAscertainCorrected[CurrentEval_EpiParamIdx]
    end

    if RiskGroupFlag == 0
        PropnInfData = PropnLowRiskInfectedCorrected_EvalRunList
        PropnAscertainedData = PropnLowRiskAscertainCorrected_EvalRunList
    elseif RiskGroupFlag == 1
        PropnInfData = PropnAtRiskInfectedCorrected_EvalRunList
        PropnAscertainedData = PropnAtRiskAscertainCorrected_EvalRunList
    end

    #--------------------------------------------------------------------------
    ### Pass options & input data into RunHealthEconEval.
    ### Perform multiple evaluations!
    #--------------------------------------------------------------------------
    ByYrOfAgeHealthEpsCounts,ByYrOfAgeHealthEpsCosts,ByYrOfAgeHealthEpsQALYloss,
                AllAgeHealthEpsCounts,AllAgeHealthEpsCosts,AllAgeHealthEpsQALYloss,
                VaccCosts,VaccDeployedTuple,VaccineDoseCosts =
    RunHealthEconEval(RiskGroupFlag,AgeStrucModelFlag,AgeGrpNum,StrainNum,TotalParamSetsEval,
        DiscountRate,VaccUptake,AdminChargePerVaccDose,
        PropnInfData,PropnAscertainedData,GPVisitFromAscertainProbFlag,GrpSpecificAscertainRates,
        UncertaintyDistFlag_GPAscertain,UncertaintyDistFlag_HealthEpsRelLhood,
        UncertaintyDistFlag_HealthEpsCosts,UncertaintyDistFlag_QALYloss,
        UncertaintyDistFlag_VaccCosts)

    # #--------------------------------------------------------------------------
    # ### Save outputs to file
    # #--------------------------------------------------------------------------
    # SaveID = ARGS[11]
    # if RiskGroupFlag == 0
    #     OutputFName = ".mat"
    # elseif RiskGroupFlag == 1
    #     OutputFName = ".mat"
    # end
    #
    # file = matopen(OutputFName, "w")
    # write(file, "ByYrOfAgeHealthEpsCounts", ByYrOfAgeHealthEpsCounts)
    # write(file, "ByYrOfAgeHealthEpsCosts", ByYrOfAgeHealthEpsCosts)
    # write(file, "ByYrOfAgeHealthEpsQALYloss", ByYrOfAgeHealthEpsQALYloss)
    # write(file, "AllAgeHealthEpsCounts", AllAgeHealthEpsCounts)
    # write(file, "AllAgeHealthEpsCosts", AllAgeHealthEpsCosts)
    # write(file, "AllAgeHealthEpsQALYloss", AllAgeHealthEpsQALYloss)
    # write(file, "VaccCosts", VaccCosts)
    # write(file, "VaccDeployedTuple", VaccDeployedTuple)
    # write(file, "VaccineDoseCosts", VaccineDoseCosts)
    # close(file)

end
### MAIN SCRIPT ENDS ###

end

#--------------------------------------------------------------------------
# PASS TO FUNCTION
#--------------------------------------------------------------------------
RunHealthEcon_VaccAgeSweep(ARGS)

for ii in workers()
    rmprocs(ii)
end
