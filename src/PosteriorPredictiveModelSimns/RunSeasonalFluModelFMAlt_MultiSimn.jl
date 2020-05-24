#Purpose:
#Script to multiple model simulations with parameters read from file.

#Model specifics:
#SEIR influenza deterministic transmission dynamic model ODEs
#Includes age & strain structure.
#Exposure history (for previous season only) included
#--------------------------------------------------------------------------

# Load libraries
using Combinatorics
using DifferentialEquations
using DataFrames
using XLSX
using DelimitedFiles
using LinearAlgebra
using MAT



#-------------------------------------------------------------------------------
###  VACCINATION UPTAKE FUNCTIONS
#-------------------------------------------------------------------------------

#Based on historical data
function HistoricalVaccUptake(SeasonsToSimulate)

    #BASED ON HISTORICAL DATA

    #Cell per season
    #Per cell, 2D array
    #rows for age (0 to 90+), cols for calendar day of year
    VaccUptakeBySeason = Vector{Array}(undef,SeasonsToSimulate)

    # Import the data - Pandemic flu vacc (2009/2010 season)
    PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc.xlsx","Both","C3:NC103")

    # Collate into Array, Assign to storage cell
    PandemicFluVaccUptakeArray = Array{Float64, 2}(PandemicFluVaccUptake)
    #PandemicFluVaccUptakeArray[66:end,:] .= 0 #Set specified age groups to have no vaccination
    VaccUptakeBySeason[1] = PandemicFluVaccUptakeArray

    #Import the data - Sesonal flu vacc
    SheetNames = ["2010_2011","2011_2012","2012_2013","2013_2014",
        "2014_2015","2015_2016","2016_2017","2017_2018"]

    for ii = 1:length(SheetNames)
        HistoricalSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeBySeasonCalYr_All.xlsx","$(SheetNames[ii])","C3:NC103")

        # Collate into Array, Assign to storage cell
        VaccUptakeBySeasonArray = Array{Float64, 2}(HistoricalSeasonalFluVaccUptake)
        VaccUptakeBySeason[ii+1] = VaccUptakeBySeasonArray  #Add 1 to idx as first entry is for 2009/2010 season
    end

return VaccUptakeBySeason
end

#Based on historical data, at-risk vaccine uptake only
function HistoricalAtRiskVaccOnlyUptake(SeasonsToSimulate)

    #BASED ON HISTORICAL DATA

    #Cell per season
    #Per cell, 2D array
    #rows for age (0 to 90+), cols for calendar day of year
    VaccUptakeBySeason = Vector{Array}(undef,SeasonsToSimulate)

    # Import the data - Pandemic flu vacc (2009/2010 season)
    PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc.xlsx","Both","C3:NC103")

    # Collate into Array, Assign to storage cell
    PandemicFluVaccUptakeArray = Array{Float64, 2}(PandemicFluVaccUptake)
    #PandemicFluVaccUptakeArray[66:end,:] .= 0 #Set specified age groups to have no vaccination
    VaccUptakeBySeason[1] = PandemicFluVaccUptakeArray

    #Import the data - Sesonal flu vacc
    SheetNames = ["2010_2011","2011_2012","2012_2013","2013_2014",
        "2014_2015","2015_2016","2016_2017","2017_2018"]

    for ii = 1:length(SheetNames)
        HistoricalSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_All_HistoricalAtRiskOnly.xlsx","$(SheetNames[ii])","C3:NC103")

        # Collate into Array, Assign to storage cell
        VaccUptakeBySeasonArray = Array{Float64, 2}(HistoricalSeasonalFluVaccUptake)
        VaccUptakeBySeason[ii+1] = VaccUptakeBySeasonArray  #Add 1 to idx as first entry is for 2009/2010 season
    end

return VaccUptakeBySeason
end


function RunSeasonalFluModelFMAlt_MultiSimns(AscertainProbFn,ExpHistVaccType,FitAggAgeTuple,AgeGrpSuscepTuple,
                                            ParticleSets,SeedRunID,TotalRunNum,OutputFName,StoreFlag_PopnFOI,
                                            VaccUptakeFunction)
#Inputs:
#   AscertainProbFn  - Description given before fn call
#   ExpHistVaccType - See description given before fn call
#   FitAggAgeTuple - See description given before fn call
#   AgeGrpSuscepTuple - See description given before fn call
#   ParticleSets - Parameter sets to be simualted
#   SeedRunID - Index value. TO be incremented by run number giving unique
#                  save filename
#   TotalRunNum - Number of simulations, with unique parameter sets, to be
#                   performed
#   OutputFName - Save location
#   StoreFlag_PopnFOI - (Flag variable) If active, record the overall force of infection at each timestep
#   VaccUptakeFunction - (function) Module to load required vaccine uptake data

#-------------------------------------------------------------------------------
# SPECIFY TYPE OF RUN THROUGH FLAG VARIABLE
#-------------------------------------------------------------------------------
# (INFLUENCES VACCINE UPTAKE/EFFICACY, & USE OF ODEBurnIn, ODEH1N1OnlyTime, ODEAlleStrainTime)
#1 - exploratory; 2 - historical; 3 - alternative vacc. scheme
SimnRunType = 2

#------------------------------------------------------------------------------
###  LOAD CONTACT DATA
ContactArray = readdlm("../../Data/ContactData/UKAdeqContact_Over100Vers.csv",',')

#------------------------------------------------------------------------------
### IMPORT MORTALITY RATES
#MortalityFile = "C:/Users/edwar/Documents/GitHub/DoH-FluVaccine/Data/DemographicData/MortalityRateByAge_ONS90plus.txt"
MortalityFile = "../../Data/DemographicData/ModifiedData/MortalityProbPerAge0to100_EH.txt"

#------------------------------------------------------------------------------
### IMPORT INITIAL AGE DISTRIBUTION

#Set proportion of individuals initially in each age
ONSPopnDistEngland20102018 = readdlm("../../Data/DemographicData/ONSPopnDistEngland20102018_0to100.txt",',')

#------------------------------------------------------------------------------
### DISEASE DYNAMICS SIMN/FLAG PARAMETERS
SimnStartDate=9 #Month of year to start simulation on

#Run time for ODE model. Take values post burn in
if SimnRunType == 1
    ODEBurnIn = 0*365
    ODEStaticPopnTime = 20*365
    ODEInferenceTime = 7*365
    ODEForwardSimnTime = 4*365
elseif SimnRunType == 2

    SeasonsToSimulate = 9

    ODEBurnIn = 0*365
    ODEStaticPopnTime = 1*365
    ODEInferenceTime = (SeasonsToSimulate-1)*365
    ODEForwardSimnTime = 0*365
elseif SimnRunType == 3

    #Specify number of seasons that will use historical data
    HistoricalSeasonNum = 7

    ODEBurnIn = 0*365
    ODEStaticPopnTime = 20*365
    ODEInferenceTime = HistoricalSeasonNum*365
    ODEForwardSimnTime = 4*365
else
    error("Incorrect RunType entered")
end

MaxTime = ODEBurnIn + ODEStaticPopnTime + ODEInferenceTime + ODEForwardSimnTime

#Compute total time post burn-in
ODESampleTime = ODEStaticPopnTime + ODEInferenceTime + ODEForwardSimnTime

#Specify timestep to use in ode45 scheme
timestep = 1.

#Concatenate simulation variables
SimnParam=[SimnStartDate,ODEBurnIn,ODEStaticPopnTime,ODEInferenceTime,ODEForwardSimnTime,
            MaxTime,ODESampleTime,timestep,StoreFlag_PopnFOI]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
###  TRANSMISSION RELATED PARAMETERS
NumOfStrains = 4 #Specify number of strains in the system (A/H1N1, A/H3N2, two B lineages)
gamma = 1/3.8*ones(NumOfStrains) #recovery rate
sigma = [1/1.4,1/1.4,1/0.6,1/0.6] #latent rate

#Calculate beta, group infectious status parameters into vector
R_0 = [1.65,1.45,1.5,1.6]
spec_rad = maximum(abs.(eigvals(ContactArray)::Array{Float64,1})) #the spectral radius of C is the maximum of the absolute values of the eigenvalues
beta = gamma.*R_0/spec_rad
InfectionParam = [beta,sigma,gamma,R_0]

#Set proportion of population that are susceptible at end of season within
#a natural infection exposure history class to remain in that class
#(rather than move to the Naive exposure history class)
MultiSeasonImmPropn = 0

#--------------------------------------------------------------------------
### SUSCEPTIBILITY RELATED PARAMETERS
#--------------------------------------------------------------------------
M = size(ContactArray,1)

#Number of exposure history classes
#One for naive, one for vacc with no natural infection, one per strain
#(unvacc), one per strain with vacc
ExpHistNum = (NumOfStrains*2) + 2

#-------------------------------------------------------------------------------
###  VACCINATION UPTAKE
#-------------------------------------------------------------------------------
VaccUptakeBySeason = VaccUptakeFunction(SeasonsToSimulate)

#--------------------------------------------------------------------------
### VACCINATION - LEAKY TRANSMISSION SETTINGS
#--------------------------------------------------------------------------

#Set flag variable for "leaky" transmission being unactive or active
#0 - Infectiousness of infected vacc. group unmodified.
#1 - Infected vacc. group has reduced infectiousness
LeakyTransFlag = 0

#Validity check on LeakyTransFlag value
if LeakyTransFlag !=0 && LeakyTransFlag !=1
    error("Invalid value of $LeakyTransFlag for LeakyTransFlag, should be 0 or 1")
end

#--------------------------------------------------------------------------
### VACCINATION - EFFICACY
#--------------------------------------------------------------------------

#Set leaky vaccine efficacy parameters
#Row ii - Age group ii
#Columns for strains (A/H1N1, A/H3N2, two B lineages)
if SimnRunType == 1
    #SYNTHETIC DATA

    if LeakyTransFlag == 0
        LeakyVaccVarBySeason = ones(M,NumOfStrains)
    elseif LeakyTransFlag == 1
        alpha = ones(M,NumOfStrains)
        delta = zeros(M,NumOfStrains)
        LeakyVaccVarBySeason = [alpha,delta]
    end
elseif SimnRunType == 2
    #BASED ON HISTORICAL DATA

    #Cell per season
    #Per cell, 2D array
    #rows for age (0 to 90+), cols for strain
    VaccEfficacy_ModelFMAlt = Vector{Array}(undef,SeasonsToSimulate)

    #Assign data for 2009/2010 season (assuming pandemic flu vaccine)
    # 72% for all ages!
    VaccEfficacy_ModelFMAlt[1] = 0.72*ones(M,NumOfStrains)
    VaccEfficacy_ModelFMAlt[1][:,2:4] .= 0

    # Import the data - Sesonal flu vacc
    SheetNames = ["2010_2011","2011_2012","2012_2013","2013_2014",
        "2014_2015","2015_2016","2016_2017","2017_2018"]

    for ii = 1:length(SheetNames)
        VaccEffTemp = XLSX.readdata("../../Data/VaccEfficacy/VaccEfficacy_ByYrOfAge.xlsx","$(SheetNames[ii])","C3:F103")

        # Collate into Array, Assign to storage cell
        VaccEfficacy_ModelFMAlt[ii+1] = Array{Float64, 2}(VaccEffTemp)
    end

    #Assign to LeakyVaccVar variable
    #Cell for each season
    #Row for each age, column for each strain
    if LeakyTransFlag == 0
        LeakyVaccVarBySeason = VaccEfficacy_ModelFMAlt
    elseif LeakyTransFlag == 1
        alpha = VaccEfficacy_ModelFMAlt
        delta = VaccEfficacy_ModelFMAlt

        #delta = Vector(HistoricalSeasonNum+1)
        #for ii = 1:length(delta)
        #    delta[ii] = ones(M,NumOfStrains)
        #end
        LeakyVaccVarBySeason = [alpha,delta]
    end
elseif SimnRunType == 3
    #RUN ALTERNATIVE VACC. SCHEMES

    #Cell per season
    #Per cell, 2D array
    #rows for age (0 to 90+), cols for strain
    VaccEfficacy_ModelFMAlt = Vector{Array}(undef,SeasonsToSimulate)

    #Assign data for 2009/2010 season (assuming pandemic flu vaccine)
    # 72% for all ages!
    VaccEfficacy_ModelFMAlt[1] = 0.72*ones(M,NumOfStrains)
    VaccEfficacy_ModelFMAlt[1][:,2:4] = 0

    # Import the data - Sesonal flu vacc
    SheetNames = ["2010_2011","2011_2012","2012_2013","2013_2014",
        "2014_2015","2015_2016","2016_2017","2017_2018"]

    for ii = 1:length(SheetNames)
        VaccEffTemp = XLSX.readdata("./../Data/VaccEfficacy/VaccEfficacy_ByYrOfAge.xlsx","$(SheetNames[ii])","C3:F103")

        # Collate into Array, Assign to storage cell
        VaccEfficacy_ModelFMAlt[ii+1] = Array{Float64, 2}(VaccEffTemp)
    end

    #Assign to LeakyVaccVar variable
    #Cell for each season
    #Row for each age, column for each strain
    if LeakyTransFlag == 0
        LeakyVaccVarBySeason = VaccEfficacy_ModelFMAlt
    elseif LeakyTransFlag == 1
        alpha = VaccEfficacy_ModelFMAlt
        delta = VaccEfficacy_ModelFMAlt

        #delta = Vector(HistoricalSeasonNum+1)
        #for ii = 1:length(delta)
        #    delta[ii] = ones(M,NumOfStrains)
        #end
        LeakyVaccVarBySeason = [alpha,delta]
    end
else
    error("Incorrect RunType entered")
end


#--------------------------------------------------------------------------
### INITIAL INF. PROPORTION
#--------------------------------------------------------------------------
#Column i for strain i
#Row 1 for when only H1N1 infection allowed
#Row 2 when infection by any strain allowed
InfPropn_StartOfSeason = [1e-5 0 0 0;
                            2.5e-6 2.5e-6 2.5e-6 2.5e-6]

#--------------------------------------------------------------------------
### GIVE DETAILS ON LOADING INITIAL SUSCEPTIBLE CLASS CONDITIONS FROM FILE IF NEEDED
#--------------------------------------------------------------------------
ICFromFile = [[0],""]

#--------------------------------------------------------------------------
### NUMBER OF INFLUENZA SEASONS WORTH OF DATA BEING CONSIDERED
#--------------------------------------------------------------------------
RetainedSeasonNum = 6

#--------------------------------------------------------------------------
# AGGREGATE FIXED PARAMETERS
#--------------------------------------------------------------------------
FixedModelParams = [ContactArray,MortalityFile,ONSPopnDistEngland20102018,
                    SimnRunType,ExpHistVaccType,StoreFlag_PopnFOI,
                    SimnParam,NumOfStrains,ExpHistNum,M,InfectionParam,
                    MultiSeasonImmPropn,VaccUptakeBySeason,LeakyTransFlag,
                    LeakyVaccVarBySeason,InfPropn_StartOfSeason,ICFromFile,
                    RetainedSeasonNum,AscertainProbFn,
                    FitAggAgeTuple,AgeGrpSuscepTuple]

#--------------------------------------------------------------------------
### INITIALISE CELL TO STORE OUTPUTS OF INTEREST
#--------------------------------------------------------------------------
if StoreFlag_PopnFOI == 0
    OutputSimnCell = Array{Array{Float64},2}(undef,TotalRunNum,3)
elseif StoreFlag_PopnFOI == 1
    OutputSimnCell = Array{Array{Float64},2}(undef,TotalRunNum,4)
end

for ii = 1:TotalRunNum

    #----------------------------------------------------------------------
    ### Read in parameter values from file
    #----------------------------------------------------------------------
    ParticleSetVal = ParticleSets[ii,:]

    #----------------------------------------------------------------------
    ### CALL FUNCTION TO RUN MODEL
    #----------------------------------------------------------------------
    OutputSimnCell[ii,:] = RunModelFMAlt(FixedModelParams,ParticleSetVal)

end

#----------------------------------------------------------------------
### SAVE TO OUTPUT VARIABLE
#----------------------------------------------------------------------
OutputFile = matopen(OutputFName,"w")
write(OutputFile, "SimnData", OutputSimnCell)
close(OutputFile)


end
#--------------------------------------------------------------------------
### ADD FILES TO SEARCH PATH FOR ODES/MODEL RUN FUNCTION
include("../ModelFns/RunSeasonalFluModelODEs.jl")
include("../ModelFns/ExpHistUpdateJulV1.jl")
include("../ModelFns/ModelExtraFns.jl")

#--------------------------------------------------------------------------
# Declare the type of ascertainment probabiity modification to be formed
#--------------------------------------------------------------------------
AscertainProbFn = SeasonWithAgeScaleLinearAscertainmentFn
#AscertainProbFn = SeasonWithAgeScaleLinearAscertainmentFn_TypeSpecific

#--------------------------------------------------------------------------
#Select exposure history susceptibility modification form for vaccine-related states
#0 - Fixed values every season
#1 - Relate to previous season vaccine efficacy
#2 - Relate to previous season vaccine efficacy AND age band
#--------------------------------------------------------------------------
ExpHistVaccType = 1

if ExpHistVaccType != 0 && ExpHistVaccType !=1 && ExpHistVaccType !=2
    error("ExpHistVaccType must take value 0, 1 or 2. Current value is $(ExpHistVaccType)")
end

#--------------------------------------------------------------------------
# Declare bounds for age category groupings (for aggregating case counts)
#--------------------------------------------------------------------------
FitAggAgeFlag = 0 #Indicator variable.
# Set to: -> 1 to fit to aggregetated age band counts
#         -> 0 to fit to single yr age classes

#Error check
if FitAggAgeFlag != 0 && FitAggAgeFlag !=1
   error("FitAggAgeFlag must take value 0 or 1. Current value is $FitAggAgeFlag")
end

AgeBandLowerBounds = [0 2 18 65 85]
AgeBandUpperBounds = [1 17 64 84 100]
AgeBandNum = length(AgeBandLowerBounds)

#Concantenate bounds into a single array
AgeBandBounds = [AgeBandLowerBounds;AgeBandUpperBounds]
                #vertically stack age category lower bounds above upper bounds

#Place FitAggAgeFlag&AgeBandBoundsi nto a cell/tuple
FitAggAgeTuple = [FitAggAgeFlag,AgeBandBounds]

#--------------------------------------------------------------------------
# Declare bounds for age susceptibility groupings
#--------------------------------------------------------------------------
# AgeSuscepLowerBounds = [0,18,65]
# AgeSuscepUpperBounds = [17,64,100]
AgeSuscepLowerBounds = [0,18,65,85]
AgeSuscepUpperBounds = [17,64,84,100]
AgeGrpSuscepParamNum = length(AgeSuscepLowerBounds)

#Select age&strain susceptibility modification form
# AgeAndStrainDepSuscepFn - Age group & strain specific
# AgeDepWithStrainScalingSuscepFn - Age group specific, with strain modifier
# StrainDepWithAgeScalingSuscepFn - Strain specific, with age group modifier
# AgeOnlySuscepFn - Age group speicifc only (independent of strain)
AgeSuscepType = AgeOnlySuscepFn

#Concantenate into a cell tuple
AgeGrpSuscepTuple = [AgeGrpSuscepParamNum,AgeSuscepLowerBounds,AgeSuscepUpperBounds,AgeSuscepType]

#--------------------------------------------------------------------------
# Read in parameter values from file
#--------------------------------------------------------------------------
ParamInputFile = "../../Results/FMAlt_OptimFitOutput/FMAlt_OptimFitOutputFiles/OptimiserParamsTrace#14_Emp_TransContactMatrix_JobArrayCombined.txt";
# ParticleSets = readdlm(ParamInputFile,'\t')
ParticleSets = readdlm(ParamInputFile,',')

#----------------------------------------------------------------------
### SPECIFY OUTPUT FILE SAVE LOCATION
#----------------------------------------------------------------------
SeedRunID = "TransContactArrayTest_AtRiskStrat"
OutputFName = ".../../Results/FMAlt_SimnJuliaV1ModelRun_#$(SeedRunID).mat"

#----------------------------------------------------------------------
### SPECIFY POPULATION-LEVEL FORCE OF INFECTION FLAG OPTION
#----------------------------------------------------------------------
#0 - inactive, 1 - active
StoreFlag_PopnFOI = 1

#----------------------------------------------------------------------
### SPECIFY VACCINE UPTAKE FUNCTION NAME
#----------------------------------------------------------------------
VaccUptakeFunction = HistoricalAtRiskVaccOnlyUptake

#--------------------------------------------------------------------------
#Call function to run multiple simulations
#--------------------------------------------------------------------------
TotalRunNum = 10
RunSeasonalFluModelFMAlt_MultiSimns(AscertainProbFn,ExpHistVaccType,FitAggAgeTuple,AgeGrpSuscepTuple,
                                    ParticleSets,SeedRunID,TotalRunNum,OutputFName,StoreFlag_PopnFOI,
                                    VaccUptakeFunction)
