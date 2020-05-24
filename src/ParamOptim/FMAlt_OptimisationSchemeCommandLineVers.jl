#Purpose:
#Run optimisation scheme on the transmission model.
#Attempt to find minima using gradient free techniques
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# LOAD REQUIRED ENVIRONMENT (FILE PATH BASED ON CWD BEING THE LOCATION THIS FILE RESIDES)
#--------------------------------------------------------------------------
using Pkg
Pkg.activate("../")

#--------------------------------------------------------------------------
# ADD FILES TO SEARCH PATH FOR ODES/MODEL RUN FUNCTION
#--------------------------------------------------------------------------
include("../ModelFns/RunSeasonalFluModelODEs.jl")
include("../ModelFns/ExpHistUpdateJulV1.jl")
include("../ModelFns/ModelExtraFns.jl")

#-------------------------------------------------------------------------------
# LOAD REQUIRED PACKAGES
#-------------------------------------------------------------------------------
using Distributions
using StatsBase
using XLSX
using DifferentialEquations
using Random
using DelimitedFiles
using LinearAlgebra
using Optim

#-------------------------------------------------------------------------------
# DEFINE AND GROUP MODEL SIMULATION FIXED PARAMETERS
#-------------------------------------------------------------------------------
function FixedModelParamGenFn(SeasonsToSimulate,AscertainProbFn,ExpHistVaccType,
								FitAggAgeTuple,AgeGrpSuscepTuple,ObvsData)
#Inputs:
#	SeasonsToSimulate - (integer)
#	AscertainProbFn - (function name)
#	ExpHistVaccType - (integer, flag variable)
#	FitAggAgeTuple - (tuple)  [FitAggAgeFlag,AgeBandBounds]
#	AgeGrpSuscepTuple - (tuple) [AgeGrpSuscepParamNum,AgeSuscepLowerBounds,AgeSuscepUpperBounds,AgeSuscepType]
#	ObvsData - (array)

#Outputs:
#   FixedModelParams - (tuple) Variables used in model simulation

	#-------------------------------------------------------------------------------
	# SPECIFY TYPE OF RUN THROUGH FLAG VARIABLE
	#-------------------------------------------------------------------------------
	# (INFLUENCES VACCINE UPTAKE/EFFICACY, & USE OF ODEBurnIn, ODEH1N1OnlyTime, ODEAlleStrainTime)
	#1 - exploratory; 2 - historical; 3 - alternative vacc. scheme
	SimnRunType = 2

	#------------------------------------------------------------------------------
	###  LOAD CONTACT DATA
	ContactArray = readdlm("../../Data/DemographicData/UKAdeqContact_Over100Vers.csv",',')

	#------------------------------------------------------------------------------
	### IMPORT MORTALITY RATES
	MortalityFile = "../../Data/DemographicData/MortalityProbPerAge0to100_EH.txt"

	#------------------------------------------------------------------------------
	### IMPORT INITIAL AGE DISTRIBUTION

	#Set proportion of individuals initially in each age
	ONSPopnDistEngland20102018 = readdlm("../../Data/DemographicData/ONSPopnDistEngland20102018_0to100.txt",',')

	#------------------------------------------------------------------------------
	### DISEASE DYNAMICS SIMN/FLAG PARAMETERS
	SimnStartDate=9 #Month of year to start simulation on

	#Store population-level FOI flag option
	#0 - inactive, 1 - active
	StoreFlag_PopnFOI = 0

	#Run time for ODE model. Take values post burn in
	if SimnRunType == 1
	    ODEBurnIn = 0*365
	    ODEStaticPopnTime = 20*365
	    ODEInferenceTime = 7*365
	    ODEForwardSimnTime = 4*365
	elseif SimnRunType == 2
	    ODEBurnIn = 0*365
	    ODEStaticPopnTime = 1*365
		ODEInferenceTime = (SeasonsToSimulate-1)*365;
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
	timestep = 1

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
	R_0 = [1.8,1.2,1.5,1.6]
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

	#------------------------------------------------------------------------------
	###  VACCINATION UPTAKE
	#--------------------------------------------------------------------------
	if SimnRunType == 1
	    #SYNTHETIC DATA

	    #Set vaccine uptake rate
	    #Entry (ii,kk)- Age class ii, day kk
	    NumVaccDays = 91
	    TargetCovPerSeason = ones(M,NumVaccDays)*((91/365)/3)
	    VaccUptakeBySeason = zeros(M,365) #Proportion given the vaccine during each day
	    VaccUptakeBySeason[:,244:334] .= TargetCovPerSeason/NumVaccDays

	elseif SimnRunType == 2
	    #BASED ON HISTORICAL DATA

	    #Cell per season
	    #Per cell, 2D array
	    #rows for age (0 to 100+), cols for calendar day of year
		VaccUptakeBySeason = Vector{Array}(undef,SeasonsToSimulate)

	    # Import the data - Pandemic flu vacc (2009/2010 season)
	    PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc.xlsx","Both","C3:NC103")

	    # Collate into Array, Assign to storage cell
	    VaccUptakeBySeason[1] = Array{Float64, 2}(PandemicFluVaccUptake)

	    #Import the data - Sesonal flu vacc
	    SheetNames = ["2010_2011","2011_2012","2012_2013","2013_2014",
	        "2014_2015","2015_2016","2016_2017","2017_2018"]

	    for ii = 1:length(SheetNames)
	        HistoricalSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_DailyVaccUptakeBySeasonCalYr_All.xlsx","$(SheetNames[ii])","C3:NC103")

	        # Collate into Array, Assign to storage cell
	        VaccUptakeBySeason[ii+1] = Array{Float64, 2}(HistoricalSeasonalFluVaccUptake)
	            #Add 1 to idx as first entry is for 2009/2010 season
	    end

	elseif SimnRunType == 3
	    #RUN ALTERNATIVE VACC. SCHEMES

	    #Cell per season
	    #Per cell, 2D array
	    #rows for age (0 to 100+), cols for calendar day of year
		VaccUptakeBySeason = Vector{Array}(undef,SeasonsToSimulate)

	    # Import the data - Pandemic flu vacc (2009/2010 season)
	    PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc.xlsx","Both","C3:NC103")

	    # Collate into Array, Assign to storage cell
	    VaccUptakeBySeason[1] = Array{Float64, 2}(PandemicFluVaccUptake)

	    #Import the data - Sesonal flu vacc
	    SheetNames = ["2010_2011","2011_2012","2012_2013","2013_2014",
	        "2014_2015","2015_2016","2016_2017","2017_2018"]

	    for ii = 1:length(SheetNames)
	        HistoricalSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_DailyVaccUptakeBySeasonCalYr_All.xlsx","$(SheetNames[ii])","C3:NC103")

	        # Collate into Array, Assign to storage cell
	        VaccUptakeBySeason[ii+1] = Array{Float64, 2}(HistoricalSeasonalFluVaccUptake)
	            #Add 1 to idx as first entry is for 2009/2010 season
	    end
	else
	    error("Incorrect RunType entered")
	end

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
	    #rows for age (0 to 100+), cols for strain
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
	    #rows for age (0 to 100+), cols for strain
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
	RetainedSeasonNum = size(ObvsData,1)

	#--------------------------------------------------------------------------
 	### AGGREGATE FIXED PARAMETERS
	#--------------------------------------------------------------------------
	FixedModelParams = [ContactArray,MortalityFile,ONSPopnDistEngland20102018,
	                    SimnRunType,ExpHistVaccType,StoreFlag_PopnFOI,
	                    SimnParam,NumOfStrains,ExpHistNum,M,InfectionParam,
	                    MultiSeasonImmPropn,VaccUptakeBySeason,LeakyTransFlag,
	                    LeakyVaccVarBySeason,InfPropn_StartOfSeason,ICFromFile,
	                    RetainedSeasonNum,AscertainProbFn,
	                    FitAggAgeTuple,AgeGrpSuscepTuple]

	return FixedModelParams
end


#-------------------------------------------------------------------------------
# SUMMARY STATISTIC FUNCTION DEFINITION (with parameter sets failing the temporal check assigned a finite, large error value)
#-------------------------------------------------------------------------------
function  FMAlt_SummStatFun_PoissDevWithTemporalMod(ObservedData,SimnData)

#Inputs:
#   ObservedData,SimnData - (arrays)

#Outputs:
#   AmendedSummStatVal - (float) Error value

    #Disaggregate simulation data
    x = SimnData[1]
    Total_I = SimnData[2]

	#Compute poisson deviance
	SummStatIterNum = length(ObservedData)
	TempSum = 0.0
    for ii = 1:SummStatIterNum
        if ObservedData[ii].!=0
            TempSum = TempSum + (ObservedData[ii]*log(ObservedData[ii]/x[ii])) - (ObservedData[ii] - x[ii])
        end
    end
    SummStatVal_Overall = 2*(TempSum + sum(x[ObservedData.==0]))

    #Perform temporal check
    NumOfSeasons = convert(Int64,size(Total_I,1)/366) #%Number of seasons obtained by dividing number of daily records by days in yr (+1 to account for day 0 recording1)
	AgeGrpNum = convert(Int64,size(Total_I,2)) #Number of age groups, number of columns of Total_I
	Total_I_Array = zeros(NumOfSeasons,AgeGrpNum,366)
	Total_I_Transpose = Total_I'
    for jj = 1:NumOfSeasons
       StartIdx = ((jj-1)*366) + 1
       EndIdx = jj*366

	   #Allocate seasons worth of records to new array, slice per season
       #Take transpose of Total_I so array dimensions agree
       Total_I_Array[jj,:,:] = Total_I_Transpose[:,StartIdx:EndIdx]
    end

	#Find day of seasonal year in which infection is at peak
	#Over 3rd dimension as that is day of season!
	MaxInfValBySeason,inds = findmax(Total_I_Array[4:end,:,:],dims=3)
	MaxInfIdxBySeason = map(x->x[3], inds)

    #If peak outside Sep-Feb, set amended to be very large!
    if sum(MaxInfIdxBySeason.>182) == 0 #181 days September-February. Plus account for initial value
        AmendedSummStatVal = SummStatVal_Overall
    else
		AmendedSummStatVal = Inf
    end

    return AmendedSummStatVal

end


#-------------------------------------------------------------------------------
#FUNCTION TO BE OPTIMISED
#-------------------------------------------------------------------------------
function OptimError(OptimParamSet,ObservedData,SummStatFn::Function,ModelSimnFn::Function,
                FixedModelParams)

#Inputs:
#   OptimParamSet - (vector, floats) Parameters being optimised. Take present value & get error value
#	ObservedData - (array)
#	SummStatFn::Function
#	ModelSimnFn::Function
#	FixedModelParams - (tuple)

#Outputs:
#   err - (float) Error value for tested particle

	#Initialise indicator variable that parameter constraints are met
	ParamConstraintFlag = 1

	#---------------------------------------------------------------------------
	#Check if suscepibility values meet constraints
	#---------------------------------------------------------------------------

	#Get number of susceptibility groups in use
	AgeGrpSuscepTuple = FixedModelParams[end]
	AgeGrpSuscepParamNum = AgeGrpSuscepTuple[1]

	#Get susceptibility construct in use
	AgeGrpSuscepType = AgeGrpSuscepTuple[end]

	#If three groups: 0-17, 18-64, 65+
	#If four groups: 0-17, 18-64, 65-84, 85+
	#Enforce that 18-64 is <= other age categories
	if AgeGrpSuscepType == AgeAndStrainDepSuscepFn
		if AgeGrpSuscepParamNum == 3

			#Get relevant parameters
			AdultSuscepValsIdx = [9 12 15 18]
			AdultSuscepVals = OptimParamSet[AdultSuscepValsIdx]

			#Specify values of other suscep param elements
			OtherSuscepValsIdx = [8 10; 11 13; 14 16; 17 19] #Get relevant indices
			OtherSuscepVals = OptimParamSet[OtherSuscepValsIdx]

			#For each strain, check if any adult suscep vals are greater than suscep values for other age grou
			SuscepValConstraintCheck = 0
			for ii = 1:length(AdultSuscepValsIdx)
				SuscepValConstraintCheck = SuscepValConstraintCheck + sum(AdultSuscepVals[ii].>OtherSuscepVals[ii,:])
			end

			#Assign value of ParamConstraintFlag. To be used to run model simn or not.
			if SuscepValConstraintCheck > 0
				ParamConstraintFlag = 0
			end

		elseif AgeGrpSuscepParamNum == 4

			#Get relevant parameters
			AdultSuscepValsIdx = [9 13 17 21]
			AdultSuscepVals = OptimParamSet[AdultSuscepValsIdx]

			#Specify values of other suscep param elements
			OtherSuscepValsIdx = [8 10 11; 12 14 15; 16 18 19; 20 22 23] #Get relevant indices
			OtherSuscepVals = OptimParamSet[OtherSuscepValsIdx]

			#For each strain, check if any adult suscep vals are greater than suscep values for other age grou
			SuscepValConstraintCheck = 0
			for ii = 1:length(AdultSuscepValsIdx)
				SuscepValConstraintCheck = SuscepValConstraintCheck + sum(AdultSuscepVals[ii].>OtherSuscepVals[ii,:])
			end

			#Assign value of ParamConstraintFlag. To be used to run model simn or not.
			if SuscepValConstraintCheck > 0
				ParamConstraintFlag = 0
			end
		else
			error("AgeGrpSuscepParamNum is $AgeGrpSuscepParamNum. Invalid value.")
		end
	elseif AgeGrpSuscepType == AgeDepWithStrainScalingSuscepFn || AgeGrpSuscepType == AgeOnlySuscepFn
		if AgeGrpSuscepParamNum == 3

			#Get relevant parameters
			AdultSuscepVal = OptimParamSet[9]

			#Specify values of other suscep param elements
			OtherSuscepValsIdx = [8 10] #Get relevant indices
			OtherSuscepVals = OptimParamSet[OtherSuscepValsIdx]

			#Check if any adult suscep vals are greater than suscep values for other age groups
			SuscepValConstraintCheck = sum(AdultSuscepVal.>OtherSuscepVals)

			#Assign value of ParamConstraintFlag. To be used to run model simn or not.
			if SuscepValConstraintCheck > 0
				ParamConstraintFlag = 0
			end

		elseif AgeGrpSuscepParamNum == 4

			#Get relevant parameters
			AdultSuscepVal = OptimParamSet[9]

			#Specify values of other suscep param elements
			OtherSuscepValsIdx = [8 10 11] #Get relevant indices
			OtherSuscepVals = OptimParamSet[OtherSuscepValsIdx]

			#Check if any adult suscep vals are greater than suscep values for other age groups
			SuscepValConstraintCheck = sum(AdultSuscepVal.>OtherSuscepVals)

			#Assign value of ParamConstraintFlag. To be used to run model simn or not.
			if SuscepValConstraintCheck > 0
				ParamConstraintFlag = 0
			end
		else
			error("AgeGrpSuscepParamNum is $AgeGrpSuscepParamNum. Invalid value.")
		end
	elseif AgeGrpSuscepType == PiecewiseLinearAgeSuscepFn

		#Get relevant parameters
		Age18SuscepVal = OptimParamSet[10]
		Age41SuscepVal = OptimParamSet[11]
		Age65SuscepVal = OptimParamSet[12]

		#Specify values of other suscep param elements
		OtherSuscepValsIdx = [8,9,13,14] #Get relevant indices
		OtherSuscepVals = OptimParamSet[OtherSuscepValsIdx]

		#Check if any adult suscep vals (ages 18, 41, 65) are greater than suscep values for other age groups
		SuscepValConstraintCheck = sum(Age18SuscepVal.>OtherSuscepVals) + sum(Age41SuscepVal.>OtherSuscepVals) +
										sum(Age65SuscepVal.>OtherSuscepVals)


		#If parameter values non-compatible, assign zero prior probability.
		if SuscepValConstraintCheck > 0
			ParamConstraintFlag = 0
		end

		#Check if age 45 sucep values exceeds values at ages 18 or 65
		Age41SuscepValConstraintCheck = sum(Age41SuscepVal.>Age18SuscepVal) + sum(Age41SuscepVal.>Age65SuscepVal)

		#If parameter values non-compatible, assign zero prior probability.
		if Age41SuscepValConstraintCheck > 0
			ParamConstraintFlag = 0
		end

	else
		error("AgeGrpSuscepType is $AgeGrpSuscepType. Invalid input value.")
	end

	#---------------------------------------------------------------------------
	#Check if ascertainment values meet constraints
	#---------------------------------------------------------------------------

	#Only needed if using ascertainment with age modifier
	#Get ascertainment type value.
	AscertainProbFlag = FixedModelParams[end-2]

	#Will only enter loop if AscertainProbFlag value is satisfied
	if AscertainProbFlag == SeasonWithAgeScaleStepAscertainmentFn

		#Get other relevant parameters from FixedModelParams
		FitAggAgeTuple = FixedModelParams[end-1]
		AgeBandBounds = FitAggAgeTuple[2]::Array{Int64,2}
		AscertainmentGrpNum = convert(Int64,size(AgeBandBounds,2)) #Number of columns of age bound array equiv. to number of age bands in use

		RetainedSeasonNum = FixedModelParams[end-3]::Int64

		#Get season-by-season base values
		AscertainmentAgeScalingParamNum = AscertainmentGrpNum - 1
		AccessIdxRelEnd =  RetainedSeasonNum + AscertainmentAgeScalingParamNum - 1
		AscertainSeasonVals = OptimParamSet[end-AccessIdxRelEnd:end-AscertainmentAgeScalingParamNum]


		#Get age scaling values.
		AscertainmentAgeScalingParamNum = AscertainmentGrpNum - 1
		AscertainProbAgeScaleVals = OptimParamSet[(end-AscertainmentAgeScalingParamNum+1):end]

		#Take product of age scaling values with base values
		ProductAscertainVals = AscertainSeasonVals.*AscertainProbAgeScaleVals'


		#Check if any exceed 1. If so, constraint check is failed
		AscertainValConstraintCheck = sum(ProductAscertainVals.>1)
		if AscertainValConstraintCheck > 0
			ParamConstraintFlag = 0
		end

	end

	#---------------------------------------------------------------------------
	#Assign error value of current particle set
	#---------------------------------------------------------------------------
	if ParamConstraintFlag == 1
	    #Run model with designated parameter set
	    SimnData = ModelSimnFn(FixedModelParams,OptimParamSet)

	    #Generate summary statistic for current prior set
	    err = SummStatFn(ObservedData,SimnData)
	else #Parameter set not valid, failed constraints check
		err = Inf
	end
	println("err: $err")
	println("OptimParamSet: $OptimParamSet")

    return err
end


#-------------------------------------------------------------------------------
#OPTIMISER SCHEME PARAMETER BOUNDS & INITIAL CONDITIONS
#-------------------------------------------------------------------------------

function OptimSchemeInputs_EmpDataRun1()

	#Box Constrained Optimization (18 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, #Age suscep params
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability
					 	0, 0]  #Ascertainment age band scaling
	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, #Age suscep params
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability
						10, 10] #Ascertainment age band scaling

	#Initial conditions
  InitialVals = [2.2196788520495394, 1.8113084405939324, 1.7658002345748325, 2.003630388400491,
              0.39457461519503473, 0.6268167869850989, 0.6795522747683778,
              0.9594467863493714, 0.4969590511676417, 0.7626444031824217,
              0.00728360853708827, 0.004389609292259963, 0.012773241006523885, 0.00959389273863017, 0.01264897554460983, 0.012682308198086523,
                3.5188403155380343, 4.283950925238483]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun2()

	#Box Constrained Optimization (20 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, 0, #Age suscep params
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability
					 	0, 0, 0]  #Ascertainment age band scaling
	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, 1, #Age suscep params
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability
						10, 10, 10] #Ascertainment age band scaling

	#Initial conditions
  InitialVals = [2.2, 1.8, 1.75, 2.,
              0.4, 0.6, 0.7,
              0.95, 0.5, 0.75, 0.8,
              0.007, 0.004, 0.01, 0.001, 0.01, 0.015,
              3,4,4]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun3()

	#Box Constrained Optimization (21 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0,  #Age suscep params
            0, 0, 0,    #Suscep strain scaling
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability
					 	0, 0]  #Ascertainment age band scaling
	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1,  #Age suscep params
            1, 1, 1,    #Suscep strain scaling
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability
						10, 10] #Ascertainment age band scaling

	#Initial conditions
  InitialVals = [2.2, 1.8, 1.75, 2.,
              0.4, 0.6, 0.7,
              0.95, 0.5, 0.75,
              1, 1, 1,
              0.007, 0.004, 0.01, 0.001, 0.01, 0.015,
              3,4]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun4()

	#Box Constrained Optimization (23 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, 0, #Age suscep params
            0, 0, 0,    #Suscep strain scaling
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability
					 	0, 0, 0]  #Ascertainment age band scaling
	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, 1, #Age suscep params
            1, 1, 1,    #Suscep strain scaling
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability
						10, 10, 10] #Ascertainment age band scaling

	#Initial conditions
  InitialVals = [2.2, 1.8, 1.75, 2.,
              0.4, 0.6, 0.7,
              0.95, 0.5, 0.75, 0.8,
              1, 1, 1,
              0.007, 0.004, 0.01, 0.001, 0.01, 0.015,
              3,4,4]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun5()

	#Box Constrained Optimization (27 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0,  #Suscep params A(H1N1)
					 	0, 0, 0,  #Suscep params A(H3N2)
            0, 0, 0,  #Suscep params B/Yamagata
            0, 0, 0,  #Suscep params B/Victoria
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability
					 	0, 0]  #Ascertainment age band scaling
	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
					 	1, 1, 1,  #Suscep params A(H1N1)
					 	1, 1, 1,  #Suscep params A(H3N2)
            1, 1, 1,  #Suscep params B/Yamagata
            1, 1, 1,  #Suscep params B/Victoria
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability
						10, 10] #Ascertainment age band scaling

	#Initial conditions
  InitialVals = [2.2, 1.8, 1.75, 2.,
              0.4, 0.6, 0.7,
              1, 1, 1,
              1, 1, 1,
              1, 1, 1,
              1, 1, 1,
              0.007, 0.004, 0.01, 0.001, 0.01, 0.015,
              3, 4]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun6()

	#Box Constrained Optimization (32 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, 0, #Suscep params A(H1N1)
					 	0, 0, 0, 0, #Suscep params A(H3N2)
            0, 0, 0, 0, #Suscep params B/Yamagata
            0, 0, 0, 0, #Suscep params B/Victoria
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability
					 	0, 0, 0]  #Ascertainment age band scaling
	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
					 	1, 1, 1, 1,  #Suscep params A(H1N1)
					 	1, 1, 1, 1,  #Suscep params A(H3N2)
            1, 1, 1, 1,  #Suscep params B/Yamagata
            1, 1, 1, 1,  #Suscep params B/Victoria
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability
						10, 10, 10] #Ascertainment age band scaling

	#Initial conditions
  InitialVals = [2.2, 1.8, 1.75, 2.,
              0.4, 0.6, 0.7,
              1, 1, 1, 1,
              1, 1, 1, 1,
              1, 1, 1, 1,
              1, 1, 1, 1,
              0.007, 0.004, 0.01, 0.001, 0.01, 0.015,
              3, 4, 4]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun7()

	#Box Constrained Optimization (16 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, #Age suscep params
					 	0, 0, 0, 0, 0, 0] #Ascertainment probability

	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, #Age suscep params
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2]  #Ascertainment probability

	#Initial conditions
  InitialVals = [2.2, 1.8, 1.75, 2.,
              0.4, 0.6, 0.7,
              0.95, 0.5, 0.75,
              0.007, 0.004, 0.01, 0.001, 0.01, 0.015]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end


function OptimSchemeInputs_EmpDataRun8()

	#Box Constrained Optimization (17 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, 0, #Age suscep params
					 	0, 0, 0, 0, 0, 0] #Ascertainment probability

	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, 1, #Age suscep params
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2]  #Ascertainment probability

	#Initial conditions
  InitialVals = [2.2, 1.8, 1.75, 2.,
              0.4, 0.6, 0.7,
              0.95, 0.5, 0.75, 0.8,
              0.007, 0.004, 0.01, 0.001, 0.01, 0.015]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun9()

	#Box Constrained Optimization (19 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0,  #Age suscep params
            0, 0, 0,    #Suscep strain scaling
					 	0, 0, 0, 0, 0, 0] #Ascertainment probability

	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1,  #Age suscep params
            1, 1, 1,    #Suscep strain scaling
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2] #Ascertainment probability

	#Initial conditions
  InitialVals = [2.2, 1.8, 1.75, 2.,
              0.4, 0.6, 0.7,
              0.95, 0.5, 0.75,
              1, 1, 1,
              0.007, 0.004, 0.01, 0.001, 0.01, 0.015]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun10()

	#Box Constrained Optimization (20 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, 0, #Age suscep params
            0, 0, 0,    #Suscep strain scaling
					 	0, 0, 0, 0, 0, 0] #Ascertainment probability

	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, 1, #Age suscep params
            1, 1, 1,    #Suscep strain scaling
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2]  #Ascertainment probability

	#Initial conditions
  InitialVals = [2.2, 1.8, 1.75, 2.,
              0.4, 0.6, 0.7,
              0.95, 0.5, 0.75, 0.8,
              1, 1, 1,
              0.007, 0.004, 0.01, 0.001, 0.01, 0.015]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun11()

	#Box Constrained Optimization (25 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0,  #Suscep params A(H1N1)
					 	0, 0, 0,  #Suscep params A(H3N2)
            0, 0, 0,  #Suscep params B/Yamagata
            0, 0, 0,  #Suscep params B/Victoria
					 	0, 0, 0, 0, 0, 0] #Ascertainment probability

	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
					 	1, 1, 1,  #Suscep params A(H1N1)
					 	1, 1, 1,  #Suscep params A(H3N2)
            1, 1, 1,  #Suscep params B/Yamagata
            1, 1, 1,  #Suscep params B/Victoria
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2]  #Ascertainment probability

	#Initial conditions
  InitialVals = [2.2, 1.8, 1.75, 2.,
              0.4, 0.6, 0.7,
              1, 1, 1,
              1, 1, 1,
              1, 1, 1,
              1, 1, 1,
              0.007, 0.004, 0.01, 0.001, 0.01, 0.015]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun12()

	#Box Constrained Optimization (29 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, 0, #Suscep params A(H1N1)
					 	0, 0, 0, 0, #Suscep params A(H3N2)
            0, 0, 0, 0, #Suscep params B/Yamagata
            0, 0, 0, 0, #Suscep params B/Victoria
					 	0, 0, 0, 0, 0, 0] #Ascertainment probability

	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
					 	1, 1, 1, 1,  #Suscep params A(H1N1)
					 	1, 1, 1, 1,  #Suscep params A(H3N2)
            1, 1, 1, 1,  #Suscep params B/Yamagata
            1, 1, 1, 1,  #Suscep params B/Victoria
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2]  #Ascertainment probability

	#Initial conditions
  InitialVals = [2.2, 1.8, 1.75, 2.,
              0.4, 0.6, 0.7,
              1, 1, 1, 1,
              1, 1, 1, 1,
              1, 1, 1, 1,
              1, 1, 1, 1,
              0.007, 0.004, 0.01, 0.001, 0.01, 0.015]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end


function OptimSchemeInputs_EmpDataRun13()

	#Box Constrained Optimization (21 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, #Age suscep params
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability
					 	0, 0, 0, 0, 0]  #Ascertainment piecewise linear function, values at knot pts
	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, #Age suscep params
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability
						1, 1, 1, 1, 1] #Ascertainment piecewise linear function, values at knot pts

	#Initial conditions
 InitialVals = [1.9939664237050447,	2.1068352296781976,	1.9413227004922382,	1.9617936733090868,
                   0.5539800160919733,	1.0,	0.7815061435209223,
                   0.9602970885070373,	0.5285686182318795,	0.9999734396063796,
                   0.04315162349420877,	0.02558525863595175,	0.07465446048576288,	0.08002236173704035,	0.050169560645328765,	0.1721714771320867,
                   0.02945603686229785,	0.015478865163097714,	0.03987283196527448,	0.13715723518984838,	0.24315700851525884]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun14()

	#Box Constrained Optimization (22 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, 0, #Age suscep params
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability
						0, 0, 0, 0, 0]  #Ascertainment piecewise linear function, values at knot pts
	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, 1, #Age suscep params
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability
						1, 1, 1, 1, 1] #Ascertainment piecewise linear function, values at knot pts

	#Initial conditions
  InitialVals = [2.2214475787802463,	2.3620810513179764,	2.1827770026831135,	2.06592034263203,
                    0.6354301935236959,	1.0,	0.07500695241787407,
                    0.7696120295247926,	0.44566103545692837,	0.5387237142971943,	1.0,
                    0.034195732620690665,	0.023902267682289865,	0.05642106680687132,	0.0838380506284142,	0.03890722366472773,	0.15691021866954943,
                    0.03676016766868808,	0.01746672509809063,	0.04608295742620308,	0.19281484759471984,	0.42780903581447266]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end


function OptimSchemeInputs_EmpDataRun15()

	#Box Constrained Optimization (32 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, #Age suscep params
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability (type A)
					 	0, 0, 0, 0, 0,  #Ascertainment piecewise linear function, values at knot pts (type A)
						0, 0, 0, 0, 0, 0, #Ascertainment probability (type B)
						0, 0, 0, 0, 0]  #Ascertainment piecewise linear function, values at knot pts (type B)
	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, #Age suscep params
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability (type A)
						1, 1, 1, 1, 1, #Ascertainment piecewise linear function, values at knot pts (type A)
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability (type B)
						1, 1, 1, 1, 1] #Ascertainment piecewise linear function, values at knot pts (type B)
	#Initial conditions
	InitialVals = [1.9740321122089362,	2.326461205501085,	1.8301051838889917,	1.7462776295274212,
                    0.7002228890471642,	0.5508488638011485,	0.0,
                    0.9749941565417115,	0.562484156171011,	1.0,
                    0.011255781053653971,	0.03252936206660889,	0.02839427402403924,	0.07414005472191315,	0.020446241196739195,	0.0924674814407407,
                    0.026571720243574872,	0.016388021043706093,	0.04250158882900217,	0.12342828903593253,	0.18277171562290512,
                    0.2,	0.0008634797040275315,	0.2,	0.019194946936152802,	0.14788997638005327,	0.05831368832926678,
                    0.0943077096921879,	0.045873819950463954,	0.12380245007685785,	0.42534155264708434,	0.8689049016701581]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun16()

	#Box Constrained Optimization (33 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, 0, #Age suscep params
						0, 0, 0, 0, 0, 0, #Ascertainment probability (type A)
					 	0, 0, 0, 0, 0,  #Ascertainment piecewise linear function, values at knot pts (type A)
						0, 0, 0, 0, 0, 0, #Ascertainment probability (type B)
						0, 0, 0, 0, 0]  #Ascertainment piecewise linear function, values at knot pts (type B)
	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, 1, #Age suscep params
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability (type A)
						1., 1., 1., 1., 1., #Ascertainment piecewise linear function, values at knot pts (type A)
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability (type B)
						1., 1., 1., 1., 1.] #Ascertainment piecewise linear function, values at knot pts (type B)

	#Initial conditions
	InitialVals = [2.2310979065770447,	2.2111231698424243,	2.295252256495237,	2.2464908787273443,
                  1.0,	1.0,	0.15123761609526928,
                  0.47008657811249815,	0.4505850732723051,	1.0,	1.0,
                  0.2,	0.2,	0.2,	0.2,	0.2,	0.2,
                  1.0,	1.0,	1.0,	1.0,	1.0,
                  0.2,	0.05089665766458926,	0.2,	0.2,	0.2,	0.2,
                  1.0,	1.0,	1.0,	1.0,	1.0]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun17()

	#Box Constrained Optimization (54 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, #Age suscep params
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability for A(H1N1)
					 	0, 0, 0, 0, 0,  #Ascertainment piecewise linear function, values at knot pts for A(H1N1)
						0, 0, 0, 0, 0, 0, #Ascertainment probability for A(H3N2)
						0, 0, 0, 0, 0,  #Ascertainment piecewise linear function, values at knot pts for A(H3N2)
            0, 0, 0, 0, 0, 0, #Ascertainment probability for B/Yam
					 	0, 0, 0, 0, 0,  #Ascertainment piecewise linear function, values at knot pts for B/Yam
						0, 0, 0, 0, 0, 0, #Ascertainment probability for B/Vic
						0, 0, 0, 0, 0]  #Ascertainment piecewise linear function, values at knot pts for B/Vic

	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, #Age suscep params
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability for A(H1N1)
						1, 1, 1, 1, 1, #Ascertainment piecewise linear function, values at knot pts for A(H1N1)
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability for A(H3N2)
						1, 1, 1, 1, 1, #Ascertainment piecewise linear function, values at knot pts for A(H3N2)
            0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability for B/Yam
						1, 1, 1, 1, 1, #Ascertainment piecewise linear function, values at knot pts for B/Yam
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability for B/Vic
						1, 1, 1, 1, 1] #Ascertainment piecewise linear function, values at knot pts for B/Vic

	#Initial conditions
  InitialVals = [2.274448491888507,	2.258158017140962,	2.362446523200594,	2.3131104556272213,
						0.5606209038465599,	1.0,	0.5922058843950375,
						0.5780275395000632,	0.26054267309087853,	0.7466894184843647,
						0.2,	0.2,	0.2,	0.2,	0.14184177557354255,	0.2,
						1.0,	1.0,	1.0,	1.0,	1.0,
						0.2,	0.2,	0.2,	0.2,	0.2,	0.2,
						1.0,	1.0,	1.0,	1.0,	1.0,
						0.2,	0.046003166290085974,	0.2,	0.2,	0.2,	0.2,
						1.0,	1.0,	1.0,	1.0,	1.0,
						0.2,	0.04158775598040706,	0.1857477928912368,	0.2,	0.10994220719744964,	0.2,
						1.0,	1.0,	1.0,	1.0,	1.0]
	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun18()

	#Box Constrained Optimization (55 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, 0, #Age suscep params
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability for A(H1N1)
					 	0, 0, 0, 0, 0,  #Ascertainment piecewise linear function, values at knot pts for A(H1N1)
						0, 0, 0, 0, 0, 0, #Ascertainment probability for A(H3N2)
						0, 0, 0, 0, 0,  #Ascertainment piecewise linear function, values at knot pts for A(H3N2)
            0, 0, 0, 0, 0, 0, #Ascertainment probability for B/Yam
					 	0, 0, 0, 0, 0,  #Ascertainment piecewise linear function, values at knot pts for B/Yam
						0, 0, 0, 0, 0, 0, #Ascertainment probability for B/Vic
						0, 0, 0, 0, 0]  #Ascertainment piecewise linear function, values at knot pts for B/Vic

	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, 1, #Age suscep params
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability for A(H1N1)
						1, 1, 1, 1, 1, #Ascertainment piecewise linear function, values at knot pts for A(H1N1)
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability for A(H3N2)
						1, 1, 1, 1, 1, #Ascertainment piecewise linear function, values at knot pts for A(H3N2)
            0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability for B/Yam
						1, 1, 1, 1, 1, #Ascertainment piecewise linear function, values at knot pts for B/Yam
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability for B/Vic
						1, 1, 1, 1, 1] #Ascertainment piecewise linear function, values at knot pts for B/Vic

	#Initial conditions
  InitialVals = [2.1726125292752223,	2.128720227128342,	2.2095525202647086,	2.178573391245209,
                    1.0,	1.0,	0.6531326228822851,
                    0.6157924341671746,	0.2766656410030943,	1.0,	1.0,
                    0.2,	0.2,	0.2,	0.2,	0.10451568338239968,	0.2,
                    1.0,	1.0,	1.0,	1.0,	1.0,
                    0.2,	0.2,	0.2,	0.2,	0.2,	0.2,
                    1.0,	1.0,	1.0,	1.0,	1.0,
                    0.2,	0.058684256181240374,	0.2,	0.2,	0.2,	0.2,
                    1.0,	1.0,	1.0,	1.0,	1.0,
                    0.2,	0.03900280196461758,	0.17252535644051326,	0.2,	0.10629982554875546,	0.2,
                    1.0,	1.0,	1.0,	1.0,	1.0,]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun19()

	#Box Constrained Optimization (25 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
					 	0, 0, 0, 0, 0, 0, 0, #Age suscep params (at knot values)
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability
					 	0, 0, 0, 0, 0]  #Ascertainment piecewise linear function, values at knot pts
	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, 1, 1, 1, 1, #Age suscep params (at knot values)
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability
						1, 1, 1, 1, 1] #Ascertainment piecewise linear function, values at knot pts

	#Initial conditions
 InitialVals = [1.2144032869540302,	1.9330495952833904,	1.344640819775389,	1.2438211237737278,
                   0.5695642123154823,	0.2590510032133626,	0.8728146492386213,
                   0.8509179184838518,	0.5915134887598517,	0.5258892538367392,	0.33761177184318003,	0.457733131383213,	0.7202435616124259,	0.8999282036974463,
                   0.013930841105384702,	0.02200890482040792,	0.04944954010371311,	0.04838159033036072,	0.02683638230592086,	0.08801491076842648,
                   0.6409846072419105,	0.5081652874815161,	0.003991068820820454,	0.962359529977683,	0.09043325152336523]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

function OptimSchemeInputs_EmpDataRun20()

	#Box Constrained Optimization (36 parameters)
	ParamLowerBounds = [1., 1., 1., 1., #R_0 values
					 	0, 0, 0, #Exp hist params
						0, 0, 0, 0, 0, 0, 0, #Age suscep params (at knot values)
					 	0, 0, 0, 0, 0, 0, #Ascertainment probability (type A)
					 	0, 0, 0, 0, 0,  #Ascertainment piecewise linear function, values at knot pts (type A)
						0, 0, 0, 0, 0, 0, #Ascertainment probability (type B)
						0, 0, 0, 0, 0]  #Ascertainment piecewise linear function, values at knot pts (type B)
	ParamUpperBounds = [3, 3, 3, 3, #R_0 values
						1, 1, 1, #Exp hist params
						1, 1, 1, 1, 1, 1, 1, #Age suscep params (at knot values)
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability (type A)
						1, 1, 1, 1, 1, #Ascertainment piecewise linear function, values at knot pts (type A)
						0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  #Ascertainment probability (type B)
						1, 1, 1, 1, 1] #Ascertainment piecewise linear function, values at knot pts (type B)
	#Initial conditions
	InitialVals = [1.2144032869540302,	1.9330495952833904,	1.344640819775389,	1.2438211237737278,
                   0.5695642123154823,	0.2590510032133626,	0.8728146492386213,
                   0.8509179184838518,	0.5915134887598517,	0.5258892538367392,	0.33761177184318003,	0.457733131383213,	0.7202435616124259,	0.8999282036974463,
                   0.013930841105384702,	0.02200890482040792,	0.04944954010371311,	0.04838159033036072,	0.02683638230592086,	0.08801491076842648,
                   0.6409846072419105,	0.5081652874815161,	0.003991068820820454,	0.962359529977683,	0.09043325152336523,
            			  0.013930841105384702,	0.02200890482040792,	0.04944954010371311,	0.04838159033036072,	0.02683638230592086,	0.08801491076842648,
                   0.6409846072419105,	0.5081652874815161,	0.003991068820820454,	0.962359529977683,	0.09043325152336523]

	return ParamLowerBounds,ParamUpperBounds,InitialVals
end

#-------------------------------------------------------------------------------
#FN TO CARRY OUT OPTIMISATION
#-------------------------------------------------------------------------------
function RunOptimiser(OutputFileName,ObservedData,  #As named
					SeasonsToSimulate,AscertainProbFn,ExpHistVaccType,
					FitAggAgeTuple,AgeGrpSuscepTuple,   #Inputs for FixedModelParamGenFn
						ModelSimnFn,SummStatFn,  #Function handle equivalents
						ParamLowerBounds,ParamUpperBounds,InitialVals,OptimIterationNum)  #Optimiser variables

	#--------------------------------------------------------------------------
	# GENERATE FIXEDMODELPARAMS VARIABLE, INPUT FOR OPTIMISATION FUNCTION
	#--------------------------------------------------------------------------
	FixedModelParams = FixedModelParamGenFn(SeasonsToSimulate,AscertainProbFn,ExpHistVaccType,
											FitAggAgeTuple,AgeGrpSuscepTuple,ObservedData)

	#--------------------------------------------------------------------------
	#RUN OPTIMISATION
	#--------------------------------------------------------------------------
	res = optimize(OptimParamSet -> OptimError(OptimParamSet, ObservedData,SummStatFn::Function,ModelSimnFn::Function,
	               FixedModelParams), InitialVals,
					ParticleSwarm(; lower = ParamLowerBounds, upper = ParamUpperBounds, n_particles = 3),
					Optim.Options(iterations = OptimIterationNum,store_trace=true, extended_trace=true))
	#res = optimize(OptimParamSet -> OptimError(OptimParamSet, ObservedData,SummStatFn::Function,ModelSimnFn::Function,
	# 				                FixedModelParams), InitialVals,
	# 								NelderMead(; lower = ParamLowerBounds, upper = ParamUpperBounds), Optim.Options(iterations = 3))


	#--------------------------------------------------------------------------
	#ASSIGN OUTPUTS
	#--------------------------------------------------------------------------

	#Assign best-fit parameters to separate variable
	BestFitParams = Optim.minimizer(res)
	ParamTrace = Optim.x_trace(res)
	ErrorTrace = Optim.f_trace(res)

	#Assign to output file
	writedlm(OutputFileName[1],ParamTrace)
	writedlm(OutputFileName[2],ErrorTrace)

end

function FMAlt_OptimiseCallFn(ARGS)

#Take command line arguments, ARGS, assign to variable names
#To convert strings to numbers, use parse

#--------------------------------------------------------------------------
# LIST SIMULATION ID & SAVE NAME DETAILS
#--------------------------------------------------------------------------
OptimRunID = ARGS[1]
OutputFileName = ["FMAlt_OptimFitOutputFiles/OptimiserParamsTrace#$(OptimRunID).txt",
				   "FMAlt_OptimFitOutputFiles/OptimiserErrorTrace#$(OptimRunID).txt"]

#--------------------------------------------------------------------------
# GATHER OBSERVED DATA
#--------------------------------------------------------------------------
ObvsDataFileName = ARGS[2]
ObservedDataTemp = readdlm(ObvsDataFileName,',')
MaxSeasonNumToSimulate = size(ObservedDataTemp,1)

#--------------------------------------------------------------------------
# SPECIFY MODEL RUN AND SUMMARY STATISTIC FUNCTION
#--------------------------------------------------------------------------

#Specify function name inputs
s_1 = Symbol(ARGS[3]) #Convert string to Symbol
s_2 = Symbol(ARGS[4]) #Convert string to Symbol

#Make Symbols callable functions
ModelSimnFn = getfield(Main, s_1) #Fn specifying model simulation
SummStatFn = getfield(Main, s_2) #Fn computing error measure, compare observed to simulated data

#--------------------------------------------------------------------------
# SET UP FIXED MODEL PARAMS
#--------------------------------------------------------------------------
#FUNCTION TO OUTPUT FIXED MODEL PARAMS?

#Total simulation run time
SeasonsToSimulate = parse(Int64, ARGS[5])

#--------------------------------------------------------------------------
# Declare the type of ascertainment probabiity modification to be formed
#--------------------------------------------------------------------------
#Specify function name inputs
AscertainProbFnSymb = Symbol(ARGS[6]) #Convert string to Symbol

#Make Symbols callable functions
AscertainProbFn = getfield(Main, AscertainProbFnSymb) #Fn specifying model simulation

#--------------------------------------------------------------------------
#Select exposure history susceptibility modification form for vaccine-related states
#0 - Fixed values every season
#1 - Relate to previous season vaccine efficacy
#2 - Relate to previous season vaccine efficacy AND age band
#--------------------------------------------------------------------------
ExpHistVaccType = parse(Int64, ARGS[7])

if ExpHistVaccType != 0 && ExpHistVaccType !=1 && ExpHistVaccType !=2
	error("ExpHistVaccType must take value 0, 1 or 2. Current value is $(ExpHistVaccType)")
end

#--------------------------------------------------------------------------
# Declare bounds for age category groupings (for aggregating case counts &
# applying ascertainment probabilitiy, if required)
#--------------------------------------------------------------------------
FitAggAgeFlag = parse(Int64, ARGS[8]) #Indicator variable.
# Set to: -> 1 to fit to aggregetated age band counts
#         -> 0 to fit to single yr age classes

#Error check
if FitAggAgeFlag != 0 && FitAggAgeFlag !=1
   error("FitAggAgeFlag must take value 0 or 1. Current value is $FitAggAgeFlag")
end

#Use Meta.parse for evaluating a string into a vector
AgeBandLowerBounds = eval(Meta.parse(ARGS[9]))
AgeBandUpperBounds = eval(Meta.parse(ARGS[10]))

#Concantenate bounds into a single array
AgeBandBounds = [AgeBandLowerBounds;AgeBandUpperBounds]
				#vertically stack age category lower bounds above upper bounds

#Place FitAggAgeFlag&AgeBandBoundsi nto a cell/tuple
FitAggAgeTuple = [FitAggAgeFlag,AgeBandBounds]

#--------------------------------------------------------------------------
# Declare bounds for age susceptibility groupings
#--------------------------------------------------------------------------
AgeSuscepLowerBounds = eval(Meta.parse(ARGS[11]))
AgeSuscepUpperBounds = eval(Meta.parse(ARGS[12]))
AgeGrpSuscepParamNum = length(AgeSuscepLowerBounds)

#Select age&strain susceptibility modification form
# AgeAndStrainDepSuscepFn - Age group & strain specific
# AgeDepWithStrainScalingSuscepFn - Age group specific, with strain modifier
# StrainDepWithAgeScalingSuscepFn - Strain specific, with age group modifier
# AgeOnlySuscepFn - Age group speicifc only (independent of strain)
AgeSuscepTypeSymb = Symbol(ARGS[13]) #Convert string to Symbol

#Make Symbols callable functions
AgeSuscepType = getfield(Main, AgeSuscepTypeSymb) #Fn specifying model simulation


#Concantenate into a cell tuple
AgeGrpSuscepTuple = [AgeGrpSuscepParamNum,AgeSuscepLowerBounds,AgeSuscepUpperBounds,AgeSuscepType]

#--------------------------------------------------------------------------
# SPECIFY PORTION OF INPUT DATA TO BE FIT TO
# Put data into 3D array: Row per season, column per age group, slice (3rd
# dim.) per strain.
#--------------------------------------------------------------------------
NumOfStrains = 4

#Get number of age bands being fitted to
if FitAggAgeFlag == 0 #Fitting to single age classes
    AgeBandFitNum = 101
else #Fit to user-defined age bands
    AgeBandFitNum = length(AgeBandLowerBounds)
end

ObservedDataReshaped = reshape(ObservedDataTemp,(MaxSeasonNumToSimulate,AgeBandFitNum,NumOfStrains))

#Portion of input data corresponding to 2012/2013 onwards
ObservedData = ObservedDataReshaped[1:end,:,:]

#--------------------------------------------------------------------------
# DEFINE OPTIMISATION SCHEME INPUTS
#--------------------------------------------------------------------------

#Specify function name inputs
s_3 = Symbol(ARGS[14]) #Convert string to Symbol

#Make Symbols callable functions
OptimSchemeInputFn = getfield(Main, s_3) #Fn specifying model simulation

#Get optimsiation scheme variables from input function
ParamLowerBounds,ParamUpperBounds,InitialVals = OptimSchemeInputFn()

#Number of iterations optimisation scheme should perform
OptimIterationNum = parse(Int64, ARGS[15])

#--------------------------------------------------------------------------
# RUN OPTIMISATION
#--------------------------------------------------------------------------
RunOptimiser(OutputFileName,ObservedData,
		 	SeasonsToSimulate,AscertainProbFn,ExpHistVaccType,
		 	FitAggAgeTuple,AgeGrpSuscepTuple,
		 	ModelSimnFn,SummStatFn,
		 	ParamLowerBounds,ParamUpperBounds,InitialVals,OptimIterationNum)
end

#--------------------------------------------------------------------------
# PASS TO FUNCTION
#--------------------------------------------------------------------------
FMAlt_OptimiseCallFn(ARGS)
