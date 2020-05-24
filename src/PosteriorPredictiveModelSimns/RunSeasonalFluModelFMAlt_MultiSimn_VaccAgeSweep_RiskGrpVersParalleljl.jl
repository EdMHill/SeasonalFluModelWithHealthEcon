#Purpose:
#Script to multiple model simulations with parameters read from file,
#for portion of the populaton based on risk group status.
#Set up to run in parallel

#Model specifics:
#SEIR influenza deterministic transmission dynamic model ODEs
#Includes age & strain structure.
#Exposure history (for previous season only) included
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# LOAD REQUIRED ENVIRONMENT (FILE PATH BASED ON CWD BEING THE LOCATION THIS FILE RESIDES)
#--------------------------------------------------------------------------
@everywhere using Pkg
@everywhere Pkg.activate("../")

#-------------------------------------------------------------------------------
# LOAD REQUIRED PACKAGES
#-------------------------------------------------------------------------------
@everywhere using Combinatorics
@everywhere using DifferentialEquations
@everywhere using DataFrames
@everywhere using XLSX
@everywhere using DelimitedFiles
@everywhere using LinearAlgebra
@everywhere using MAT
@everywhere using Distributed

#-------------------------------------------------------------------------------
###  VACCINATION UPTAKE FUNCTIONS
#-------------------------------------------------------------------------------

#Revised paediatric and elder age programme age bounds
# Two years of age resolution. Strats span age 2 upwards
@everywhere function VaccUptake_ExpandSchoolAndElderly_NoUnder2s_RecovRisk(SeasonsToSimulate,JobIdx,VaccAgeBoundsData,SimnRunType)

	#Initialise tuple per season
	#Per tupel, 2D array
	#rows for age (0 to 100+), cols for calendar day of year
	VaccUptakeBySeason = Vector{Array}(undef,SeasonsToSimulate)

	# Import the data - Pandemic flu vacc (2009/2010 season)
	#Load data based on risk group type
	if SimnRunType == 2
		#AT RISK
		PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc.xlsx","At risk","C3:NC103")
	elseif SimnRunType == 3
		#LOW RISK
		PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc.xlsx","Low risk","C3:NC103")
	end
	# Collate into Array, Assign to storage cell
	PandemicFluVaccUptakeArray = Array{Float64, 2}(PandemicFluVaccUptake)
	#PandemicFluVaccUptakeArray[19:end,:] .= 0 #Set specified age groups to have no vaccination
	VaccUptakeBySeason[1] = PandemicFluVaccUptakeArray

	#Import the data - Sesonal flu vacc
	SheetNames = ["2010_2011","2011_2012","2012_2013","2013_2014",
		"2014_2015","2015_2016","2016_2017","2017_2018"]

	#Iterate over each season
	for ii = 1:length(SheetNames)

		#Load data based on risk group type
		if SimnRunType == 2
		    #AT RISK
				#Load at-risk only vaccine programme data
				AtRiskProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_AtRisk_HistoricalAtRiskOnly.xlsx","$(SheetNames[ii])","C3:NC103")

		        #Load vaccine uptake data with additinal age groups covered
				ExpandedProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_AtRisk_AllPopnHighVaccUptake.xlsx","$(SheetNames[ii])","C3:NC103")
		elseif SimnRunType == 3
		    #LOW RISK
			#Load at-risk only vaccine programme data
			AtRiskProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_LowRisk_HistoricalAtRiskOnly.xlsx","$(SheetNames[ii])","C3:NC103")

			#Load vaccine uptake data with additinal age groups covered
			ExpandedProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_LowRisk_AllPopnHighVaccUptake.xlsx","$(SheetNames[ii])","C3:NC103")
	end

		# Collate input data into Arrays
		AtRiskVaccUptakeArray = Array{Float64,2}(AtRiskProgSeasonalFluVaccUptake)
		ExpandedProgVaccUptakeArray = Array{Float64,2}(ExpandedProgSeasonalFluVaccUptake)

		#Assign rows based on JobIdx (read from text file)
		CurrentJobIdx_VaccAgeBoundsData = VaccAgeBoundsData[JobIdx,:]
		YoungClusterUB = VaccAgeBoundsData[JobIdx,1] #First entry the upper bound for young age cluster
		ElderClusterLB = VaccAgeBoundsData[JobIdx,2]   #Second entry for lower bound of elder age cluster

		#Populate vaccine uptake array for current season iteration
		VaccUptakeBySeasonArray = AtRiskVaccUptakeArray #Initial data. Match at risk only programme

		if isnan(YoungClusterUB)  #Elder cluster only
			ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
			VaccUptakeBySeasonArray[(ElderClusterLB+1):end,:] = ExpandedProgVaccUptakeArray[(ElderClusterLB+1):end,:] #Amend elder age cluster rows
		elseif isnan(ElderClusterLB) #Young ages only
			YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
			VaccUptakeBySeasonArray[3:(YoungClusterUB+1),:] = ExpandedProgVaccUptakeArray[3:(YoungClusterUB+1),:] #Amend young age cluster rows, from age 4 upwards
		else #Both ages
			YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
			ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
			VaccUptakeBySeasonArray[3:(YoungClusterUB+1),:] = ExpandedProgVaccUptakeArray[3:(YoungClusterUB+1),:] #Amend young age cluster rows, from age 4 upwards
			VaccUptakeBySeasonArray[(ElderClusterLB+1):end,:] = ExpandedProgVaccUptakeArray[(ElderClusterLB+1):end,:] #Amend elder age cluster rows
		end

		#Assign populated array to output tuple
		VaccUptakeBySeason[ii+1] = VaccUptakeBySeasonArray  #Add 1 to idx as first entry is for 2009/2010 season
	end
return VaccUptakeBySeason
end

#Revised paediatric and elder age programme age bounds
# Two years of age resolution. Strats span age 4 upwards
@everywhere function VaccUptake_ExpandSchoolAndElderly_NoUnder4s_RecovRisk(SeasonsToSimulate,JobIdx,VaccAgeBoundsData,SimnRunType)

	#Initialise tuple per season
	#Per tupel, 2D array
	#rows for age (0 to 100+), cols for calendar day of year
	VaccUptakeBySeason = Vector{Array}(undef,SeasonsToSimulate)

	# Import the data - Pandemic flu vacc (2009/2010 season)
	#Load data based on risk group type
	if SimnRunType == 2
		#AT RISK
		PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc.xlsx","At risk","C3:NC103")
	elseif SimnRunType == 3
		#LOW RISK
		PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc.xlsx","Low risk","C3:NC103")
	end
	# Collate into Array, Assign to storage cell
	PandemicFluVaccUptakeArray = Array{Float64, 2}(PandemicFluVaccUptake)
	#PandemicFluVaccUptakeArray[19:end,:] .= 0 #Set specified age groups to have no vaccination
	VaccUptakeBySeason[1] = PandemicFluVaccUptakeArray

	#Import the data - Sesonal flu vacc
	SheetNames = ["2010_2011","2011_2012","2012_2013","2013_2014",
		"2014_2015","2015_2016","2016_2017","2017_2018"]

	#Iterate over each season
	for ii = 1:length(SheetNames)

		#Load data based on risk group type
		if SimnRunType == 2
		    #AT RISK
				#Load at-risk only vaccine programme data
				AtRiskProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_AtRisk_HistoricalAtRiskOnly.xlsx","$(SheetNames[ii])","C3:NC103")

		        #Load vaccine uptake data with additinal age groups covered
				ExpandedProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_AtRisk_AllPopnHighVaccUptake.xlsx","$(SheetNames[ii])","C3:NC103")
		elseif SimnRunType == 3
		    #LOW RISK
			#Load at-risk only vaccine programme data
			AtRiskProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_LowRisk_HistoricalAtRiskOnly.xlsx","$(SheetNames[ii])","C3:NC103")

			#Load vaccine uptake data with additinal age groups covered
			ExpandedProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_LowRisk_AllPopnHighVaccUptake.xlsx","$(SheetNames[ii])","C3:NC103")
	end

		# Collate input data into Arrays
		AtRiskVaccUptakeArray = Array{Float64,2}(AtRiskProgSeasonalFluVaccUptake)
		ExpandedProgVaccUptakeArray = Array{Float64,2}(ExpandedProgSeasonalFluVaccUptake)

		#Assign rows based on JobIdx (read from text file)
		CurrentJobIdx_VaccAgeBoundsData = VaccAgeBoundsData[JobIdx,:]
		YoungClusterUB = VaccAgeBoundsData[JobIdx,1] #First entry the upper bound for young age cluster
		ElderClusterLB = VaccAgeBoundsData[JobIdx,2]   #Second entry for lower bound of elder age cluster

		#Populate vaccine uptake array for current season iteration
		VaccUptakeBySeasonArray = AtRiskVaccUptakeArray #Initial data. Match at risk only programme

		if isnan(YoungClusterUB)  #Elder cluster only
			ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
			VaccUptakeBySeasonArray[(ElderClusterLB+1):end,:] = ExpandedProgVaccUptakeArray[(ElderClusterLB+1):end,:] #Amend elder age cluster rows
		elseif isnan(ElderClusterLB) #Young ages only
			YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
			VaccUptakeBySeasonArray[5:(YoungClusterUB+1),:] = ExpandedProgVaccUptakeArray[5:(YoungClusterUB+1),:] #Amend young age cluster rows, from age 4 upwards
		else #Both ages
			YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
			ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
			VaccUptakeBySeasonArray[5:(YoungClusterUB+1),:] = ExpandedProgVaccUptakeArray[5:(YoungClusterUB+1),:] #Amend young age cluster rows, from age 4 upwards
			VaccUptakeBySeasonArray[(ElderClusterLB+1):end,:] = ExpandedProgVaccUptakeArray[(ElderClusterLB+1):end,:] #Amend elder age cluster rows
		end

		#Assign populated array to output tuple
		VaccUptakeBySeason[ii+1] = VaccUptakeBySeasonArray  #Add 1 to idx as first entry is for 2009/2010 season
	end
return VaccUptakeBySeason
end

#Revised paediatric and elder age programme age bounds
# Two years of age resolution. Strats span the age space.
@everywhere function VaccUptake_ExpandSchoolAndElderly_SpanAgeSpace_RecovRisk(SeasonsToSimulate,JobIdx,VaccAgeBoundsData,SimnRunType)

	#Initialise tuple per season
	#Per tupel, 2D array
	#rows for age (0 to 100+), cols for calendar day of year
	VaccUptakeBySeason = Vector{Array}(undef,SeasonsToSimulate)

	# Import the data - Pandemic flu vacc (2009/2010 season)
	#Load data based on risk group type
	if SimnRunType == 2
		#AT RISK
		PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc.xlsx","At risk","C3:NC103")
	elseif SimnRunType == 3
		#LOW RISK
		PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc.xlsx","Low risk","C3:NC103")
	end
	# Collate into Array, Assign to storage cell
	PandemicFluVaccUptakeArray = Array{Float64, 2}(PandemicFluVaccUptake)
	#PandemicFluVaccUptakeArray[19:end,:] .= 0 #Set specified age groups to have no vaccination
	VaccUptakeBySeason[1] = PandemicFluVaccUptakeArray

	#Import the data - Sesonal flu vacc
	SheetNames = ["2010_2011","2011_2012","2012_2013","2013_2014",
		"2014_2015","2015_2016","2016_2017","2017_2018"]

	#Iterate over each season
	for ii = 1:length(SheetNames)

		#Load data based on risk group type
		if SimnRunType == 2
		    #AT RISK
				#Load at-risk only vaccine programme data
				AtRiskProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_AtRisk_HistoricalAtRiskOnly.xlsx","$(SheetNames[ii])","C3:NC103")

		        #Load vaccine uptake data with additinal age groups covered
				ExpandedProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_AtRisk_AllPopnHighVaccUptake.xlsx","$(SheetNames[ii])","C3:NC103")
		elseif SimnRunType == 3
		    #LOW RISK
			#Load at-risk only vaccine programme data
			AtRiskProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_LowRisk_HistoricalAtRiskOnly.xlsx","$(SheetNames[ii])","C3:NC103")

			#Load vaccine uptake data with additinal age groups covered
			ExpandedProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_LowRisk_AllPopnHighVaccUptake.xlsx","$(SheetNames[ii])","C3:NC103")
	end

		# Collate input data into Arrays
		AtRiskVaccUptakeArray = Array{Float64,2}(AtRiskProgSeasonalFluVaccUptake)
		ExpandedProgVaccUptakeArray = Array{Float64,2}(ExpandedProgSeasonalFluVaccUptake)

		#Assign rows based on JobIdx (read from text file)
		CurrentJobIdx_VaccAgeBoundsData = VaccAgeBoundsData[JobIdx,:]
		YoungClusterUB = VaccAgeBoundsData[JobIdx,1] #First entry the upper bound for young age cluster
		ElderClusterLB = VaccAgeBoundsData[JobIdx,2]   #Second entry for lower bound of elder age cluster

		#Populate vaccine uptake array for current season iteration
		VaccUptakeBySeasonArray = AtRiskVaccUptakeArray #Initial data. Match at risk only programme

		if isnan(YoungClusterUB)  #Elder cluster only
			ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
			VaccUptakeBySeasonArray[(ElderClusterLB+1):end,:] = ExpandedProgVaccUptakeArray[(ElderClusterLB+1):end,:] #Amend elder age cluster rows
		elseif isnan(ElderClusterLB) #Young ages only
			YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
			VaccUptakeBySeasonArray[1:(YoungClusterUB+1),:] = ExpandedProgVaccUptakeArray[1:(YoungClusterUB+1),:] #Amend young age cluster rows
		else #Both ages
			YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
			ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
			VaccUptakeBySeasonArray[1:(YoungClusterUB+1),:] = ExpandedProgVaccUptakeArray[1:(YoungClusterUB+1),:] #Amend young age cluster rows
			VaccUptakeBySeasonArray[(ElderClusterLB+1):end,:] = ExpandedProgVaccUptakeArray[(ElderClusterLB+1):end,:] #Amend elder age cluster rows
		end

		#Assign populated array to output tuple
		VaccUptakeBySeason[ii+1] = VaccUptakeBySeasonArray  #Add 1 to idx as first entry is for 2009/2010 season
	end
return VaccUptakeBySeason

end

@everywhere function VaccUptake_ExpandSchoolAndElderly_RecovRisk(SeasonsToSimulate,JobIdx,VaccAgeBoundsData,SimnRunType)

	#Initialise tuple per season
	#Per tupel, 2D array
	#rows for age (0 to 100+), cols for calendar day of year
	VaccUptakeBySeason = Vector{Array}(undef,SeasonsToSimulate)

	# Import the data - Pandemic flu vacc (2009/2010 season)
	#Load data based on risk group type
	if SimnRunType == 2
		#AT RISK
		PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc.xlsx","At risk","C3:NC103")
	elseif SimnRunType == 3
		#LOW RISK
		PandemicFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc.xlsx","Low risk","C3:NC103")
	end
	# Collate into Array, Assign to storage cell
	PandemicFluVaccUptakeArray = Array{Float64, 2}(PandemicFluVaccUptake)
	#PandemicFluVaccUptakeArray[19:end,:] .= 0 #Set specified age groups to have no vaccination
	VaccUptakeBySeason[1] = PandemicFluVaccUptakeArray

	#Import the data - Sesonal flu vacc
	SheetNames = ["2010_2011","2011_2012","2012_2013","2013_2014",
		"2014_2015","2015_2016","2016_2017","2017_2018"]

	#Iterate over each season
	for ii = 1:length(SheetNames)

		#Load data based on risk group type
		if SimnRunType == 2
		    #AT RISK
				#Load at-risk only vaccine programme data
				AtRiskProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_AtRisk_HistoricalAtRiskOnly.xlsx","$(SheetNames[ii])","C3:NC103")

		        #Load vaccine uptake data with additinal age groups covered
				ExpandedProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_AtRisk_AtRiskPlusElderly&0-20Schedule.xlsx","$(SheetNames[ii])","C3:NC103")
		elseif SimnRunType == 3
		    #LOW RISK
			#Load at-risk only vaccine programme data
			AtRiskProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_LowRisk_HistoricalAtRiskOnly.xlsx","$(SheetNames[ii])","C3:NC103")

			#Load vaccine uptake data with additinal age groups covered
			ExpandedProgSeasonalFluVaccUptake = XLSX.readdata("../../Data/VaccUptake/ByYrOfAge_DailyVaccUptakeCalYr_LowRisk_AtRiskPlusElderly&0-20Schedule.xlsx","$(SheetNames[ii])","C3:NC103")
	end

		# Collate input data into Arrays
		AtRiskVaccUptakeArray = Array{Float64,2}(AtRiskProgSeasonalFluVaccUptake)
		ExpandedProgVaccUptakeArray = Array{Float64,2}(ExpandedProgSeasonalFluVaccUptake)

		#Assign rows based on JobIdx (read from text file)
		CurrentJobIdx_VaccAgeBoundsData = VaccAgeBoundsData[JobIdx,:]
		YoungClusterUB = VaccAgeBoundsData[JobIdx,1] #First entry the upper bound for young age cluster
		ElderClusterLB = VaccAgeBoundsData[JobIdx,2]   #Second entry for lower bound of elder age cluster

		#Populate vaccine uptake array for current season iteration
		VaccUptakeBySeasonArray = AtRiskVaccUptakeArray #Initial data. Match at risk only programme

		if isnan(YoungClusterUB)  #Elder cluster only
			ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
			VaccUptakeBySeasonArray[(ElderClusterLB+1):end,:] = ExpandedProgVaccUptakeArray[(ElderClusterLB+1):end,:] #Amend elder age cluster rows
		elseif isnan(ElderClusterLB) #Young ages only
			YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
			VaccUptakeBySeasonArray[1:(YoungClusterUB+1),:] = ExpandedProgVaccUptakeArray[1:(YoungClusterUB+1),:] #Amend young age cluster rows
		else #Both ages
			YoungClusterUB = convert(Int64,YoungClusterUB) #Convert to interger
			ElderClusterLB = convert(Int64,ElderClusterLB) #Convert to interger
			VaccUptakeBySeasonArray[1:(YoungClusterUB+1),:] = ExpandedProgVaccUptakeArray[1:(YoungClusterUB+1),:] #Amend young age cluster rows
			VaccUptakeBySeasonArray[(ElderClusterLB+1):end,:] = ExpandedProgVaccUptakeArray[(ElderClusterLB+1):end,:] #Amend elder age cluster rows
		end

		#Assign populated array to output tuple
		VaccUptakeBySeason[ii+1] = VaccUptakeBySeasonArray  #Add 1 to idx as first entry is for 2009/2010 season
	end
return VaccUptakeBySeason

end

#-------------------------------------------------------------------------------
# FUNCTION TO BUILD EXPOSURE HISTORY ARRAY
#-------------------------------------------------------------------------------
@everywhere function BuildExpHistArray(ExpHistVaccType,NumOfStrains,M,ExpHistArrayParams,AgeGrpSuscepTuple)
#PURPOSE: Construct exposure history array based on current parameter set
#Inputs:
# ExpHistVaccType - Flag variable to specify how susceptibility will be modified for
#                   vaccine-related exposure history classes
#                  vacc. groups is reduced (0 - Not reduced, 1 - Reduced (age agnostic), 2 - Reduced (age dependent)
#   NumOfStrains - As titled
#   M - Number of age classes
#   ExpHistArrayParams - File containing parameter sets to be run
# AgeGrpSuscepTuple - Cell: Three entries
#                   --> Cell 1 - AgeGrpSuscepParamNum
#                   --> Cell 2 - AgeSuscepLowerBounds
#                   --> Cell 3 - AgeSuscepUpperBounds

#Outputs:
#   ExpHistArray - interaction array between exposure history and susceptibility to
#   the current season strain variant

    #Number of exposure history classes
    ExpHistNum = (NumOfStrains*2) + 2

    #Column per exposure history. Row per strain.
    ExpHistArray = ones(NumOfStrains,ExpHistNum,M)
    HalfExpHistNum = convert(Int64,ExpHistNum/2)

    #Column 1 - no exposure in previous season
    #Columns 2-5 - natural infection by one of the strains
    #Column 6  - vacc. in previous season, no natural infection
    #Columns 7-10 - Natural infection by ones of the sratins AND vaccinated

	if ExpHistVaccType == 1 #Age agnostic exposure history parameters

	    #Right hand columns, vaccinated previous season
	    ExpHistArray[:,HalfExpHistNum+1:end,:] .= ExpHistArrayParams[3]

	    for ii = 2:HalfExpHistNum
	        #Update entries: natural infection
	        ExpHistArray[ii-1,ii,:] .= ExpHistArrayParams[1]
	        ExpHistArray[ii-1,ii+HalfExpHistNum,:] .= ExpHistArrayParams[1]
	    end

	    #Update entries: influenza B corss-reactivity
	    ExpHistArray[NumOfStrains-1,[HalfExpHistNum end],:] .= ExpHistArrayParams[2]
	    ExpHistArray[NumOfStrains,[HalfExpHistNum-1 end-1],:] .= ExpHistArrayParams[2]
	elseif ExpHistVaccType == 2 #Age-dependent exposure history parameters

		#Disaggregate AgeGrpSuscepTuple
		AgeGrpSuscepParamNum = AgeGrpSuscepTuple[1]::Int64
		AgeSuscepLowerBounds = AgeGrpSuscepTuple[2]::Array{Int64,1}
		AgeSuscepUpperBounds = AgeGrpSuscepTuple[3]::Array{Int64,1}

		#Iterate through each age band (groupings match those used for susceptibility parameters)
		for jj = 1:AgeGrpSuscepParamNum

			#Array indexing variables
			StartIdx = AgeSuscepLowerBounds[jj] + 1
			EndIdx = AgeSuscepUpperBounds[jj] + 1

			#Right hand columns, vaccinated previous season
			ExpHistArray[:,HalfExpHistNum+1:end,StartIdx:EndIdx] .= ExpHistArrayParams[jj,3]

			for ii = 2:HalfExpHistNum
				#Update entries: natural infection
				ExpHistArray[ii-1,ii,StartIdx:EndIdx] .= ExpHistArrayParams[jj,1]
				ExpHistArray[ii-1,ii+HalfExpHistNum,StartIdx:EndIdx] .= ExpHistArrayParams[jj,1]
			end

			#Update entries: influenza B corss-reactivity
			ExpHistArray[NumOfStrains-1,[HalfExpHistNum end],StartIdx:EndIdx] .= ExpHistArrayParams[jj,2]
			ExpHistArray[NumOfStrains,[HalfExpHistNum-1 end-1],StartIdx:EndIdx] .= ExpHistArrayParams[jj,2]
		end
	end
    return ExpHistArray
end

#-------------------------------------------------------------------------------
# FUNCTIONS TO CONSTRUCT SUSCEPTIBILITY ARRAY
#-------------------------------------------------------------------------------

#TEMPLATE FOR ALL FUNCTIONS
#function  AgeSuscepTypeFn(x,AgeGrpSuscepParamNum,TransmissExpHistParamNum,NumOfStrains)

#	return AgeSuscepByGroup, AgeSuscepParamNum
#end

#Inputs:
#   x - (vector) Current paramter set
#   AgeGrpSuscepParamNum - (integer)
#   TransmissExpHistParamNum - (integer)
#   NumOfStrains - (integer)

#Outputs:
#   AgeSuscepByGroup - (2D array) Row by age group. Column by strain.
#	AgeSuscepParamNum - (scalar, integer)

#Age group speicifc only (independent of strain)
@everywhere function AgeOnlySuscepFn(x,AgeGrpSuscepParamNum,TransmissExpHistParamNum,NumOfStrains)

	AgeSuscepParamNum = AgeGrpSuscepParamNum::Int64
	AgeSuscepParamsOnly = x[TransmissExpHistParamNum+1:TransmissExpHistParamNum+AgeGrpSuscepParamNum]::Array{Float64,1}

	#Construct AgeSuscepByGroup from AgeSuscepParamsOnly
	#Row per susceptibility grouping, column per strain
	AgeSuscepByGroup = repeat(AgeSuscepParamsOnly,outer = (1,NumOfStrains)) #Populate each column with age-specific suscep.

	return AgeSuscepByGroup, AgeSuscepParamNum
end


#-------------------------------------------------------------------------------
# FUNCTIONS TO APPLY ASCERTAINMENT PROBABILITY TO INFECTION CASE LOAD
#-------------------------------------------------------------------------------

#TEMPLATE FOR ALL FUNCTIONS
#function  AscertainProbFn(AscertainProb,SeasonRatePerAgeStrain,
# 										AgeBandedCaseRateTotalSum_FitCheckSeasons,
# 										FitAggAgeFlag,AgeBandNum,
# 										AgeBandBounds,M)

#	return SeasonRatePerAgeStrain
#end

#Inputs:
#   AscertainProb - Ascertainment related parameters
#   SeasonRatePerAgeStrain - (3D array) Row per season. Column per age. Slice per age. Ascertained cases per 100,000.
#   AgeBandedCaseRateTotalSum_FitCheckSeasons - (3D array) Model simulated values for proprotion of age class infected. Row per season, column per age, slice per strain.
#   FitAggAgeFlag - (indicator flag) 0 for fitting single yr age groups. 1 for aggregated age groups.
#   AgeBandNum - (integer) Number of age classes in use for strucuring ascertainment profile.
#   AgeBandBounds - (2D array) First row for lower bounds. Second row for upper bounds.
# 	M - (interger) Total number of single year age classes in use
#Outputs:
#   SeasonRatePerAgeStrain - (3D array) Populated array


# Seasons specific baseline, with age group scaling. Piecewise linear function.
@everywhere function  SeasonWithAgeScaleLinearAscertainmentFn(AscertainProb,SeasonRatePerAgeStrain,
										AgeBandedCaseRateTotalSum_FitCheckSeasons,
										FitAggAgeFlag,AgeBandNum,
										AgeBandBounds,M)

	#Error check on parameter vector length
	RetainedSeasonNum = size(SeasonRatePerAgeStrain,1)
	if length(AscertainProb) != (RetainedSeasonNum+AgeBandNum)
		error("Number of parameters allocated to AscertainProb is $(length(AscertainProb)). Should be $(RetainedSeasonNum+AgeBandNum).")
	end

	#Note, only invoked if fitting to single year age classes
	#Fitting to binned data, throws error message & ends programme
	if FitAggAgeFlag == 0

		### First, scale values by max seasonal ascertainment probability (attained by those aged 100yrs) ###
		MaxAscertainProbBySeason = AscertainProb[1:RetainedSeasonNum]
		FluToGPconversionBySeason = 1*MaxAscertainProbBySeason #Multiply by 1, keep as a proportion of population unit!
		for jj = 1:RetainedSeasonNum
			SeasonRatePerAgeStrain[jj,:,:] = AgeBandedCaseRateTotalSum_FitCheckSeasons[jj,:,:].*FluToGPconversionBySeason[jj]
		end

		### Second, scale age groups with age-specific ascertainment modification ###

		#Assign ascertainment age band end value scaling factors to variable
		AscertainProbAgeScale = AscertainProb[RetainedSeasonNum+1:end]

		#End endpoint value of 1 (max age has highest ascertainment, unscaled)
		AscertainProbAgeScaleIncEndPoint = [AscertainProbAgeScale;1]

		#Disaggregate AgeBandBounds
		AgeBandLowerBounds = AgeBandBounds[1,:]::Array{Int64,1}
		AgeBandUpperBounds = AgeBandBounds[2,:]::Array{Int64,1}

		#Initialise array - Ascertain prob. by yr of age
		AscertainProbScaleByAge = zeros(M)

		#Iterate over each year of age ii
		for ii = 1:M

			#Convert loop index to year of age
			YrOfAgeVal = ii - 1

		  	#Check upper bound, first <= than.
		  	AgeBandVectorIdx = findfirst(YrOfAgeVal .<= AgeBandUpperBounds)

		  	#Identify values at edge of bins, Entry idx and idx+1.
		  	LeftEdgeVal = AscertainProbAgeScaleIncEndPoint[AgeBandVectorIdx]
		  	RightEdgeVal = AscertainProbAgeScaleIncEndPoint[AgeBandVectorIdx+1]

		  	#For age ii, ascertain scale value is
		  	#AscertainScale[idx] + (IdxInBin/NumYrsInBin)*(AscertainScale[idx+1] - AscertainScale[idx])

		  	#AgeBandSpan = (AgeBandUpperBounds[AgeBandVectorIdx] - AgeBandLowerBounds[AgeBandVectorIdx]) + 1 #Add 1 to reach start age of next age bin!
			if AgeBandVectorIdx == length(AgeBandLowerBounds) #FOr last age category, get age years up to and including 100
				AgeBandSpan = 100 - AgeBandLowerBounds[AgeBandVectorIdx]
			else
			  AgeBandSpan = AgeBandLowerBounds[AgeBandVectorIdx+1] - AgeBandLowerBounds[AgeBandVectorIdx]
			end

		  	IdxInAgeBand = (ii-1) - AgeBandLowerBounds[AgeBandVectorIdx] #Subtract 1 as ages begin from zero!
		  	AscertainProbScaleByAge[ii] = LeftEdgeVal + ((IdxInAgeBand/AgeBandSpan)*(RightEdgeVal - LeftEdgeVal))

			#Error check
			if AgeBandSpan == 0
				error("An age band contains only a single year of age. Incompatible with the selected piecewise linear ascertainment function.")
			end

		end

		#Apply individual age ascertainment scaling factor to previously calculated rate assuming max ascertainment (value for 100 yr old)
		for jj = 1:length(AscertainProbScaleByAge)
			SeasonRatePerAgeStrain[:,jj,:] = SeasonRatePerAgeStrain[:,jj,:].*AscertainProbScaleByAge[jj]
		end

	else
		error("Banded age data (FitAggAgeFlag set to 1) not compatible with choice of ascertainment function.")
	end

	return SeasonRatePerAgeStrain
end

#-------------------------------------------------------------------------------
# FUNCTION TO PASS INPUTS TO ODE SOLVER AND PROCESS OUTPUTS
#-------------------------------------------------------------------------------


# Replicate FMAlt_NonCalibVersion.
@everywhere function  RunModelFMAlt_RecoverRiskGrp(FixedModelParams,x,CurrentParamSetPopnFOI)
#Inputs:
#   FixedModelParams - Parameters with consistent values across all runs
#   x - Values of parameters undergoing consideration
#	CurrentParamSetPopnFOI - Force of infection at the population level (tied to current parameter set undergoing consideration)
#Outputs:
#   SimnData - Model output to be fed into summary statistic function

#--------------------------------------------------------------------------
### Outline of steps
#--------------------------------------------------------------------------
# (i) Assign FixedModelParams to variables
# (ii) For parameters being inferred, update associated parameters
# (iii) Run ODE model and get end of seasons case counts
# (iv) For each age band, sum across relevant ages
# (v) Declare output variable names

	#----------------------------------------------------------------------
	### (i) Disaggregate FixedModelParams inputs
	#----------------------------------------------------------------------
	ContactArray = FixedModelParams[1]
	MortalityFile = FixedModelParams[2]
	ONSPopnDistEngland = FixedModelParams[3]
	SimnRunType = FixedModelParams[4]
	ExpHistVaccType = FixedModelParams[5]
	StoreFlag_PopnFOI = FixedModelParams[6]
	SimnParam = FixedModelParams[7]
	NumOfStrains = FixedModelParams[8]
	ExpHistNum = FixedModelParams[9]
	M = FixedModelParams[10] #Number of single year age classes
	InfectionParam = FixedModelParams[11]
	MultiSeasonImmPropn = FixedModelParams[12]
	VaccUptakeBySeason = FixedModelParams[13]
	LeakyTransFlag = FixedModelParams[14]
	LeakyVaccVarBySeason = FixedModelParams[15]
	InfPropn_StartOfSeason = FixedModelParams[16]
	ICFromFile = FixedModelParams[17]
	RetainedSeasonNum = FixedModelParams[18]
	FitAggAgeTuple = FixedModelParams[19]
	AgeGrpSuscepTuple = FixedModelParams[20]
    RiskGrpSpecificAgeDist = FixedModelParams[21]
	AscertainProbFn = FixedModelParams[22]

	#----------------------------------------------------------------------
	### (ii) For parameters being inferred, update associated parameters
	#----------------------------------------------------------------------

	### Get R0, beta updated later! ###
    R_0 = x[1:4]::Array{Float64,1}
    #gamma = InfectionParam[3]::Array{Float64,1} #rate of loss of infectiousness
	#beta = (gamma.*R_0)./spec_rad ##ecover strain transmission rates from contact structure & R_0
	InfectionParam[4] = R_0 #Assign basic reproduction number vector to tuple

	TransmissParamNum = length(R_0)

	### Disaggregate susceptibility related variables ###
    #Role in altering transmission & exposure history associated
    #quantities!
	AgeGrpSuscepParamNum = AgeGrpSuscepTuple[1]::Int64
    AgeSuscepLowerBounds = AgeGrpSuscepTuple[2]::Array{Int64,1}
    AgeSuscepUpperBounds = AgeGrpSuscepTuple[3]::Array{Int64,1}
	AgeSuscepTypeFn = AgeGrpSuscepTuple[4]

	### Update exposure history ###
    #Build exposure history array. Assign to variable
	if ExpHistVaccType == 1 #Related to previous season vaccine efficacy (age agnostic)
        ExpHistArrayParamNum = 3

        #Pick out range of particle set vector, x, corresponding to exposure
        #history parameters
        TransmissExpHistParamNum = TransmissParamNum + ExpHistArrayParamNum
        ExpHistArrayParams = x[TransmissParamNum+1:TransmissExpHistParamNum]::Array{Float64,1}

    elseif ExpHistVaccType == 2 #Related to previous season vaccine efficacy AND age band
        ExpHistArrayParamNum = AgeGrpSuscepParamNum*3

        #Pick out range of particle set vector, x, corresponding to exposure
        #history parameters
        TransmissExpHistParamNum = TransmissParamNum + ExpHistArrayParamNum
        ExpHistArrayParamsVector = x[TransmissParamNum+1:TransmissExpHistParamNum]::Array{Float64,1}

        #Reshape ExpHistArrayParamsVector into 2D array
        # -> Row per age grouping (with differing suscepibility between groups)
        # -> Column per exposure history parameter type
        #       - Col1: Suscep. following natural infection
        #       - Col2: Type B cross-reactivity carry over
        #       - Col3: Propn. prior season vaccine efficacy carry over
        ExpHistArrayParams = reshape(ExpHistArrayParamsVector,AgeGrpSuscepParamNum,3)
    end

    ExpHistArray = BuildExpHistArray(ExpHistVaccType,NumOfStrains,M,
										ExpHistArrayParams,AgeGrpSuscepTuple)
	ExpHistArrayFnInputs = [ExpHistArray,ExpHistArrayParams]

	### Populate AgeSuscep array ###
	AgeSuscepByGroup, AgeSuscepParamNum = AgeSuscepTypeFn(x,AgeGrpSuscepParamNum,TransmissExpHistParamNum,NumOfStrains)

    #Age Suscep by yr of age to be passed to flu model run function
    #Row per age band, column per strain
    AgeSuscep = zeros(M,NumOfStrains)

    #Populate AgeSuscep, iterate through each age band
    #Assign susceptibilities to relevant single year of age classess
    for ii = 1:AgeGrpSuscepParamNum
        AgeSuscepStartIdx = AgeSuscepLowerBounds[ii] + 1 #Ages begin at 0. Add 1 to align with array indexing
        AgeSuscepEndIdx = AgeSuscepUpperBounds[ii] + 1

		#Populate each row of AgeSuscep array
		#Assign values in AgeSuscepByGroup[ii,:] to each single age class within age bandinterval
		if ii == AgeGrpSuscepParamNum #If eldest age suscep category, sum to final column.
			NumSingleYrAgeClasses = M-AgeSuscepStartIdx+1
        else
			NumSingleYrAgeClasses = AgeSuscepEndIdx-AgeSuscepStartIdx+1
        end

		for jj = 0:(NumSingleYrAgeClasses-1)
			AgeSuscep[AgeSuscepStartIdx+jj,:] = AgeSuscepByGroup[ii,:]
		end
    end

	#----------------------------------------------------------------------
    ###(iii) Run ODE model and get end of seasons case counts
    #----------------------------------------------------------------------

	#Run the model!
	StoreArrays = RunSeasonalFluModelFMAlt_RecoverRiskGrp(ContactArray,MortalityFile,ONSPopnDistEngland,RiskGrpSpecificAgeDist,
											SimnParam,InfectionParam,CurrentParamSetPopnFOI,AgeSuscep,AgeGrpSuscepTuple,
											VaccUptakeBySeason,LeakyVaccVarBySeason,LeakyTransFlag,
											ExpHistArrayFnInputs,ExpHistNum,ExpHistVaccType,MultiSeasonImmPropn,
											InfPropn_StartOfSeason,ICFromFile,SimnRunType)

	#Get timestep value & indicator value StoreFlag_PopnFOI
	timestep = SimnParam[8]

	#Disaggregate StoreArrays
	T = StoreArrays[1]::Array{Float64,1}
	C = StoreArrays[10]::Array{Float64,3}
	E_NotV = StoreArrays[4]::Array{Float64,3}
	E_V = StoreArrays[5]::Array{Float64,3}
	I_NotV = StoreArrays[6]::Array{Float64,3}
	I_V = StoreArrays[7]::Array{Float64,3}
	PopnDist = StoreArrays[11]::Array{Float64,2}

	#Visits in week (or month) t, denoted c_{t}, difference in the cumulative proportion
	#of consultations over the previous timestep
	#i.e. c_{t} = p(C(t)-C(t-1)),
	SliceVarIdx = convert(Int64,(365/timestep)+1) #Increment to access arrays at to obtain end of influenza season values
	CumulCaseCount = cat(C[1:SliceVarIdx:end,:,:],C[end:end,:,:]; dims=1) #Get cumul. case count at specified intervals
														#Concatenate over 1st dimension, add final cumulative case count
														#C[end:end,:,:], notation end:end used to return a three dimensional array
	PopnDistBySeason = PopnDist[SliceVarIdx:SliceVarIdx:end,:] #cat(1,PopnDist[366:366:end,:],PopnDist[end:end,:])
												#Get popn. distribution at specified intervals
														#Concatenate over 1st dimension, add final population distribution

	UnNormCaseCountPerSeason_YrAgeGrps = CumulCaseCount[2:end,:,:]-CumulCaseCount[1:end-1,:,:]

	#----------------------------------------------------------------------
    ### (iv) Get totals for desired age bands
    #----------------------------------------------------------------------
	FitAggAgeFlag = FitAggAgeTuple[1]::Int64 #Diasaggregate FitAggAgeTuple
    AgeBandBounds = FitAggAgeTuple[2]::Array{Int64,2}
	AgeBandNum = convert(Int64,size(AgeBandBounds,2)) #Number of columns of age bound array equiv. to number of age bands in use
    AgeBandLowerBounds = AgeBandBounds[1,:]  #isaggregate AgeBandBounds array into lower and upper bounds
    AgeBandUpperBounds = AgeBandBounds[2,:]

	#Normalise case count by propotion of popn contained in each age band
	if FitAggAgeFlag == 0  #Use single year age classes
		StrainSeasonRateTotalSum_YrAgeGrps = UnNormCaseCountPerSeason_YrAgeGrps./PopnDistBySeason #Per single yr age class, produce estimated flu+ GP visit per 100,000 popn
		AgeBandedCaseRateTotalSum_FitCheckSeasons = StrainSeasonRateTotalSum_YrAgeGrps[4:end,:,:] #Get model simulated values for 2012/2013 onward
	elseif FitAggAgeFlag == 1  #Use intervals/age buckets designated by AgeBandBounds
		AgeBandedCaseRateTotalSum_FitCheckSeasons = zeros(RetainedSeasonNum,AgeBandNum,NumOfStrains)
	    for ii = 1:AgeBandNum
	        StartIdx = AgeBandLowerBounds[ii] + 1 #Ages begin at 0. Add 1 to align with array indexing
	        EndIdx = AgeBandUpperBounds[ii] + 1

	        #For each flu season, sum across ages within age band and store
	        #Flu season per row (dimension 1), so for each row & strain slice sum across
	        #stated range of columns
	        if ii == AgeBandNum #If eldest age category, sum to final column.

				AgeBandInfAsPropnOvPopn = sum(UnNormCaseCountPerSeason_YrAgeGrps[:,StartIdx:end,:],2)
				AgeBandPopnAsPropnOvPopn = sum(PopnDistBySeason[:,StartIdx:end,:],2)


				#AgeBandedCaseRateTotalSum_FitCheckSeasons[:,ii,:] =
	             #   sum(StrainSeasonRateTotalSum_YrAgeGrps_FitCheckSeasons[:,StartIdx:end,:],2)
	        else #If not eldest age category, sum up to EndIdx.
				AgeBandInfAsPropnOvPopn = sum(UnNormCaseCountPerSeason_YrAgeGrps[:,StartIdx:EndIdx,:],2)
				AgeBandPopnAsPropnOvPopn = sum(PopnDistBySeason[:,StartIdx:EndIdx,:],2)
	        end

			#Carry out normalisation
			StrainSeasonRateTotalSum_AgeBandGrps = AgeBandInfAsPropnOvPopn./AgeBandPopnAsPropnOvPopn
			AgeBandedCaseRateTotalSum_FitCheckSeasons[:,ii,:] = StrainSeasonRateTotalSum_AgeBandGrps[4:end,:,:]
	    end
	end


	#----------------------------------------------------------------------
	### (v) Get ascertainable cases, based on ascertainment prob type
	#----------------------------------------------------------------------
	AscertainParamBeginIdx = TransmissExpHistParamNum+AgeSuscepParamNum+1
    AscertainProb = x[AscertainParamBeginIdx:end] #Proposed ascertainment probabilty values to test

	#Initialise SeasonRatePerAgeStrain storage array
    if FitAggAgeFlag == 0  #Use single year age classes
        SeasonRatePerAgeStrain = zeros(RetainedSeasonNum,M,NumOfStrains)
    elseif FitAggAgeFlag == 1 #Use intervals/age buckets designated by AgeBandBounds
        SeasonRatePerAgeStrain = zeros(RetainedSeasonNum,AgeBandNum,NumOfStrains)
    end

	#Apply ascertainment function
	AscertainedCases = AscertainProbFn(AscertainProb,SeasonRatePerAgeStrain,
											AgeBandedCaseRateTotalSum_FitCheckSeasons,
											FitAggAgeFlag,AgeBandNum,
											AgeBandBounds,M)

	#----------------------------------------------------------------------
	### (vi) Declare output variable names
	#----------------------------------------------------------------------
	#Compute infected temporal profile
    #Stratified by age group, but aggregated over strains
	Total_Infected_AllDimn = sum(I_NotV,dims=3) + sum(I_V,dims=3) + sum(E_NotV,dims=3) + sum(E_V,dims=3) #Note, returns a three dimensional array, n x n x 1!

	#Remove singleton dimensions from Total_I_AllDimn
	Total_Infected = dropdims(Total_Infected_AllDimn; dims=3)

	#Assign outputs that may be used for health economic evaluation
	#Output 1: Proportion of specified age bandings that were infected
	#Output 2: Proportion of specified age bandings that were ascertained cases
	#Output 3: Timeseries of all infected
	#Output 4: Population distribution info
	#SimnData = [AgeBandedCaseRateTotalSum_FitCheckSeasons,PopnDist]

	SimnData = [AgeBandedCaseRateTotalSum_FitCheckSeasons,AscertainedCases]
	#SimnData = [AgeBandedCaseRateTotalSum_FitCheckSeasons,AscertainedCases,Total_Infected,PopnDist]

	return SimnData::Array{Array{Float64,3},1}
end


@everywhere function RunSeasonalFluModelFMAlt_MultiSimns_RiskGrpVers(SimnRunType,timestep,popnFOIarray,ExpHistVaccType,
											FitAggAgeTuple,AgeGrpSuscepTuple,AscertainProbFn,
                                            ParticleSets,TotalRunNum,OutputFile,
											VaccUptakeFunction,JobIdx,VaccAgeBoundsFileName)
#Inputs:
#	SimnRunType - 1 - exploratory; 2 - At risk; 3 - Low risk (NFLUENCES VACCINE UPTAKE/EFFICACY)
#	timestep - Gives the temporal resolution of outputs from the ODE scheme
# 	popnFOIarray - Force of infection at the population level. Entry per simulation set
#   ExpHistVaccType - See description given before fn call
#   FitAggAgeTuple - See description given before fn call
#   AgeGrpSuscepTuple - See description given before fn call
#	AscertainProbFn - See description given before fn call
#   ParticleSets - Parameter sets to be simulated
#   TotalRunNum - Number of simulations, with unique parameter sets, to be
#                   performed
#   OutputFile - Save location (MAT object)
#   VaccUptakeFunction - (function) Module to load required vaccine uptake data
#   JobIdx - (scalar, integer) An ID used for modify the vaccine uptake data
#   VaccAgeBoundsFileName - (string) File containing age bounds used in each strategy

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

    SeasonsToSimulate = 9

    ODEBurnIn = 0*365
    ODEStaticPopnTime = 1*365
    ODEInferenceTime = (SeasonsToSimulate-1)*365
    ODEForwardSimnTime = 0*365


    #Get age distrubtion data within at-risk group
    RiskGrpSpecificAgeDist = readdlm("AgeDistWithinRiskGrpWork/AtRiskPopn_CountPerAge.txt",',')
elseif SimnRunType == 3

	SeasonsToSimulate = 9

	ODEBurnIn = 0*365
    ODEStaticPopnTime = 1*365
    ODEInferenceTime = (SeasonsToSimulate-1)*365
    ODEForwardSimnTime = 0*365

    #Get age distrubtion data within low risk group
    RiskGrpSpecificAgeDist = readdlm("AgeDistWithinRiskGrpWork/LowRiskPopn_CountPerAge.txt",',')
else
    error("Incorrect RunType entered")
end

MaxTime = ODEBurnIn + ODEStaticPopnTime + ODEInferenceTime + ODEForwardSimnTime

#Compute total time post burn-in
ODESampleTime = ODEStaticPopnTime + ODEInferenceTime + ODEForwardSimnTime

#Concatenate simulation variables
SimnParam=[SimnStartDate,ODEBurnIn,ODEStaticPopnTime,ODEInferenceTime,ODEForwardSimnTime,MaxTime,ODESampleTime,timestep,StoreFlag_PopnFOI]
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

#------------------------------------------------------------------------------
###  VACCINATION UPTAKE
#--------------------------------------------------------------------------
VaccAgeBoundsData = readdlm(VaccAgeBoundsFileName,',')
@time VaccUptakeBySeason = VaccUptakeFunction(SeasonsToSimulate,JobIdx,VaccAgeBoundsData,SimnRunType)

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
RetainedSeasonNum = 6

#--------------------------------------------------------------------------
# AGGREGATE FIXED PARAMETERS
#--------------------------------------------------------------------------
FixedModelParams = [ContactArray,MortalityFile,ONSPopnDistEngland20102018,
                    SimnRunType,ExpHistVaccType,StoreFlag_PopnFOI,
                    SimnParam,NumOfStrains,ExpHistNum,M,InfectionParam,
                    MultiSeasonImmPropn,VaccUptakeBySeason,LeakyTransFlag,
                    LeakyVaccVarBySeason,InfPropn_StartOfSeason,ICFromFile,
                    RetainedSeasonNum,FitAggAgeTuple,AgeGrpSuscepTuple,RiskGrpSpecificAgeDist,
					AscertainProbFn]

#--------------------------------------------------------------------------
### INITIALISE CELL TO STORE OUTPUTS OF INTEREST
#--------------------------------------------------------------------------
OutputSimnCell = Array{Array{Float64},2}(undef,TotalRunNum,2)

for ii = 1:TotalRunNum

    #----------------------------------------------------------------------
    ### Read in parameter values from file
    #----------------------------------------------------------------------
    ParticleSetVal = ParticleSets[ii,:]
    println("ParticleSetVal: $ParticleSetVal")

	  CurrentParamSetPopnFOI = popnFOIarray[ii]
	  #CurrentParamSetPopnFOI = popnFOIarray

    #----------------------------------------------------------------------
    ### CALL FUNCTION TO RUN MODEL
    #----------------------------------------------------------------------
    OutputSimnCell[ii,:] = RunModelFMAlt_RecoverRiskGrp(FixedModelParams,ParticleSetVal,CurrentParamSetPopnFOI)

    println("Run $ii complete")

end

#----------------------------------------------------------------------
### SAVE TO OUTPUT VARIABLE
#----------------------------------------------------------------------
println("JobIdx: $JobIdx")
write(OutputFile, "SimnData_JobID$JobIdx", OutputSimnCell) #Label output with Job ID

end

function FMAlt_MultiSimnRecoverRiskVers_RunVaccAgeSweep(ARGS)

#--------------------------------------------------------------------------
### ADD FILES TO SEARCH PATH FOR ODES/MODEL RUN FUNCTION
@everywhere include("RunSeasonalFluModelFMAlt_RecRiskFns_Julia.jl")
@everywhere include("../ModelFns/ExpHistUpdateJulV1.jl")

#--------------------------------------------------------------------------
### Set ID for overall population data input file
#--------------------------------------------------------------------------
OverallPopnDataFileID = ARGS[1]

#-------------------------------------------------------------------------------
# SPECIFY TYPE OF RUN THROUGH FLAG VARIABLE
#-------------------------------------------------------------------------------
# (INFLUENCES VACCINE UPTAKE/EFFICACY)
RiskGroupFlag = parse(Int64, ARGS[2])
if RiskGroupFlag == 0
	SimnRunType = 3
elseif RiskGroupFlag == 1
	SimnRunType = 2
end
#-------------------------------------------------------------------------------
# DECLARE REGULARITY OF OUTPUT FROM ODE SCHEME
#-------------------------------------------------------------------------------
timestep = 1.

#--------------------------------------------------------------------------
#Select exposure history susceptibility modification form for vaccine-related states
#0 - Fixed values every season
#1 - Relate to previous season vaccine efficacy
#2 - Relate to previous season vaccine efficacy AND age band
#--------------------------------------------------------------------------
ExpHistVaccType = parse(Int64, ARGS[3])

if ExpHistVaccType != 0 && ExpHistVaccType !=1 && ExpHistVaccType !=2
    error("ExpHistVaccType must take value 0, 1 or 2. Current value is $(ExpHistVaccType)")
end

#--------------------------------------------------------------------------
# Declare bounds for age category groupings (for aggregating case counts &
# applying ascertainment probabilitiy, if required)
#--------------------------------------------------------------------------
FitAggAgeFlag = parse(Int64, ARGS[4]) #Indicator variable.
# Set to: -> 1 to fit to aggregetated age band counts
#         -> 0 to fit to single yr age classes

#Error check
if FitAggAgeFlag != 0 && FitAggAgeFlag !=1
   error("FitAggAgeFlag must take value 0 or 1. Current value is $FitAggAgeFlag")
end

#Use Meta.parse for evaluating a string into a vector
AgeBandLowerBounds = eval(Meta.parse(ARGS[5]))
AgeBandUpperBounds = eval(Meta.parse(ARGS[6]))

#Concantenate bounds into a single array
AgeBandBounds = [AgeBandLowerBounds;AgeBandUpperBounds]
                #vertically stack age category lower bounds above upper bounds

#Place FitAggAgeFlag&AgeBandBounds into a cell/tuple
FitAggAgeTuple = [FitAggAgeFlag,AgeBandBounds]

#--------------------------------------------------------------------------
# Declare bounds for age susceptibility groupings
#--------------------------------------------------------------------------
AgeSuscepLowerBounds = eval(Meta.parse(ARGS[7]))
AgeSuscepUpperBounds = eval(Meta.parse(ARGS[8]))
AgeGrpSuscepParamNum = length(AgeSuscepLowerBounds)

#Select age&strain susceptibility modification form
# AgeAndStrainDepSuscepFn - Age group & strain specific
# AgeDepWithStrainScalingSuscepFn - Age group specific, with strain modifier
# StrainDepWithAgeScalingSuscepFn - Strain specific, with age group modifier
# AgeOnlySuscepFn - Age group speicifc only (independent of strain)
AgeSuscepTypeSymb = Symbol(ARGS[9]) #Convert string to Symbol
AgeSuscepType = getfield(Main, AgeSuscepTypeSymb)

#Concantenate into a cell tuple
AgeGrpSuscepTuple = [AgeGrpSuscepParamNum,AgeSuscepLowerBounds,AgeSuscepUpperBounds,AgeSuscepType]

#--------------------------------------------------------------------------
# Specify ascertainment function to be employed
#--------------------------------------------------------------------------
AscertainProbSymb = Symbol(ARGS[10]) #Convert string to Symbol
AscertainProbFn = getfield(Main, AscertainProbSymb)

#--------------------------------------------------------------------------
# Read in parameter values from file
#--------------------------------------------------------------------------
ParamInputFile = "../../../../Results/FMAlt/FMAlt_OptimFitOutput/FMAlt_OptimFitOutputFiles/OptimiserParamsTrace#14_Emp_TransContactMatrix_JobArrayCombined.txt";
# ParticleSets = readdlm(ParamInputFile,'\t')
ParticleSets = readdlm(ParamInputFile,',')

# ParticleSets = [2.01062 2.1363 1.95428 1.82481 0.667407 1.0 0.0 0.883729 0.488437 0.601442 0.976659 0.00764689 0.00563268 0.0130081 0.0187802 0.00920678 0.0379274 0.152626 0.0690459 0.174723 0.841427 0.5083;
#                     1.99096 2.04144 1.9072 1.82224 0.70022 1.0 0.0 0.791175 0.567555 0.731654 1.0 0.00655942 0.00458805 0.0104459 0.0141421 0.00961335 0.0296937 0.207495 0.107754 0.225128 0.839037 0.56751;
#                     2.05301 2.10649 1.95661 1.85392 0.727368 1.0 0.0 0.761494 0.596216 0.78163 1.0 0.00436862 0.00314509 0.00763382 0.0102824 0.00513015 0.0199153 0.265364 0.153105 0.296588 0.956122 0.676962;
#                     2.2107 2.28835 2.12554 2.02315 0.705386 1.0 0.0 0.70467 0.510179 0.719491 1.0 0.00662774 0.00447091 0.0107637 0.0158623 0.00817594 0.0286303 0.200995 0.108641 0.234815 0.789829 0.453041;
#                     2.07687 2.17054 1.78996 1.85545 0.753033 1.0 0.0 0.823868 0.611319 0.787901 1.0 0.00485143 0.00348207 0.00844311 0.0113242 0.00560439 0.0228705 0.171516 0.105677 0.233172 0.653591 0.457677]

#----------------------------------------------------------------------
### SPECIFY OUTPUT FILE SAVE LOCATION
#----------------------------------------------------------------------
SaveFileID = ARGS[11]

#----------------------------------------------------------------------
### SPECIFY VACCINE UPTAKE FUNCTION NAME
#----------------------------------------------------------------------
VaccUptakeSymb = Symbol(ARGS[12]) #Convert string to Symbol
VaccUptakeFunction = getfield(Main, VaccUptakeSymb) #Make Symbol callable functions

#--------------------------------------------------------------------------
#Call function to run multiple simulations
#--------------------------------------------------------------------------
JobNum = parse(Int64, ARGS[13])
TotalRunNumPerJob = parse(Int64, ARGS[14])
VaccAgeBoundsFileName = ARGS[15]

#Loop over each job
@sync @distributed for JobIdx = 1:JobNum

	#--------------------------------------------------------------------------
	### LOAD RELEVANT POPULATION LEVEL INFECTIOUS DATA
	MATfile = matopen("../../Results/FMAlt_SimnJuliaV1ModelRun_#$(OverallPopnDataFileID)_JobID#$(JobIdx).mat")
	SimnDataAllPopn = read(MATfile,"SimnData_JobID$JobIdx")
	close(MATfile)

	#Pick out portions of data file related to population force of infection measurements
	popnFOIarray = SimnDataAllPopn[:,4]

	#Set output filename
	if RiskGroupFlag == 0
		OutputFName = "../../Results/FMAlt_SimnJuliaV1ModelRun_#$(SaveFileID)_JobID#$(JobIdx)_LowRisk.mat"
	else RiskGroupFlag == 1
		OutputFName = "../../Results/FMAlt_SimnJuliaV1ModelRun_#$(SaveFileID)_JobID#$(JobIdx)_AtRisk.mat"
	end

	#Set up MAT file
	OutputFile = matopen(OutputFName,"w")

	#Run risk group simulation
	RunSeasonalFluModelFMAlt_MultiSimns_RiskGrpVers(SimnRunType,timestep,popnFOIarray,ExpHistVaccType,
									FitAggAgeTuple,AgeGrpSuscepTuple,AscertainProbFn,
                                    ParticleSets,TotalRunNumPerJob,OutputFile,
									VaccUptakeFunction,JobIdx,VaccAgeBoundsFileName)

	#Close MAT file
	close(OutputFile)
end



end

#--------------------------------------------------------------------------
# PASS TO FUNCTION
#--------------------------------------------------------------------------
FMAlt_MultiSimnRecoverRiskVers_RunVaccAgeSweep(ARGS)

for ii in workers()
    rmprocs(ii)
end
