#Purpose:
#Functions and script to calculate number of episodes, monetary costs & health costs
#based on predicted influenza case load

#Section overview:

# 1: Calculate counts of infected & GP visits
# 2: Set up health count output variables
# 3: Compute counts of all remaining clinical events
# 4: Discount rate defined
# 5: Compute monetary costs for all clinical event types (but not vaccination)
# 6: QALY loss calculation I (non-fatal clinical events)
# 7: QALY loss calculation II (fatal, in-hospital)
# 8: QALY loss calculation III (fatal, out of hospital)
# 9: QALY loss calculation IV (combined QALY losses)
# 10: Vaccination costs
# 11: Vaccines deployed tuple assignment
#--------------------------------------------------------------------------


function CombinedHealthEpsStat(ByYrofAge_CombinedHealthEpsStat,ByYrOfAge_IndHealthEpsStat,HealthEpsEntityNum)
#Inputs:
#ByYrofAge_CombinedHealthEpsStat - (Tuples)
#ByYrOfAge_IndHealthEpsStat - (scalar, integer) Number of health episodes to be iterated over
#HealthEpsEntityNum - (scalar, integer) Number of health episodes to be iterated over

#Outputs:
#ByYrOfAgeHealthEpsData,AllAgeHealthEpsData - (Tuples)

for HealthEpsIdx = 1:HealthEpsEntityNum
    #Contribution towards combined cost of all health episodes
    ByYrofAge_CombinedHealthEpsStat = ByYrofAge_CombinedHealthEpsStat + ByYrOfAge_IndHealthEpsStat[HealthEpsIdx]
end

return ByYrofAge_CombinedHealthEpsStat

end

function DoHealthEconEval(PropnInfData,PropnAscertainedData,GPVisitFromAscertainProbFlag,
                            PopnData,FebrileCaseParam,GPFluAscertain,
                           HealthEpsRelLhood,HealthEpsCosts,HealthEpsNonFatalQALYLoss,
                           WeightsQoL,RemainLifeExpAllSeasons,DiscountRate,VaccCostsData)
#Inputs from single simulation run.
#Outputs are strain-statified. Can be post-processed to get data pertaining to all influenza strains.

#Inputs:
#    PropnInfData - (3D array) Row per season, column per age, slice per strain.
#    PropnAscertainedData - (3D array) Row per season, column per age, slice per strain.
#    GPVisitFromAscertainProbFlag - (flag variable) Determines GP visits should be obtained from ascertainment cases or ALL infected cases
#    PopnData - (2D array) Population size. Row per season. Column per age.
#    FebrileCaseParam - (scalar) Febrile case parameter value (age and strain independent)
#    GPFluAscertain - (column vector) GP visit proportion. Entry per age.
#    HealthEpsRelLhood - (3D array)
#                       -> Row per single year age group,
#                       -> Column per episode type (non-fatal inpatient, outpatient, fatal case in-hospital, fatal case out of hospital)
#                       -> Slice per strain type (Slice one&two - Influenza A; Slice three&four - Influenza B)
#    HealthEpsCosts - (3D array)
#                       -> Row per single year age group,
#                       -> Column per episode type (GP consultation, non-fatal inpatient, outpatient, fatal case in-hospital, fatal case out of hospital)
#                       -> Slice per strain type (Slice one&two - Influenza A; Slice three&four - Influenza B)
#    HealthEpsNonFatalQALYLoss - (3D array)
#                       -> Row per single year age group,
#                       -> Column per non-fatal episode type (non-hospitalised, hospitalised)
#                       -> Slice per strain type (Slice one&two - Influenza A; Slice three&four - Influenza B)
#     FatalEventOutOfHospPropn - (vector)  Entry per year of age
#     WeightsQoL - (vector) Adjusted quality of life. Entry per year of age (0-100+)
#     RemainLifeExpAllSeasons - (2D array) Row per age. Column per reference year
#     DiscountRate - (two element vector,float) [CostDiscountRate,QALYDiscountRate]
#     VaccCostsData - (Tuple) Entry 1: Vaccine uptake ->Row per age, ->Column per season.
#                              Entry 2: Administration charge per vaccine dose
#                              Entry 3: Cost per vaccine dose

#Outputs:
#   ByYrOfAgeHealthEpsCounts, AllAgeHealthEpsCounts (cells, 7 entries)
#       - Counts by year of age AND overall: (i) GP visit; (ii) non-fatal hospital admissions, (iii) outpatient
#           visits, (iv) fatal cases (in hospital); (v) fatal cases (out of hosptial);
#           (vi) all infected cases (including asymptomatic); (vii) symptomatic cases only.
#   ByYrOfAgeHealthEpsCosts, AllAgeHealthEpsCounts (cells, 6 entries)
#       -   Monetary costs by year of age AND overall: (i) GP visit; (ii) non-fatal hospital admissions, (iii) outpatient
#           visits, (iv) fatal cases (in hospital only) (v) fatal cases (out of hospital)
#           (vi) Total cost by age (all treatments) and over all ages (all treatments)
#    ByYrOfAgeHealthEpsQALYloss ,AllAgeHealthEpsQALYloss (cells, 5 entries)
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
### SECTION 1: GET TOTAL INFECED CASE COUNT. COMPUTE GP CONSULTATIONS
#--------------------------------------------------------------------------
InfectedCaseCount = PropnInfData.*PopnData
AscertainedCaseCount = PropnAscertainedData.*PopnData

#Get relevant statistic values from dimensions of InfectedCaseCount
GPVisitCountTempDims = size(InfectedCaseCount)

#Check number of dimensions of InfectedCaseCount
#If 2, using non-age structured model
#If 3, using age-structured model
if ndims(InfectedCaseCount) == 2
    SeasonNum = GPVisitCountTempDims[1]
    ByYrOfAgeGrpNum = 1
    StrainNum = GPVisitCountTempDims[2]

    #Reshape so first dimension based on YrOfAge
    InfectedCaseCountReshaped = reshape(InfectedCaseCount,(ByYrOfAgeGrpNum,SeasonNum,StrainNum))
    AscertainedCaseCountReshaped = reshape(AscertainedCaseCount,(ByYrOfAgeGrpNum,SeasonNum,StrainNum))
else
    SeasonNum = GPVisitCountTempDims[1]
    ByYrOfAgeGrpNum = GPVisitCountTempDims[2]
    StrainNum = GPVisitCountTempDims[3]

    #Have a 3D array. Permute so first dimension based on YrOfAge
    InfectedCaseCountReshaped = permutedims(InfectedCaseCount,[2 1 3])
    AscertainedCaseCountReshaped = permutedims(AscertainedCaseCount,[2 1 3])
end

#Get GP visit count
#Note: NO ROUNDING IN USE!
if GPVisitFromAscertainProbFlag == 1 #Ascertained cases correspond to GP cases!
    GPVisitCount = AscertainedCaseCountReshaped
elseif GPVisitFromAscertainProbFlag == 2 #Ascertained cases correspond to GP cases! (risk group specific probability)
    GPFluAscertainReshaped = GPFluAscertain'  #Note, GPFluAscertain here is a 2D array, row by flu season, column by age
    GPVisitCount = InfectedCaseCountReshaped.*GPFluAscertainReshaped  #GPFluAscertainReshaped. Now row by age, column by flu season.
else #Multiply each infected case count by GP visit probability
    GPVisitCount = InfectedCaseCountReshaped.*GPFluAscertain
end

#--------------------------------------------------------------------------
### SECTION 2: AMALGAMATE HEALTH EPISODE COUNT OUTPUT VARIABLES
#--------------------------------------------------------------------------

#Output: Counts by year of age AND overall: (i) GP visit; (ii) non-fatal hospital admissions, (iii) outpatient
#visits, (iv) fatal cases (in hospital); (v) fatal cases (out of hosptial);
#(vi) all infected cases; (vii) symptomatic cases only.
ByYrOfAgeHealthEpsCounts = Array{Array{Float64}}(undef,7)
AllAgeHealthEpsCounts = Array{Array{Float64}}(undef,7)

#GP visits
ByYrOfAgeHealthEpsCounts[1] = GPVisitCount
AllAgeHealthEpsCounts[1] = sum(GPVisitCount,dims=1)

#Total infected cases (including asymptomatic)
ByYrOfAgeHealthEpsCounts[6] = InfectedCaseCountReshaped
AllAgeHealthEpsCounts[6] = sum(InfectedCaseCountReshaped,dims=1)

#Total infected cases (symptomatic only)
SympCaseCount = FebrileCaseParam*InfectedCaseCountReshaped;
ByYrOfAgeHealthEpsCounts[7] = SympCaseCount
AllAgeHealthEpsCounts[7] = sum(SympCaseCount,dims=1)

#--------------------------------------------------------------------------
### SECTION 3:
### COMPUTE NUMBER OF EACH HEALTH EPISODE TYPE AS COUNTS, NO DISCOUNTING PERFORMED
#--------------------------------------------------------------------------
#Input: Relative likelihood of health episode (relative to GP
#consultation).
#   - Array with row per single year age group, four columns (non-fatal
#   inpatient, outpatient, fatal case in-hospital, fatal case out of hospital).
# Slice per influenza type (Slice one&two - Influenza A; Slice three&four - Influenza B)

#Call function.
#Popualtes counts for non-fatal admissions, outpatient visits, in-hospital deaths
for HealthEpsIdx = 1:3
    OutputCellIdx = HealthEpsIdx + 1

    #For designated health episode, scale GP count using relative likelihood data
    ByYrOfAge_HealthEpsCount = GPVisitCount.*copy(selectdim(HealthEpsRelLhood,2,[HealthEpsIdx])) #Slice equivalent to HealthEpsRelLhood[:,HealthEpsIdx,:]
    #use of [] retains singleton dimension

    #Assign to output cell
    ByYrOfAgeHealthEpsCounts[OutputCellIdx] = ByYrOfAge_HealthEpsCount
    AllAgeHealthEpsCounts[OutputCellIdx] = sum(ByYrOfAge_HealthEpsCount,dims=1)
end

#Obtain counts for out-of-hospital deaths. Apply scaling factor to in-hospital deaths
InHospFatalCaseCount = ByYrOfAgeHealthEpsCounts[4]
FatalCasePropnOutOfHosp = copy(selectdim(HealthEpsRelLhood,2,[4])) #Slice equivalent to HealthEpsRelLhood[:,4,:]
OutOfHospFatalCaseCount = InHospFatalCaseCount.*FatalCasePropnOutOfHosp

#Assign to output tuple
OutOfHospFatalCaseCountTupleIdx = 5
ByYrOfAgeHealthEpsCounts[OutOfHospFatalCaseCountTupleIdx] = OutOfHospFatalCaseCount
AllAgeHealthEpsCounts[OutOfHospFatalCaseCountTupleIdx] = sum(OutOfHospFatalCaseCount,dims=1)

#--------------------------------------------------------------------------
### SECTION 4: DEFINE DISCOUNT RATE, GET DISCOUNT SCALING TERM FOR EACH TIMESTEP
#--------------------------------------------------------------------------

#Disaggregate discount rate
CostDiscountRate = DiscountRate[1]
QALYDiscountRate = DiscountRate[2]

TimeElapsedFromRefPoint = 0:1:SeasonNum-1 #Note initial season will have no discounting!
CostDiscountTerm = (1/(1+CostDiscountRate)).^TimeElapsedFromRefPoint
QALYDiscountTerm = (1/(1+QALYDiscountRate)).^TimeElapsedFromRefPoint

#--------------------------------------------------------------------------
### SECTION 5: COMPUTE MONETARY COST OF EACH HEALTH EPISODE TYPE
#--------------------------------------------------------------------------
#Input (HealthEpsCosts): Monetary cost of health episode
#   - Array with row per single year age group, five columns (GP consultation, non-fatal inpatient, outpatient, fatal case in-hospital, fatal case out of hospital)

#Output: Monetary costs by year of age AND overall: (i) GP visit; (ii) non-fatal hospital admissions, (iii) outpatient
#visits, (iv) fatal cases (in hospital only) (v) fatal cases (out of hospital)
#(vi) Total cost by age (all treatments) and over all ages (all treatments)
ByYrOfAgeHealthEpsCosts = Array{Array{Float64}}(undef,6)
AllAgeHealthEpsCosts = Array{Array{Float64}}(undef,6)

ByYrofAge_CombinedHealthEpsCost = zeros(ByYrOfAgeGrpNum,SeasonNum,StrainNum)
for HealthEpsIdx = 1:5
    #For designated health episode, scale number of events by cost per episode
    ByYrOfAge_HealthEpsCost = ByYrOfAgeHealthEpsCounts[HealthEpsIdx].*copy(selectdim(HealthEpsCosts,2,[HealthEpsIdx])) #.*HealthEpsCosts[:,HealthEpsIdx]

    #Carry out discounting in later seasons
    ByYrOfAge_HealthEpsCostDiscounted = zeros(ByYrOfAgeGrpNum,SeasonNum,StrainNum)
    for SeasonIdx = 1:SeasonNum
        ByYrOfAge_HealthEpsCostDiscounted[:,SeasonIdx,:] = ByYrOfAge_HealthEpsCost[:,SeasonIdx,:]*CostDiscountTerm[SeasonIdx]
    end

    #Assign to output cell
    ByYrOfAgeHealthEpsCosts[HealthEpsIdx] = ByYrOfAge_HealthEpsCostDiscounted
    AllAgeHealthEpsCosts[HealthEpsIdx] = sum(ByYrOfAge_HealthEpsCostDiscounted,dims=1)

    #Contribution towards combined cost of all health episodes
    ByYrofAge_CombinedHealthEpsCost = ByYrofAge_CombinedHealthEpsCost + ByYrOfAge_HealthEpsCostDiscounted
end

#Assign final output to Costs cell.
#Total cost by age (all treatments) and over all ages (all treatments)
ByYrOfAgeHealthEpsCosts[6] = ByYrofAge_CombinedHealthEpsCost
AllAgeHealthEpsCosts[6] = sum(ByYrofAge_CombinedHealthEpsCost,dims=1)

#--------------------------------------------------------------------------
### SECTION 6: QALY LOSS CALCULATION I (FOR NON-FATAL EVENTS)
#--------------------------------------------------------------------------
#Input (HealthEpsNonFatalQALYLoss): QALY losses (non-fatal)
#   - 3D Array. Row per age group. Two columns (GP consultation,
#   hospitalisation), Slice for influenza type (Slice one&two - Influenza A; Slice three&four - Influenza B)

#Output: QALY losses by year of age AND overall: (i) Symptomatic, Non-hospitalised; (ii) non-fatal hospital admissions,
#(iii) fatal cases (in hospital only), (iv) fatal cases (out of hospital)
#(v) Total cost by age (all treatments) and over all ages (all treatments)
ByYrOfAgeHealthEpsQALYloss = Array{Array{Float64}}(undef,5)
AllAgeHealthEpsQALYloss = Array{Array{Float64}}(undef,5)

#Get non-hospitalised events by subtracting from overall number of cases,
#cumulative sum of non-fatal hospital admissions (ByYrOfAgeHealthEpsCounts{2})
#and fatal cases (in-hospital & out-of-hospital, ByYrOfAgeHealthEpsCounts{4} & ByYrOfAgeHealthEpsCounts{5} respecitvely)
HospAndFatalCases = ByYrOfAgeHealthEpsCounts[2] + ByYrOfAgeHealthEpsCounts[4] + ByYrOfAgeHealthEpsCounts[5]

#Get symptomatic, non-hospitalised cases
println("FebrileCaseParam: $FebrileCaseParam")
SympNonHospCases = SympCaseCount - HospAndFatalCases

#Put into tuple non-fatal cases, straified by non-hospitalsied and
#hospitalsied
HospNonFatalCases = ByYrOfAgeHealthEpsCounts[2]
NonFatalCaseBreakdown = [SympNonHospCases,HospNonFatalCases]

#println("HealthEpsNonFatalQALYLoss: $(HealthEpsNonFatalQALYLoss[1:2,:,:])")
for NonFatalCaseTypeIdx = 1:2

    #For designated health episode, scale number of events by QALY loss per
    #episode
    ByYrOfAge_HealthEpsQALYloss = NonFatalCaseBreakdown[NonFatalCaseTypeIdx].*copy(selectdim(HealthEpsNonFatalQALYLoss,2,[NonFatalCaseTypeIdx])) #.*HealthEpsNonFatalQALYLoss[:,NonFatalCaseTypeIdx,:]

    #Carry out discounting in later seasons
    ByYrOfAge_HealthEpsQALYlossDiscounted = zeros(ByYrOfAgeGrpNum,SeasonNum,StrainNum);
    for SeasonIdx = 1:SeasonNum
        ByYrOfAge_HealthEpsQALYlossDiscounted[:,SeasonIdx,:] = ByYrOfAge_HealthEpsQALYloss[:,SeasonIdx,:]*QALYDiscountTerm[SeasonIdx]
    end

    #Assign to output cell
    ByYrOfAgeHealthEpsQALYloss[NonFatalCaseTypeIdx] = ByYrOfAge_HealthEpsQALYlossDiscounted
    AllAgeHealthEpsQALYloss[NonFatalCaseTypeIdx] = sum(ByYrOfAge_HealthEpsQALYlossDiscounted,dims=1)
end

#--------------------------------------------------------------------------
### SECTION 7: QALY LOSS CALCULATION II (FATAL, IN HOSPITAL)
#--------------------------------------------------------------------------

#Get death count data from Health Episode Count array
DeathCountsInHosp = ByYrOfAgeHealthEpsCounts[4]

#Iterate through each season, get combined fatal QALY losses
#Uses QoL weights and Expectancy of life estimates
CombinedFatalQALYlossPerSeason = zeros(1,SeasonNum,StrainNum)
ByYrOfAgeFatalQALYlossPerSeason = zeros(ByYrOfAgeGrpNum,SeasonNum,StrainNum)
for SeasonIdx = 1:SeasonNum

    #Specify SimnTime
    SimnTime = SeasonIdx - 1  #Reference point. First season simulated has value 0

    #Death counts occuring in SeasonIdx
    DeathCountsInHosp_BySeasonAndStrain = copy(selectdim(DeathCountsInHosp,2,[SeasonIdx])) #MATLAB equiv: DeathCountsInHosp[:,SeasonIdx,:]

    #Extract remaining life expectancy associated with SeasonIdx
    RemainLifeExpSingleSeason = RemainLifeExpAllSeasons[:,SeasonIdx]

    #Call function to calulate QALY losses resulting from fatal cases
    CombinedFatalQALYlossPerSeason[:,SeasonIdx,:],ByYrOfAgeFatalQALYlossPerSeason[:,SeasonIdx,:] = FatalQALYLoss(SimnTime,DeathCountsInHosp_BySeasonAndStrain,WeightsQoL,RemainLifeExpSingleSeason,QALYDiscountRate)
end

#--------------------------------------------------------------------------
### SECTION 8: QALY LOSS CALCULATION III (FATAL, OUT OF HOSPITAL)
#--------------------------------------------------------------------------

#Get death count data from Health Episode Count array
DeathCountsOutOfHosp = ByYrOfAgeHealthEpsCounts[5]

#Iterate through each season, get combined fatal QALY losses
Combined_OutOfHospFatalQALYlossPerSeason = zeros(1,SeasonNum,StrainNum)
ByYrOfAge_OutOfHospFatalQALYlossPerSeason = zeros(ByYrOfAgeGrpNum,SeasonNum,StrainNum)
for SeasonIdx = 1:SeasonNum

    #Speciy SimnTime
    SimnTime = SeasonIdx - 1 #Reference point. First season simulated has value 0

    #Death counts occuring in SeasonIdx
    DeathCountsOutOfHosp_BySeasonAndStrain = copy(selectdim(DeathCountsOutOfHosp,2,[SeasonIdx])) #MATLAB equiv: DeathCountsOutOfHosp[:,SeasonIdx,:]

    #Extract remaining life expectancy associated with SeasonIdx
    RemainLifeExpSingleSeason = RemainLifeExpAllSeasons[:,SeasonIdx]

    #Call function to calulate QALY losses resulting from fatal cases
    Combined_OutOfHospFatalQALYlossPerSeason[:,SeasonIdx,:],ByYrOfAge_OutOfHospFatalQALYlossPerSeason[:,SeasonIdx,:] = FatalQALYLoss(SimnTime,DeathCountsOutOfHosp_BySeasonAndStrain,WeightsQoL,RemainLifeExpSingleSeason,QALYDiscountRate)
end

#--------------------------------------------------------------------------
### SECTION 9: QALY LOSS CALCULATION IV (combined QALY losses)
#--------------------------------------------------------------------------

#In-hospital fatalities
ByYrOfAgeHealthEpsQALYloss[3] = ByYrOfAgeFatalQALYlossPerSeason
AllAgeHealthEpsQALYloss[3] = CombinedFatalQALYlossPerSeason

#Out of hospital mortality events
ByYrOfAgeHealthEpsQALYloss[4] = ByYrOfAge_OutOfHospFatalQALYlossPerSeason
AllAgeHealthEpsQALYloss[4] = Combined_OutOfHospFatalQALYlossPerSeason

#Assign final output to QALY losses cell.
#Total QALY loss by age (all episodes) and over all ages (across all episode types)
ByYrofAge_CombinedHealthEpsQALYloss = zeros(ByYrOfAgeGrpNum,SeasonNum,StrainNum)
HealthEpsEntityNum = length(ByYrOfAgeHealthEpsQALYloss) - 1 #Subtract one as final entry of ByYrOfAgeHealthEpsQALYloss is cumulative statistics!
ByYrofAge_CombinedHealthEpsQALYloss = CombinedHealthEpsStat(ByYrofAge_CombinedHealthEpsQALYloss,ByYrOfAgeHealthEpsQALYloss,HealthEpsEntityNum)
# for QALYLossEpsIdx = 1:4
#     #Contribution towards combined cost of all health episodes
#     ByYrofAge_CombinedHealthEpsQALYloss = ByYrofAge_CombinedHealthEpsQALYloss + ByYrOfAgeHealthEpsQALYloss[QALYLossEpsIdx]
# end

#Total QALY loss by age (all episodes) and over all ages (across all episode types)
ByYrOfAgeHealthEpsQALYloss[5] = ByYrofAge_CombinedHealthEpsQALYloss
AllAgeHealthEpsQALYloss[5] = sum(ByYrofAge_CombinedHealthEpsQALYloss,dims=1)

#--------------------------------------------------------------------------
### SECTION 10: CALCULATE VACCINE COSTS
#--------------------------------------------------------------------------

#Disaggregate VaccCostsData tuple
VaccUptake = VaccCostsData[1]::Array{Float64}
AdminChargePerVaccDose = VaccCostsData[2]::Float64
CostPerVaccineDose = VaccCostsData[3]::Float64

#Compute cost of vaccines used (apply discounting term to each season/column vector)
#Vaccine costs tuple output: (i) Administration charge; (ii) per dose cost; (iii) overall vaccine cost
VaccCosts = Array{Array{Float64}}(undef,3);
if ndims(InfectedCaseCount) == 2
    VaccDeployed = round.(VaccUptake.*vec(PopnData'))
    VaccCosts[1] = AdminChargePerVaccDose*VaccDeployed.*CostDiscountTerm
    VaccCosts[2] = CostPerVaccineDose*VaccDeployed.*CostDiscountTerm
    VaccCosts[3] = (CostPerVaccineDose + AdminChargePerVaccDose)*VaccDeployed.*CostDiscountTerm

else
    VaccDeployed = round.(VaccUptake.*PopnData')
    VaccCosts[1] = AdminChargePerVaccDose*VaccDeployed.*CostDiscountTerm'
    VaccCosts[2] = CostPerVaccineDose*VaccDeployed.*CostDiscountTerm'
    VaccCosts[3] = (CostPerVaccineDose + AdminChargePerVaccDose)*VaccDeployed.*CostDiscountTerm'

end

#--------------------------------------------------------------------------
### SECTION 11: CALCULATE VACCINE DEPLOYED (discounted) & ASSIGN TO VaccDeployedTuple
#--------------------------------------------------------------------------

#Apply discounting term to each season/column vector)
if ndims(InfectedCaseCount) == 2
    VaccDeployedDiscountedCount = VaccDeployed.*CostDiscountTerm
else
    VaccDeployedDiscountedCount = VaccDeployed.*CostDiscountTerm'
end

VaccDeployedTuple = [VaccDeployed,VaccDeployedDiscountedCount]

return ByYrOfAgeHealthEpsCounts,ByYrOfAgeHealthEpsCosts,ByYrOfAgeHealthEpsQALYloss,
    AllAgeHealthEpsCounts,AllAgeHealthEpsCosts,AllAgeHealthEpsQALYloss,VaccCosts,VaccDeployedTuple
end
