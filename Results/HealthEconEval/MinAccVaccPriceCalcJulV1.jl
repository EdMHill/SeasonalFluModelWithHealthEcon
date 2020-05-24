#Purpose:
#Compute the allowable vaccination price to not exceed target WTP threshold

#Steps:
# #1: Specify relevant data from reference strategy and alternative strategy
# #2: For each considered strategy, get releveant costs, total QALY losses
# #3: For each alternative strategy, relative to refrence policy, get (non-vacc associated) incremental costs and QALY gains
# #4: Get maximum allowable incremental cost, based on QALYS gained
# #5: Per alternative strategy, get permitable spend
# #6: Per alternative strategy, get allowable spend on vaccines & discounted vaccine count over time horizon
# #7: Compute absolute & relative vaccine cost to attain WTP thershold

#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#Load required packages
#--------------------------------------------------------------------------
using MAT

#--------------------------------------------------------------------------
### Include required function files
#--------------------------------------------------------------------------
include("HealthEconEvalTotalCalcJulV1.jl")

#--------------------------------------------------------------------------
### Step 1: Specify relevant data from reference strategy and alternative strategy
#--------------------------------------------------------------------------

#Collate input filenames
# HealthEconEvalFNames = ["HealthEconEvalOutputData/HealthEconOutputs_LowRisk_JulV1Run#OptimFit13G_NoVacc.mat" "HealthEconEvalOutputData/HealthEconOutputs_LowRisk_JulV1Run#OptimFit13G_HistoricalVacc.mat";
#                             "HealthEconEvalOutputData/HealthEconOutputs_AtRisk_JulV1Run#OptimFit13G_NoVacc.mat" "HealthEconEvalOutputData/HealthEconOutputs_AtRisk_JulV1Run#OptimFit13G_HistoricalVacc.mat"]

# HealthEconEvalFNames = ["HealthEconEvalOutputData/HealthEconOutputs_LowRisk_JulV1Run#OptimFit13G_HistAtRiskOnlyVacc.mat" "HealthEconEvalOutputData/HealthEconOutputs_LowRisk_JulV1Run#OptimFit13G_HistAtRiskPlusChildVacc.mat" "HealthEconEvalOutputData/HealthEconOutputs_LowRisk_JulV1Run#OptimFit13G_HistAtRiskPlusElderlyVacc.mat" "HealthEconEvalOutputData/HealthEconOutputs_LowRisk_JulV1Run#OptimFit13G_HistoricalVacc.mat";
#                         "HealthEconEvalOutputData/HealthEconOutputs_AtRisk_JulV1Run#OptimFit13G_HistAtRiskOnlyVacc.mat" "HealthEconEvalOutputData/HealthEconOutputs_AtRisk_JulV1Run#OptimFit13G_HistAtRiskPlusChildVacc.mat" "HealthEconEvalOutputData/HealthEconOutputs_AtRisk_JulV1Run#OptimFit13G_HistAtRiskPlusElderlyVacc.mat" "HealthEconEvalOutputData/HealthEconOutputs_AtRisk_JulV1Run#OptimFit13G_HistoricalVacc.mat"]

HealthEconEvalFNames = ["HealthEconEvalOutputData/HealthEconOutputs_LowRisk_JulV1Run#TransContactArrayTest_AtRiskStrat.mat" "HealthEconEvalOutputData/HealthEconOutputs_LowRisk_JulV1Run#VaccAgeSweepTest_JobID#500.mat";
                            "HealthEconEvalOutputData/HealthEconOutputs_AtRisk_JulV1Run#TransContactArrayTest_AtRiskStrat.mat" "HealthEconEvalOutputData/HealthEconOutputs_AtRisk_JulV1Run#VaccAgeSweepTest_JobID#500.mat"]

#Number of vaccine strategies considered, columns of HealthEconEvalFNames
NumVaccStratsConsidered = size(HealthEconEvalFNames,2)

#Flag variable. Declare if reference strategy should include vaccine costs
#or not
RefStratIncludeVaccCostsFlag = 0 #0 for not including vaccine costs. 1 to include vaccine costs.

#Set baseline reference price for administration charge per vaccine dose
AdminChargePerVaccDose = [10.;10.]
AdminChargeOffsetToRefStrat = AdminChargePerVaccDose[2:end] .- AdminChargePerVaccDose[1] #Differences between alternate strat admin charge to reference strat

if length(AdminChargePerVaccDose) != NumVaccStratsConsidered
    error("Amount of AdminChargePerVaccDose quantities does not match number of strategies under consideration.")
end

#Specify vaccine wastage amount
VaccWastePercentage = [10.;10.]

if (sum(VaccWastePercentage.<0)>0) || (sum(VaccWastePercentage.>=100)>0)
   error("VaccWastePercentage current value is $VaccWastePercentage. VaccWastePercentage is a percetange, must take value between 0 and 100")
end

if length(VaccWastePercentage) != NumVaccStratsConsidered
    error("Amount of vaccine wastage quantities does not match number of strategies under consideration.")
end

VaccDoseScaleForWastage = 100 ./(100 .- VaccWastePercentage)

#--------------------------------------------------------------------------
### Step 2: For each considered strategy, get relevant costs, total QALY losses
#--------------------------------------------------------------------------
TotalMonetaryCosts,TotalQALYloss,TotalMonetaryCostsWithoutVacc,UnwastedDiscountVaccDeployed =
                 HealthEconEvalTotalCalc(HealthEconEvalFNames)

#Acquire the number of replicates used
HealthEconEvalReplicateNum = size(TotalMonetaryCosts,2)

#Scale vacc deployed to account for vaccine wastage
TotalDiscountVaccDeployed = VaccDoseScaleForWastage.*UnwastedDiscountVaccDeployed

#--------------------------------------------------------------------------
### Step 3: For each alternative strategy, relative to refrence policy,
###         get (non-vacc associated) incremental costs and QALY gains
#--------------------------------------------------------------------------
ReferenceVaccStratID = 1 #Define which row of CostsTupleNum, QALYlossTupleNum corresponds to reference strategy
AlternativeVaccStratID = [2] #Vector of IDs to be considered

#Iterate through alternative programmes
#Calculate required statistics
QALYgained = zeros(NumVaccStratsConsidered-1,HealthEconEvalReplicateNum)
NonVaccIncrementalCosts = zeros(NumVaccStratsConsidered-1,HealthEconEvalReplicateNum)
for ii = 1:NumVaccStratsConsidered-1
    AlternativeVaccStratArrayIdx = AlternativeVaccStratID[ii]

    #Compute drop/increase in QALYs relative to reference policy.
    QALYgained[ii,:] = TotalQALYloss[1,:] - TotalQALYloss[AlternativeVaccStratArrayIdx,:]

    #Compute non-vaccine associated incremental costs (relative to
    #reference policy)
    #If negative, alternative strategy has lower (non-vacc) costs than baseline strategy
    NonVaccIncrementalCosts[ii,:] = TotalMonetaryCostsWithoutVacc[AlternativeVaccStratArrayIdx,:] -
                                TotalMonetaryCostsWithoutVacc[1,:]

end

#--------------------------------------------------------------------------
### Step 4: Get maximum allowable incremental cost, based on QALYS gained
#--------------------------------------------------------------------------
#Specify Willingness-to-pay (WTP) threshold
WTP = 0:5000:40000
WTPValsTested = length(WTP)

#For each WTP value, get max incremental cost
MaxIncrementalCost = zeros(NumVaccStratsConsidered-1,HealthEconEvalReplicateNum,WTPValsTested)
for ii = 1:WTPValsTested
    MaxIncrementalCost[:,:,ii] = QALYgained.*WTP[ii]
    #3D array. Row by alternative vacc strategy. Column by replicate.
    #           Slice by WTP threshold value
end

#--------------------------------------------------------------------------
### Step 5: Per alternative strategy, get permitable spend
#--------------------------------------------------------------------------

#Get total monetary costs for reference strategy (first row of TotalMonetaryCosts)
if RefStratIncludeVaccCostsFlag == 1
    RefStratTotalSpend = TotalMonetaryCosts[1,:]
else
    RefStratTotalSpend = TotalMonetaryCostsWithoutVacc[1,:] #Vaccine costs not included in reference strategy
end

#Repeat RefStratTotalSpend to match dimensions of MaxIncrementalCost array
AltVaccStratConsidered = NumVaccStratsConsidered-1
RefStratTotalSpendRepeated = repeat(RefStratTotalSpend',outer=(AltVaccStratConsidered,1,WTPValsTested))
#Note, take transpose of RefStratTotalSpend (1D column vector) to get a row vector, which is then repeated.

#Alternate policy permitted spend (to reference policy costs, add allowed incremental cost)
AltVaccStratPermittedSpend = RefStratTotalSpendRepeated + MaxIncrementalCost

#--------------------------------------------------------------------------
### Step 6: Per alternative strategy, get allowable spend on vaccines
###         & discounted vaccine count over time horizon
#--------------------------------------------------------------------------

#Get allowable vaccine related incremental cost, per WTP threshold
#If NonVaccIncrementalCosts negative, permits additional spend on vaccine
#before threshold is reached
AltVaccSchemeTotalMonetaryCostsWithoutVacc = TotalMonetaryCostsWithoutVacc[2:end,:]
AltVaccSchemeTotalMonetaryCostsWithoutVaccRepeated = repeat(AltVaccSchemeTotalMonetaryCostsWithoutVacc,outer=(1,1,WTPValsTested))
AltVaccSchemePermittedVaccCost = AltVaccStratPermittedSpend - AltVaccSchemeTotalMonetaryCostsWithoutVaccRepeated

#Discounted number of vaccines deployed per strategy
#Get vaccine deployed total per alternative vaccine policy
AltSchemeVaccsDeployed = TotalDiscountVaccDeployed[2:end,:]
RefSchemeVaccsDeployed = TotalDiscountVaccDeployed[1,:]

#Get amount of ADDITIONAL vaccines deployed per alternative vaccine policy (relative to reference strategy)
AltVaccStratNum = NumVaccStratsConsidered - 1
AltSchemeExtraVaccsDeployed = zeros(AltVaccStratNum,HealthEconEvalReplicateNum)
for AltVaccStratIdx = 1:length(AlternativeVaccStratID)
    AltSchemeExtraVaccsDeployed[AltVaccStratIdx,:] = AltSchemeVaccsDeployed[AltVaccStratIdx,:] .- TotalDiscountVaccDeployed[1,:]

    #If using a more costly admin fee, the permitted spend on vaccine units
    #will go down!
    AltVaccSchemePermittedVaccCost[AltVaccStratIdx,:,:] = AltVaccSchemePermittedVaccCost[AltVaccStratIdx,:,:] .-
                                                       (AdminChargeOffsetToRefStrat[AltVaccStratIdx]*RefSchemeVaccsDeployed)
end
#--------------------------------------------------------------------------
### Step 7A: Compute vaccine cost to attain WTP thershold
#--------------------------------------------------------------------------
AltVaccStratConsidered = NumVaccStratsConsidered-1
MaxVaccCost = zeros(HealthEconEvalReplicateNum,WTPValsTested,AltVaccStratConsidered)

#Also output cost per vaccine dose, not including administration charge
#Weight administration charge by proportion of vaccines unwasted
AltStratVaccDeployedPropn = (1 .- (VaccWastePercentage[2:end]./100))
MaxVaccCostExclAdminCharge = zeros(HealthEconEvalReplicateNum,WTPValsTested,AltVaccStratConsidered)

#Loop over each alternative strategy. Get WTP values
for VaccStratIdx = 1:AltVaccStratConsidered

    #Check if with VaccStratIdx strategy the amount of vaccines used
    #differs from baseline strategy.
    if AltSchemeExtraVaccsDeployed[VaccStratIdx,1] == 0
        #If the same, want to divide by amount of vaccs deployed, get additional spend per vaccine
        StratIdx_MaxVaccCost = AltVaccSchemePermittedVaccCost[VaccStratIdx,:,:]/AltSchemeVaccsDeployed[VaccStratIdx,1]
        MaxVaccCost[:,:,VaccStratIdx] = StratIdx_MaxVaccCost
        MaxVaccCostExclAdminCharge[:,:,VaccStratIdx] = StratIdx_MaxVaccCost #NO ADMIN FEE TO CONSIDER! No extra vaccines used. So informing extra cost relative to ref strategy overall uni price (dose + admin charge)

    else
        #Otherwise, divide by number of EXTRA vaccines, get threshold price per
        #extra dose
        StratIdx_MaxVaccCost = AltVaccSchemePermittedVaccCost[VaccStratIdx,:,:]/AltSchemeExtraVaccsDeployed[VaccStratIdx,1]
        MaxVaccCost[:,:,VaccStratIdx] = StratIdx_MaxVaccCost
        MaxVaccCostExclAdminCharge[:,:,VaccStratIdx] = StratIdx_MaxVaccCost .- (AltStratVaccDeployedPropn[VaccStratIdx].*AdminChargePerVaccDose[VaccStratIdx+1])
    end
end

#--------------------------------------------------------------------------
### Step 7B: Compute relative vaccine cost (to reference strategy) to attain WTP thershold
#--------------------------------------------------------------------------

#Need to access one risk group file per vaccine strategy considered
#Vaccine cost is consistent across all risk groups, therefore only one
#grp needs to be used
vars = matread(HealthEconEvalFNames[1,1]) #For simplicity, use first row

#Sum both administration charge (column 1 of VaccineDoseCosts) and dose price (column 2 of VaccineDoseCosts)
RefVaccStratVaccDoseCosts = vars["VaccineDoseCosts"]
RefVaccStratAdminChargePerDose = RefVaccStratVaccDoseCosts[:,1]
RefVaccStratDoseCost = RefVaccStratVaccDoseCosts[:,2]
RefVaccStratCostPerVacc = RefVaccStratAdminChargePerDose + RefVaccStratDoseCost

#Replicate. Match array dimensions of MaxVaccCost
RefVaccStratCostPerVaccReplicated = repeat(RefVaccStratCostPerVacc,outer=(1,WTPValsTested,AltVaccStratConsidered))

#Get relative cost
#Subtract vaccine cost used in reference strategy runs from max. absolute unit cost
RelAllowableVaccCost = MaxVaccCost - RefVaccStratCostPerVaccReplicated

#Maximum relative cost accounting for altered administration charge
RelAdminChargeCost = 0
RelAllowableVaccCostExlAdminCharge = RelAllowableVaccCost .- RelAdminChargeCost
