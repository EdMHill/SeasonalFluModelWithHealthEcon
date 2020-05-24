#Purpose:
#Compute the allowable vaccination price to not exceed target WTP threshold

#Differences to vanilla MinAccVaccPriceCalcJulV1.jl:
#Iterate over each alternate strategy compared to reference strategy.
#Store threshold vaccine price in a centralised array.

#Steps:
# #0: Specify relevant data from reference strategy
# #1: Specify relevant data from alternative strategy
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
using Statistics

#--------------------------------------------------------------------------
### Include required function files
#--------------------------------------------------------------------------
include("HealthEconEvalTotalCalcJulV1.jl")

#--------------------------------------------------------------------------
### Initialisation and variable declaration
#--------------------------------------------------------------------------

#Declare number of runs per strategy for health economic evaluation
HealthEconEvalRunNum = 100

#Specify number of alternative stratigies being analysed
JobNum = 10

#Initialise centralised storage array for threshold price per vaccine dose
ThresholdPricePerDoseData_ByAltStrat = zeros(HealthEconEvalRunNum,JobNum)

#--------------------------------------------------------------------------
### Step 0: Specify relevant data from reference strategy
#--------------------------------------------------------------------------
HealthEconRefStratFileName_LowRisk = "HealthEconEvalOutputData/HealthEconOutputs_LowRisk_JulV1Run#VaccAgeSweepRef.mat"
HealthEconRefStratFileName_AtRisk = "HealthEconEvalOutputData/HealthEconOutputs_AtRisk_JulV1Run#VaccAgeSweepRef.mat"

#Iterate over each alternate strategy
for JobIdx = 1:JobNum
    println("JobIdx: $JobIdx")
    #--------------------------------------------------------------------------
    ### Step 1: Specify relevant data from reference strategy and alternative strategy
    #--------------------------------------------------------------------------

    #Get current job run health economic file names
    HealthEconAltStratFileName_LowRisk = "E:/GitHub/DoH-FluVaccine/Results/HealthEconEval/HealthEconEvalOutputData/HealthEconOutputs_LowRisk_JulV1Run#VaccAgeSweep_JobID#$JobIdx.mat"
    HealthEconAltStratFileName_AtRisk = "E:/GitHub/DoH-FluVaccine/Results/HealthEconEval/HealthEconEvalOutputData/HealthEconOutputs_AtRisk_JulV1Run#VaccAgeSweep_JobID#$JobIdx.mat"

    #Append current job run of interest to file name vector
    HealthEconEvalFNames = [HealthEconRefStratFileName_LowRisk HealthEconAltStratFileName_LowRisk;
        HealthEconRefStratFileName_AtRisk HealthEconAltStratFileName_AtRisk]

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
    ### Step 8 - Data processing
    #--------------------------------------------------------------------------

    #Take relevant data from allowable vaccine cost data array.
    WTPidx = 7 #Column for specified WTP (column 5 for £20,000, column 7 for £30,000)
    ThresholdPricePerDoseData = dropdims(MaxVaccCostExclAdminCharge[:,WTPidx,:]; dims=2) #Remove singleton dimension

    #Output ThresholdPricePerDoseData to centralised storage array
    ThresholdPricePerDoseData_ByAltStrat[:,JobIdx] = ThresholdPricePerDoseData

end

#Get tenth percentile per strategy
QuantileVal = 0.1
PrctileDimn = 1 #Operate along the columns
ThresholdPricePerDose_ByAltStrat_Prctile = zeros(JobNum)

for AltStratIdx = 1:JobNum
    ThresholdPricePerDose_ByAltStrat_Prctile[AltStratIdx] = quantile(ThresholdPricePerDoseData_ByAltStrat[:,AltStratIdx],QuantileVal)
end

#Rank across strategies in descending order (based on tenth percentile)
#Those with highest willingness to pay are preferred!
StratIdx_DescendThresholdPricePerDose = sortperm(ThresholdPricePerDose_ByAltStrat_Prctile,rev=true)
DescendThresholdPricePerDose_ByAltStrat_Prctile = ThresholdPricePerDose_ByAltStrat_Prctile[StratIdx_DescendThresholdPricePerDose]
