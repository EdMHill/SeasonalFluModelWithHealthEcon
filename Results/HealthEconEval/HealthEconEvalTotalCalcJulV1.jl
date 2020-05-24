#Purpose:
#Function to output total costs and total QALY losses per replicate

#Compatible with Julia V1

#--------------------------------------------------------------------------


function HealthEconEvalTotalCalc(HealthEconEvalFNames)
#Inputs:
#   HealthEconEvalFNames - (tuple, strings) Filenames of healht economic evaluation data

#Outputs:
#   TotalMonetaryCosts, TotalQALYloss, TotalMonetaryCostsWithoutVacc, TotalDiscountVaccDeployed - (2D arrays)
#          -> row per vacciniation strat
#          -> column per health econ eval replicate


#--------------------------------------------------------------------------
### Import data from files, assign to variables
#--------------------------------------------------------------------------

#Number of vaccine strategies considered, columns of HealthEconEvalFNames
NumVaccStratsConsidered = size(HealthEconEvalFNames,2)

#Number of risk groups considered, rows of HealthEconEvalFNames
RiskGrpNum = size(HealthEconEvalFNames,1)

#Initialise cells to store all-age related data outputs
AllAgeHealthEpsCostsTuple =  Array{Array{Any}}(undef,RiskGrpNum,NumVaccStratsConsidered)
AllAgeHealthEpsQALYlossTuple = Array{Array{Any}}(undef,RiskGrpNum,NumVaccStratsConsidered)
VaccCostsTuple = Array{Array{Any}}(undef,RiskGrpNum,NumVaccStratsConsidered)
VaccDeployedTupleAllRuns = Array{Array{Any}}(undef,RiskGrpNum,NumVaccStratsConsidered)

#Iterate over risk groups and considered vaccine strategies
for FNameIdx = 1:NumVaccStratsConsidered

    RiskGrpIdx = 1
    while RiskGrpIdx <= RiskGrpNum

        vars = matread(HealthEconEvalFNames[RiskGrpIdx,FNameIdx])
        AllAgeHealthEpsCostsTuple[RiskGrpIdx,FNameIdx] = vars["AllAgeHealthEpsCosts"]
        AllAgeHealthEpsQALYlossTuple[RiskGrpIdx,FNameIdx] = vars["AllAgeHealthEpsQALYloss"]
        VaccCostsTuple[RiskGrpIdx,FNameIdx] = vars["VaccCosts"]
        VaccDeployedTupleAllRuns[RiskGrpIdx,FNameIdx] = vars["VaccDeployedTuple"]

        #Increment RiskGrpIdx
        RiskGrpIdx = RiskGrpIdx + 1
    end
end

#Acquire the number of replicates used
HealthEconEvalReplicateNum = length(VaccCostsTuple[1,1])

#--------------------------------------------------------------------------
### Compute monetary costs and QALY losses for each replicate
#--------------------------------------------------------------------------

#Get number of cell arrays in Costs tuple and QALYloss tuple
#Final entry contains sum across all health episodes under consideration
CostsTupleNum = length(AllAgeHealthEpsCostsTuple[1,1][1])
QALYlossTupleNum = length(AllAgeHealthEpsQALYlossTuple[1,1][1])

#Get number of tuple arrays in VaccCosts tuple
#Final entry contains sum across administration charge and dose cose
VaccCostsTupleNum = length(VaccCostsTuple[1,1][1])

#Initialise storage vectors
TotalMonetaryCosts = zeros(NumVaccStratsConsidered,HealthEconEvalReplicateNum)
TotalQALYloss = zeros(NumVaccStratsConsidered,HealthEconEvalReplicateNum)
TotalMonetaryCostsWithoutVacc = zeros(NumVaccStratsConsidered,HealthEconEvalReplicateNum)
TotalDiscountVaccDeployed = zeros(NumVaccStratsConsidered,HealthEconEvalReplicateNum)

#Per health economic evaluation replicate
for ii = 1:HealthEconEvalReplicateNum

    jj = 1;
    while jj<=NumVaccStratsConsidered

        #Calculate total costs (across time horizon) in each risk group
        #Calculate total QALY losses (across time horizon) in each risk
        #group
        MonetaryCostsCumulativeSum = 0
        QALYlossCumulativeSum = 0
        MonetaryCostsWithoutVaccCumulativeSum = 0
        VaccDeployedCumulativeSum = 0
        kk = 1
        while kk <= RiskGrpNum
            #Sum vacc costs
            #Sum Health episode costs
            MonetaryCostsCumulativeSum = MonetaryCostsCumulativeSum +
                                        sum(VaccCostsTuple[kk,jj][ii][VaccCostsTupleNum]) +
                                            sum(AllAgeHealthEpsCostsTuple[kk,jj][ii][CostsTupleNum])

            #Add QALY losses for risk group kk
            QALYlossCumulativeSum = QALYlossCumulativeSum +
                                        sum(AllAgeHealthEpsQALYlossTuple[kk,jj][ii][QALYlossTupleNum])


            #Add non-vaccine costs for risk group kk
            MonetaryCostsWithoutVaccCumulativeSum = MonetaryCostsWithoutVaccCumulativeSum +
                                                        sum(AllAgeHealthEpsCostsTuple[kk,jj][ii][CostsTupleNum])

            #Add vaccines deployed (discounted) for risk group kk
            #Index 2 of relevant entry of VaccDeployedTuple
            VaccDeployedCumulativeSum = VaccDeployedCumulativeSum +
                                            sum(VaccDeployedTupleAllRuns[kk,jj][ii][2])

            #Increment index
            kk = kk + 1
        end

        #Calculate total costs (across time horizon)
        #Sum vacc costs
        #Sum Health episode costs
        TotalMonetaryCosts[jj,ii] = MonetaryCostsCumulativeSum

        #Calculate total QALY losses (across time horizon)
        TotalQALYloss[jj,ii] = QALYlossCumulativeSum

        #Calculate total non-vaccine related monetary costs (across time
        #horizon)
        TotalMonetaryCostsWithoutVacc[jj,ii] = MonetaryCostsWithoutVaccCumulativeSum

        #Get number of vaccines used
        TotalDiscountVaccDeployed[jj,ii] = VaccDeployedCumulativeSum

        #Increment index
        jj = jj + 1
    end
end

return TotalMonetaryCosts,TotalQALYloss,TotalMonetaryCostsWithoutVacc,TotalDiscountVaccDeployed
end
