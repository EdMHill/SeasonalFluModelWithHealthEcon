#Purpose:
#File to store functions used to compute, for a fixed vaccine cost, the optimal age coverage of a
#low-risk vaccination programme (on top of an at-risk programme)

# Also has a function (OptCovFixedVaccPriceFn) finding optimum stratgies whilst
# mandating coverage of all those aged 65 and overs

#Outline:
# Stage 0: Import data & assign to variables
# Stage 1: Get difference in vaccines bought between alternate & reference strategy
# Stage 2: Iterate over fixed vaccine values, get cost relative to threshold value
# Stage 3: Per strategy, Multiply difference in threshold & fixed vaccine price by number of vaccines required
# Stage 4: Get quantile of costs relative to threshold value
# Stage 5: Find strategy retaining greatest value (compared to reference strategy)
# Stage 6: Map optimal strategy to age coverage values

#--------------------------------------------------------------------------


function  OptCovFixedVaccPriceFn(AltVaccStratCoveragesFile::String,
                                VaccinesBoughtFile::String,
                                ThresholdPricePerDoseFile::String,
                                QALYSliceIdx::Int64,
                                FixedVaccPriceVals::StepRange{Int64,Int64},
                                QuantileVal::Float64,
                                OutOfRangeYouth::Int64,
                                OutOfRangeElder::Int64)
#Inputs:
# AltVaccStratCoveragesFile, VaccinesBought_File, ThresholdPricePerDoseFile - Data input files
# QALYSliceIdx - Slice of threshold price per dose array to access, corresponding to a particular worth per QALY
#                  Slice 1: 0k, Slice 2: 5k,... , Slice 5: 20k, Slice 7: 30k.
# FixedVaccPriceVals - Vector of fixed price vaccine costs to use.
# QuantileVal - Desired percentile to be computed
# OutOfRangeYouth, OutOfRangeElder - For plotting purposes, if young-age centric and/or elder-age centric programme not deemed optimal return these values instead.

#Outputs:
# OptimVaccStrat - Array of age coverage. Row per fixed vacc price. Cols: [YouthAge ElderAge]

#--------------------------------------------------------------------------
### STAGE 0: Import data & assign to variables
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
### Load paired strategy data
#--------------------------------------------------------------------------
AltVaccStratCoverages = readdlm(AltVaccStratCoveragesFile,',')

#--------------------------------------------------------------------------
### Import vaccines deployed data
#--------------------------------------------------------------------------
VaccinesBought_File = matopen(VaccinesBoughtFile)

#Assign discounted vaccine totals to variables
VaccinesBought_RefStrat = read(VaccinesBought_File,"VaccDeployedPerStrat_RefStrat")
VaccinesBought_AltStrat = read(VaccinesBought_File,"VaccDeployedPerStrat")

#--------------------------------------------------------------------------
### Import threshold vaccine dose price data
#--------------------------------------------------------------------------
ThresholdPricePerDose_File = matopen(ThresholdPricePerDoseFile)

#--------------------------------------------------------------------------
### Assign selected threshold vaccine dose price data to variables
#--------------------------------------------------------------------------
ThresholdPricePerDoseData_ByAltStrat = read(ThresholdPricePerDose_File,"ThresholdPricePerDoseData_ByAltStrat")

ThresholdPricePerDose = ThresholdPricePerDoseData_ByAltStrat[:,QALYSliceIdx,:]

#Get number of alternate strategies being analysed
ReplicateNum = size(ThresholdPricePerDose,1)
AltStratNum = size(ThresholdPricePerDose,2)

#--------------------------------------------------------------------------
### STAGE 1: Get difference in vaccines bought between alternate & ref strategy
#--------------------------------------------------------------------------
VaccinesBought_AltVsRefStrat = VaccinesBought_AltStrat .- VaccinesBought_RefStrat

#--------------------------------------------------------------------------
### STAGE 2: Iterate over fixed vaccine values, get cost relative to threshold value
#--------------------------------------------------------------------------
FixedVaccPriceNum = length(FixedVaccPriceVals)

#Initialise array to store difference between fixed price and threshold price
ThresholdVsFixedPrice = zeros(ReplicateNum,AltStratNum,FixedVaccPriceNum)

for FixedVaccValItr= 1:FixedVaccPriceNum
    ThresholdVsFixedPrice[:,:,FixedVaccValItr] = ThresholdPricePerDose .- FixedVaccPriceVals[FixedVaccValItr]
end


#--------------------------------------------------------------------------
### STAGE 3: Per strategy, multiply difference in threshold & fixed vaccine price by number of vaccines required
#--------------------------------------------------------------------------

#Initialise storage array
OverallPriceDiff = zeros(ReplicateNum,AltStratNum,FixedVaccPriceNum)

#Iterate over each fixed vaccine price
for FixedVaccValItr= 1:FixedVaccPriceNum
    #Multiply each column (strategy) by number of additional vaccines needed for carrying out that strategy
    #To do this, use transpose of VaccinesBought_StandardDR_AltVsRefStrat, returning a row vector of 1xAltStratNum
    OverallPriceDiff[:,:,FixedVaccValItr] = VaccinesBought_AltVsRefStrat'.*ThresholdVsFixedPrice[:,:,FixedVaccValItr]
end

#--------------------------------------------------------------------------
### STAGE 4: Get quantile of costs relative to threshold value
#--------------------------------------------------------------------------

#Get desired quantile value across 100 replicates (across first dimension)
ThresholdVsFixedPriceQuantile = zeros(AltStratNum,FixedVaccPriceNum)
for FixedVaccPriceItr = 1:FixedVaccPriceNum
    for AltStratItr = 1:AltStratNum
        ThresholdVsFixedPriceQuantile[AltStratItr,FixedVaccPriceItr] = quantile(OverallPriceDiff[:,AltStratItr,FixedVaccPriceItr],QuantileVal)
    end
end

#ThresholdVsFixedPriceQuantile = dropdims(median(OverallPriceDiff,dims=1),dims=1)

#--------------------------------------------------------------------------
### STAGE 5: Find strategy retaining greatest value (compared to reference strategy)
#--------------------------------------------------------------------------

#Have 2D arrays: Column by fixed vaccine price, row by alternate strategy
#Take maximum across each column & return indexes
MaxVals, MaxIdxs = findmax(ThresholdVsFixedPriceQuantile,dims=1)

#Check for any maximums that are below zero
#Reassign values as NaN
MaxVals[MaxVals.<0] .= NaN

#--------------------------------------------------------------------------
### STAGE 6: Map optimal strategy to age coverage values
#--------------------------------------------------------------------------

#Initialise storage arrays
OptimVaccStrat = zeros(FixedVaccPriceNum,2)

#Iterate over fixed vaccine cost
for FixedVaccValItr = 1:FixedVaccPriceNum

    #Get fixed vaccine price for current iteration
    SelectedVaccPrice = FixedVaccPriceVals[FixedVaccValItr]

    #Check if max value is NaN, if so,
    #map young-age centric to 0yrs; elder-age centric to 102
    if isnan(MaxVals[FixedVaccValItr])
        OptimVaccStrat[FixedVaccValItr,:] = [OutOfRangeYouth OutOfRangeElder]
    else #Positive max value. Assign optimal strategy to storage array.
        #Access cartesian coordinate, row index.
        OptimStrat = MaxIdxs[FixedVaccValItr][1]

        #Map strategy number to young and elder age strategy.
        OptimVaccStrat[FixedVaccValItr,:] = AltVaccStratCoverages[OptimStrat,:]

        #Check for young-age centric only programmes. Map elder age coverage to 102.
        YoungAgeOnlyProgIdx = isnan.(OptimVaccStrat[:,2])
        OptimVaccStrat[YoungAgeOnlyProgIdx,2] .= OutOfRangeElder

        #Check for elder-age centric only programmes. Map young age coverage to two less than lowest age that's part of programme
        ElderAgeOnlyProgIdx = isnan.(OptimVaccStrat[:,1])
        OptimVaccStrat[ElderAgeOnlyProgIdx,1] .= OutOfRangeYouth
    end
end
    return OptimVaccStrat::Array{Float64,2}
end


# Finding optimum stratgies whilst
# mandating coverage of all those aged 65 and overs
function  OptCovFixedVaccPriceFn_IncOver65s(AltVaccStratCoveragesFile::String,
                                                VaccinesBoughtFile::String,
                                                ThresholdPricePerDoseFile::String,
                                                QALYSliceIdx::Int64,
                                                FixedVaccPriceVals::StepRange{Int64,Int64},
                                                QuantileVal::Float64,
                                                OutOfRangeYouth::Int64,
                                                OutOfRangeElder::Int64)
#Inputs:
# AltVaccStratCoveragesFile, VaccinesBought_File, ThresholdPricePerDoseFile - Data input files
# QALYSliceIdx - Slice of threshold price per dose array to access, corresponding to a particular worth per QALY
#                  Slice 1: 0k, Slice 2: 5k,... , Slice 5: 20k, Slice 7: 30k.
# FixedVaccPriceVals - Vector of fixed price vaccine costs to use.
# QuantileVal - Desired percentile to be computed
# OutOfRangeYouth, OutOfRangeElder - For plotting purposes, if young-age centric and/or elder-age centric programme not deemed optimal return these values instead.

#Outputs:
# OptimVaccStrat - Array of age coverage. Row per fixed vacc price. Cols: [YouthAge ElderAge]

#--------------------------------------------------------------------------
### STAGE 0: Import data & assign to variables
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
### Set up array indexes, corresponding to rows with elder age vaccine coverage
### including those 65 and above
#--------------------------------------------------------------------------
array_idx_65andabove = [1:31;51:81;100:129;148:176;195:222;241:267;286:311;330:354;
                            373:396;415:437;456:477;496:516;535:554;573:591;610:627;
                            646:662;681:696;715:729;748:761;780:792;811:822;841:851;
                            870:879;898:906;925:932;951:957;976:981;1000:1004;
                            1023:1026;1045:1047;1066:1067;1086;1326:1357]
n_65andabove_strat = length(array_idx_65andabove)

#--------------------------------------------------------------------------
### Load paired strategy data
#--------------------------------------------------------------------------
AltVaccStratCoverages = readdlm(AltVaccStratCoveragesFile,',')

#--------------------------------------------------------------------------
### Import vaccines deployed data
#--------------------------------------------------------------------------
VaccinesBought_File = matopen(VaccinesBoughtFile)

#Assign discounted vaccine totals to variables
VaccinesBought_RefStrat = read(VaccinesBought_File,"VaccDeployedPerStrat_RefStrat")
VaccinesBought_AltStrat = read(VaccinesBought_File,"VaccDeployedPerStrat")

#--------------------------------------------------------------------------
### Import threshold vaccine dose price data
#--------------------------------------------------------------------------
ThresholdPricePerDose_File = matopen(ThresholdPricePerDoseFile)

#--------------------------------------------------------------------------
### Assign selected threshold vaccine dose price data to variables
#--------------------------------------------------------------------------
ThresholdPricePerDoseData_ByAltStrat = read(ThresholdPricePerDose_File,"ThresholdPricePerDoseData_ByAltStrat")

ThresholdPricePerDose = ThresholdPricePerDoseData_ByAltStrat[:,QALYSliceIdx,:]

#Get number of alternate strategies being analysed
ReplicateNum = size(ThresholdPricePerDose,1)
AltStratNum = size(ThresholdPricePerDose,2)

#--------------------------------------------------------------------------
### STAGE 1: Get difference in vaccines bought between alternate & ref strategy
#--------------------------------------------------------------------------
VaccinesBought_AltVsRefStrat = VaccinesBought_AltStrat .- VaccinesBought_RefStrat

#--------------------------------------------------------------------------
### STAGE 2: Iterate over fixed vaccine values, get cost relative to threshold value
#--------------------------------------------------------------------------
FixedVaccPriceNum = length(FixedVaccPriceVals)

#Initialise array to store difference between fixed price and threshold price
ThresholdVsFixedPrice = zeros(ReplicateNum,AltStratNum,FixedVaccPriceNum)

for FixedVaccValItr= 1:FixedVaccPriceNum
    ThresholdVsFixedPrice[:,:,FixedVaccValItr] = ThresholdPricePerDose .- FixedVaccPriceVals[FixedVaccValItr]
end


#--------------------------------------------------------------------------
### STAGE 3: Per strategy, multiply difference in threshold & fixed vaccine price by number of vaccines required
#--------------------------------------------------------------------------

#Initialise storage array
OverallPriceDiff = zeros(ReplicateNum,AltStratNum,FixedVaccPriceNum)

#Iterate over each fixed vaccine price
for FixedVaccValItr= 1:FixedVaccPriceNum
    #Multiply each column (strategy) by number of additional vaccines needed for carrying out that strategy
    #To do this, use transpose of VaccinesBought_StandardDR_AltVsRefStrat, returning a row vector of 1xAltStratNum
    OverallPriceDiff[:,:,FixedVaccValItr] = VaccinesBought_AltVsRefStrat'.*ThresholdVsFixedPrice[:,:,FixedVaccValItr]
end

#--------------------------------------------------------------------------
### STAGE 4: Get quantile of costs relative to threshold value
#--------------------------------------------------------------------------

#Get desired quantile value across 100 replicates (across first dimension)
ThresholdVsFixedPriceQuantile = zeros(n_65andabove_strat,FixedVaccPriceNum)
for FixedVaccPriceItr = 1:FixedVaccPriceNum
    for AltStratItr = 1:n_65andabove_strat
        alt_strat_idx = array_idx_65andabove[AltStratItr]
        ThresholdVsFixedPriceQuantile[AltStratItr,FixedVaccPriceItr] = quantile(OverallPriceDiff[:,alt_strat_idx,FixedVaccPriceItr],QuantileVal)
    end
end

#ThresholdVsFixedPriceQuantile = dropdims(median(OverallPriceDiff,dims=1),dims=1)

#--------------------------------------------------------------------------
### STAGE 5: Find strategy retaining greatest value (compared to reference strategy)
#--------------------------------------------------------------------------

#Have 2D arrays: Column by fixed vaccine price, row by alternate strategy
#Take maximum across each column & return indexes
MaxVals, MaxIdxs = findmax(ThresholdVsFixedPriceQuantile,dims=1)

#Check for any maximums that are below zero
#Reassign values as NaN
MaxVals[MaxVals.<0] .= NaN

#--------------------------------------------------------------------------
### STAGE 6: Map optimal strategy to age coverage values
#--------------------------------------------------------------------------

#Initialise storage arrays
OptimVaccStrat = zeros(FixedVaccPriceNum,2)

#Iterate over fixed vaccine cost
for FixedVaccValItr = 1:FixedVaccPriceNum

    #Get fixed vaccine price for current iteration
    SelectedVaccPrice = FixedVaccPriceVals[FixedVaccValItr]

    #Check if max value is NaN, if so,
    #map young-age centric to 0yrs; elder-age centric to 102
    if isnan(MaxVals[FixedVaccValItr])
        OptimVaccStrat[FixedVaccValItr,:] = [OutOfRangeYouth OutOfRangeElder]
    else #Positive max value. Assign optimal strategy to storage array.
        #Access cartesian coordinate, row index.
        OptimStrat_idx = MaxIdxs[FixedVaccValItr][1]

        #Map strategy number to young and elder age strategy.
        OptimStrat = array_idx_65andabove[OptimStrat_idx]
        OptimVaccStrat[FixedVaccValItr,:] = AltVaccStratCoverages[OptimStrat,:]

        #Check for young-age centric only programmes. Map elder age coverage to 102.
        YoungAgeOnlyProgIdx = isnan.(OptimVaccStrat[:,2])
        OptimVaccStrat[YoungAgeOnlyProgIdx,2] .= OutOfRangeElder

        #Check for elder-age centric only programmes. Map young age coverage to two less than lowest age that's part of programme
        ElderAgeOnlyProgIdx = isnan.(OptimVaccStrat[:,1])
        OptimVaccStrat[ElderAgeOnlyProgIdx,1] .= OutOfRangeYouth
    end
end
    return OptimVaccStrat::Array{Float64,2}
end
