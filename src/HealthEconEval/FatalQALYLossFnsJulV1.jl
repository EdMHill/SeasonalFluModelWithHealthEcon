#Purpose:
#House functions used to compute QALY losses due to fatalities
#--------------------------------------------------------------------------


function  DiscountedRemainLifeExp(age,RemainLifeExp,WeightsQoL,DiscountRate)
#Iteratively calculate the quality adjusted value of remaining life
#expectancy (and discount ing also taken into account)

#Inputs:
#   age - (integer)
#   WeightsQoL - (Vector) Entry per year of age
#   RemainLifeExp - (Vector) Entry per year of age (with relevant year extracted before input to this function
#   DiscountRate - (float) To be applied to all monteary and health valued outcomes
#Outputs:
#   DiscountedLifeExp - For input age, returna djusted remaining life expectancy

#Access relevant column of RemainLifeExp
AgeArrayIdx = convert(Int64,(age + 1))
RemainLifeExpForCurrentAge = round(RemainLifeExp[AgeArrayIdx]) #Round to full years

#Iterate over remaining number of expected life years. Discount quality
#adjusted value of remaining life expectancy
DiscountedLifeExp = 0

if RemainLifeExpForCurrentAge > 0
    for ii = 1:RemainLifeExpForCurrentAge

        #Deal with case where adjusted life expectancy exceeds age 100!
        IncrementedAge = convert(Int64,age + ii)
        if IncrementedAge > 100
            QoLContribution = WeightsQoL[end]/((1 + DiscountRate)^ii)
        else #Otherwise, can access relevant row of WeightsQoL as normal
            QoLContribution = WeightsQoL[IncrementedAge]/((1 + DiscountRate)^ii);
        end
        DiscountedLifeExp = DiscountedLifeExp + QoLContribution
    end
end

return DiscountedLifeExp

end

function FatalQALYLoss(SimnTime,DeathCounts,WeightsQoL,RemainLifeExp,DiscountRate)
#Loss of QALYS from fatal cases. Calculated within this function for each age

#Inputs:
#   SimnTime - (Scalar, interger) Years/Seasons beyond initial season (first season simulated has value 0!)
#   DeathCounts - (Array) Entry per year of age. Slice per strain
#   WeightsQoL - (Vector) Entry per year of age
#   RemainLifeExp - (Vector) Entry per year of age (with relevant year extracted before input to this function
#   DiscountRate - (float) To be applied to all monteary and health valued outcomes

#Outputs:
#   CombinedFatalQALYloss - Loss of QALYS from all fatal cases

#Assign number of age classes to variable
SingleYrOfAgeNum = size(DeathCounts,1);

#Check number of age classes in use
#Greater than 1, using age-sructured model
#Equal to 1, using non-age structured model
if SingleYrOfAgeNum > 1
    #Iterate through each age, compute discounted quality adjusted
    #value of remianing life expectancy for age a
    DiscountedLifeExp = zeros(SingleYrOfAgeNum,1)
    for age = 0:100
        AgeArrayIdx = age + 1
        DiscountedLifeExp[AgeArrayIdx] = DiscountedRemainLifeExp(age,RemainLifeExp,WeightsQoL,DiscountRate)
    end
else #Non-age structured model
    #For all deaths, take age of death to be 40
    age = 40
    DiscountedLifeExp = DiscountedRemainLifeExp(age,RemainLifeExp,WeightsQoL,DiscountRate)
end

#Fatal case QALY loss per age group
ByYrOfAge_FatalQALYloss = DeathCounts.*DiscountedLifeExp.*((1/(1+DiscountRate))^SimnTime)

#Sum over all age groups
CombinedFatalQALYloss = sum(ByYrOfAge_FatalQALYloss,dims=1)

return CombinedFatalQALYloss, ByYrOfAge_FatalQALYloss

end
