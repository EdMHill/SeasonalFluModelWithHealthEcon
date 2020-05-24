#Purpose:
#Functions to populate health economic parameter tables (by year of age)
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
### HealthEconParamByAgeConstruction
#--------------------------------------------------------------------------
function HealthEconParamByAgeConstruction(ParamVals,UpperAgeBound,EntityNum,SingleYrAgeGrpNum,StrainNum)

#Inputs:
#   ParamVals - (3D array)  Row by age band. Column per episode type; Slice by influenza type (Influenza A & Influenza B)
#           -> Cols: non-fatal inpatient; outpatient; fatal case in-hospital; fatal case out of hospital (ratio to in hospital).
#   UpperAgeBound - (vector, integers) Ages pinpointing a change in QALY loss value
#   EntityNum - (scalar, integer) Number of events being parameterised
#   SingleYrAgeGrpNum - (scalar, integer) Number of single year age classes in use (typically 101, 0 to 100+)
#   StrainNum - (scalar, interger)

#Outputs:
#    HealthEpsParamByAge - (3D array)
#          -> Row per single year age group,
#          -> Column per episode type (GP consultation, non-fatal inpatient, outpatient, fatal case in-hospital, fatal case out of hospital)
#          -> slice per strain

#Get number of groupings in use
AgeBandNum = length(UpperAgeBound)

#Compute lower age boundaries from UpperAgeBound input
#First age class begins at 0.
LowerAgeBound = [0;UpperAgeBound[1:end-1].+1]

#Intiialise and populate HealthEpsNonFatalQALYLoss array
HealthEpsParamByAge = zeros(SingleYrAgeGrpNum,EntityNum,StrainNum)
for AgeBandIdx = 1:AgeBandNum

    #Based on age bands, get indexes to access HealthEpsLhood at
    StartArrayAccessIdx = LowerAgeBound[AgeBandIdx] + 1
    EndArrayAccessIdx = UpperAgeBound[AgeBandIdx] + 1
    EntriesInCurrentAgeBand = EndArrayAccessIdx - StartArrayAccessIdx + 1

    #Input value across all  strains
    HealthEpsParamByAge[StartArrayAccessIdx:EndArrayAccessIdx,:,:] =
            repeat(copy(selectdim(ParamVals,1,[AgeBandIdx])),outer=(EntriesInCurrentAgeBand,1,1))
end
return HealthEpsParamByAge
end
