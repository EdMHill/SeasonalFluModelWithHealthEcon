#Purpose:
#Function to update exposure history array values
#--------------------------------------------------------------------------

function ExpHistUpdate_FM(ExpHistVaccType,ExpHistArray,ExpHistArrayParams,VaccEfficacy,M,
                                            ExpHistNum,NumOfStrains,AgeGrpSuscepTuple)

#Inputs
# ExpHistVaccType - Flag variable to specify how susceptibility will be modified for
#                   vaccine-related exposure history classes
#                  vacc. groups is reduced (0 - Not reduced, 1 - Reduced)
# ExpHistArray -
# ExpHistArrayParams - additional ExpHistArray relevant paramter inputs as extra cell entries
# VaccEfficacy - For last season, effectiveness of vaccine
# M - number of age classes in use
# ExpHistNum - Total number of exposure histories
# NumOfStrains
# AgeGrpSuscepTuple - Cell: Three entries
#                   --> Cell 1 - AgeGrpSuscepParamNum
#                   --> Cell 2 - AgeSuscepLowerBounds
#                   --> Cell 3 - AgeSuscepUpperBounds

#Outputs
# ExpHistArray - Updated version of exposure history array

if ExpHistVaccType == 1

    #Values to access particular array indexes
    HalfExpHistNum = convert(Int64,ExpHistNum/2)

    #Assign previous season efficacy scaling value to new variable name
    ExpHistVaccScaling = ExpHistArrayParams[3]

    #Update entry in ExpHistArray corresponding to vaccianted in previous season
    ExpHistVaccScaling_TempArray_0 = repeat(ExpHistVaccScaling*VaccEfficacy, outer = (1, 1, 1))
    ExpHistVaccScaling_TempArray_1 = repeat(ExpHistVaccScaling_TempArray_0, outer = (1, 1, HalfExpHistNum))
    ExpHistVaccScaling_TempArray_2 = permutedims(ExpHistVaccScaling_TempArray_1,[2,3,1])

    ExpHistArray[:,HalfExpHistNum+1:end,:] = ones(NumOfStrains,HalfExpHistNum,M) - ExpHistVaccScaling_TempArray_2

    for ii = 2:HalfExpHistNum
        #Update entries: natural infection
        ExpHistArray[ii-1,ii+HalfExpHistNum,:] = min.(ExpHistArrayParams[1],ExpHistArray[ii-1,ii+HalfExpHistNum,:])
    end

    #Update entries: influenza B corss-reactivity
    ExpHistArray[NumOfStrains-1,end,:] = min.(ExpHistArrayParams[2],ExpHistArray[NumOfStrains-1,end,:])
    ExpHistArray[NumOfStrains,end-1,:] = min.(ExpHistArrayParams[2],ExpHistArray[NumOfStrains,end-1,:])

elseif ExpHistVaccType == 2 #Age-dependent exposure history parameters

    #Disaggregate AgeGrpSuscepTuple
    AgeGrpSuscepParamNum = AgeGrpSuscepTuple[1]::Int64
    AgeSuscepLowerBounds = AgeGrpSuscepTuple[2]::Array{Int64,1}
    AgeSuscepUpperBounds = AgeGrpSuscepTuple[3]::Array{Int64,1}

    #Values to access particular array indexes
    HalfExpHistNum = convert(Int64,ExpHistNum/2)

    #Assign previous season efficacy scaling value to new variable name
    ExpHistVaccScaling = ExpHistArrayParams[:,3]

    #Iterate through each age band (groupings match those used for susceptibility parameters)
    for jj = 1:AgeGrpSuscepParamNum

        #Array indexing variables
        StartIdx = AgeSuscepLowerBounds[jj] + 1
        EndIdx = AgeSuscepUpperBounds[jj] + 1

        #Right hand columns, vaccinated previous season
        ExpHistArray[:,HalfExpHistNum+1:end,StartIdx:EndIdx] .= ExpHistArrayParams[jj,3]

        #Update entry in ExpHistArray corresponding to vaccianted in previous season
        ExpHistVaccScaling_TempArray_0 = repeat(ExpHistVaccScaling[jj]*VaccEfficacy[StartIdx:EndIdx,:], outer = (1, 1, 1))
        ExpHistVaccScaling_TempArray_1 = repeat(ExpHistVaccScaling_TempArray_0, outer = (1, 1, HalfExpHistNum))
        ExpHistVaccScaling_TempArray_2 = permutedims(ExpHistVaccScaling_TempArray_1,[2,3,1])

        ExpHistArray[:,HalfExpHistNum+1:end,StartIdx:EndIdx] = ones(NumOfStrains,HalfExpHistNum,EndIdx-StartIdx+1) - ExpHistVaccScaling_TempArray_2

        for ii = 2:HalfExpHistNum
            #Update entries: natural infection
            ExpHistArray[ii-1,ii+HalfExpHistNum,StartIdx:EndIdx] = min.(ExpHistArrayParams[jj,1],ExpHistArray[ii-1,ii+HalfExpHistNum,StartIdx:EndIdx])
        end

        #Update entries: influenza B corss-reactivity
        ExpHistArray[NumOfStrains-1,end,StartIdx:EndIdx] = min.(ExpHistArrayParams[jj,2],ExpHistArray[NumOfStrains-1,end,StartIdx:EndIdx])
        ExpHistArray[NumOfStrains,end-1,StartIdx:EndIdx] = min.(ExpHistArrayParams[jj,2],ExpHistArray[NumOfStrains,end-1,StartIdx:EndIdx])
    end
end

return ExpHistArray::Array{Float64,3}

end
