#Purpose:
#Functions to run model with leaky vaccine

#Model specifics:
#Influenza deterministic transmission dynamic model
#Multi-strain, age structured
#Vaccination status & exposure history (for previous season only) included
#Demography not explicitly included. Population distribution updated at end of season
#-------------------------------------------------------------------------------

# Define ODE equations
function SeasonalFlu_ODEModelFMAlt_LeakyVacc(dPop,pop,p,t)
     #Inputs of form (dPop,Pop,p,t)
     # dPop - rate of change for each entity
     # pop - current value for each entity
     # t - current time

     # p comprises following    # AgeRiskParam - Vector: M - Number of age classes; RiskGrpNum - Number of risk groups
    # ContactArray - Contact matrix
    # infection_param - Vector: beta - Transmission rate, sigma - rate of loss of latency
    # gamma - rate of loss of infectiousness
    # AgeSuscep - age-dependent susceptibility
    # vacc_per_day - rate of immunisation (age specific)
    # LeakyVaccVar - Two column array. Row i for age group i
    #                   Column 1 - Vaccination efficacy (on susceptibility),
    #                   Column 2 - Infectiousness reduction for
    #                              those vaccinated
    # LeakyTransFlag - Flag variable to specify if infectiousness of infected
    #                  vacc. groups is reduced (0 - Not reduced, 1 - Reduced)    # DayofYear - Actual day of year, with 0-1 corresponding to Jan 1st, 1-2 Jan 2nd etc
    # ExpHistArray - Interaction between exposure history and susceptibility to
    #                  the current season strain variant
    # NumOfStrains - Number of strains

    #Outputs
    # dPop_vec - Return values of derivatives at time t for values pop

    #Disaggregate p inputs
    DemogParam = p[1]::Array{Float64}
    ContactArray= p[2]::Array{Float64,2}
    InfectionParam = p[3]::Array
    AgeSuscep = p[4]::Array{Float64,2}
    vacc_per_day = p[5]::Array{Float64,2}
    LeakyVaccVar = p[6]::Array
    LeakyTransFlag = p[7]::Int64
    DayOfYear = t + p[8]::Float64
    ExpHistArray = p[9]::Array{Float64,3}
    NumOfStrains = p[10]::Int64

    #Disaggregate Age and Risk class number parameters
    M = DemogParam[1]
    M = convert(Int64,M) #Assert as integer
    TotalPopn = DemogParam[2]::Float64

    #Disaggregate infection param
    beta = InfectionParam[1]::Array{Float64,1} #Transmission rate
    sigma = InfectionParam[2]::Array{Float64,1} #rate of loss of latency
    gamma = InfectionParam[3]::Array{Float64,1} #rate of loss of infectiousness

    #Get current day of year, may be used for seasonality, vaccine rate etc
    currentday = Int(ceil(mod(DayOfYear,365))) #Cast as integer, otherwise will be float & throw error
    if currentday == 0
        currentday = 365 #Deal with case where at precise end of calendar year
    end
    v = vacc_per_day[:,currentday] #Get vaccine rate for current day

    #Get number of exposure histories to be accounted for
    ExpHistNum = size(ExpHistArray,2)::Int64

    #--------------------------------------------------------------------------
    ### PUT SUSCEPTIBLE COMPARTMENT VALUES INTO ARRAY
    SusClassTotal = M*ExpHistNum #per risk group, number of equations for S^N or S^V
    S_NotV = reshape(pop[1:SusClassTotal],M,ExpHistNum)
    S_V = reshape(pop[SusClassTotal+1:SusClassTotal*2],M,ExpHistNum)

    #--------------------------------------------------------------------------
    ### PUT I, E, R COMPARTMENT VALUES INTO ARRAY
    #Initialise numbers in each class and ODE variables
    SuscepClassEndIdx = SusClassTotal*2
    GroupsPerNonSusState = M*NumOfStrains
    E_NotV = reshape(pop[SuscepClassEndIdx+1:SuscepClassEndIdx+GroupsPerNonSusState],M,NumOfStrains)
    I_NotV = reshape(pop[SuscepClassEndIdx+(1*GroupsPerNonSusState)+1:SuscepClassEndIdx+(2*GroupsPerNonSusState)],M,NumOfStrains)
    R_NotV = reshape(pop[SuscepClassEndIdx+(2*GroupsPerNonSusState)+1:SuscepClassEndIdx+(3*GroupsPerNonSusState)],M,NumOfStrains)
    E_V = reshape(pop[SuscepClassEndIdx+(3*GroupsPerNonSusState)+1:SuscepClassEndIdx+(4*GroupsPerNonSusState)],M,NumOfStrains)
    I_V = reshape(pop[SuscepClassEndIdx+(4*GroupsPerNonSusState)+1:SuscepClassEndIdx+(5*GroupsPerNonSusState)],M,NumOfStrains)
    R_V = reshape(pop[SuscepClassEndIdx+(5*GroupsPerNonSusState)+1:SuscepClassEndIdx+(6*GroupsPerNonSusState)],M,NumOfStrains)

    #Initialise ODE indexing variables
    SusODENum = 2*ExpHistNum*M
    ODENumPerNonSusGrp = NumOfStrains*M
    #--------------------------------------------------------------------------
    ### INITIALISE ODE VARIABLES
    dPop_SusNotV = zeros(M,ExpHistNum)
    dPop_SusV = zeros(M,ExpHistNum)
    dPop_ENotV = zeros(M,NumOfStrains)
    dPop_INotV = zeros(M,NumOfStrains)
    dPop_RNotV = zeros(M,NumOfStrains)
    dPop_EV = zeros(M,NumOfStrains)
    dPop_IV = zeros(M,NumOfStrains)
    dPop_RV = zeros(M,NumOfStrains)
    dPop_C = zeros(M,NumOfStrains)

    #Combine susceptible compartment sizes by age group
    S = sum(S_NotV + S_V,dims=2) #Sum across exposure histories
    S_NotVTotal=sum(S_NotV,dims=2)

    #Combine non-vaccinated and vaccinated class popn sizes
    E = sum(E_NotV + E_V,dims=2); I_total = sum(I_NotV + I_V,dims=2);
    R = sum(R_NotV + R_V,dims=2);

    #Get proportion of population unvaccinated, stratified by age
    UnvaccPropn = S_NotVTotal + sum(E_NotV,dims=2) + sum(I_NotV,dims=2) + sum(R_NotV,dims=2)

    #---------------------------------------------------------------------------
    ### CALCULATE BIRTH RATE
    n = S+E+I_total+R #Population in each age class

    propn_NotV_popn_ByAge = UnvaccPropn./n

    #---------------------------------------------------------------------------
    ### RATE OF IMMUNISATION
    mu = v./propn_NotV_popn_ByAge[:]
    mu[UnvaccPropn[:].<=0] .= 0.0

    #--------------------------------------------------------------------------
    ### LEAKY VACCINE VARIABLES
    #Disaggragte LeakyVaccVar, scale force of infection if appropriate
    if LeakyTransFlag == 0 #Unmodfiied infectiousness
        alpha = LeakyVaccVar[1]::Array{Float64,2}
        #Get total number of infecteds. Array -> row per age, column per strain
        I = I_NotV + I_V
    elseif LeakyTransFlag == 1 #Reduced infectiousness active
        alpha = LeakyVaccVar[1]::Array{Float64,2}
        delta = LeakyVaccVar[2]::Array{Float64,2}

        #Get scaled force of infection. Array -> row per age, column per strain
        I = zeros(M,NumOfStrains)
        for ii = 1:M
             for jj = 1:NumOfStrains
                  I[ii,jj] = I_NotV[ii,jj] + (1.0-delta[ii,jj])*I_V[ii,jj]
             end
        end
        #I = I_NotV + (1.0-delta).*I_V #Vectorised version
    else
        error("Invalid value of $LeakyTransFlag for LeakyTransFlag, should be 0 or 1")
    end

    #---------------------------------------------------------------------------
    ### ODE COMPUTATIONS
    #Evaluate the infectious contacts per age group
    InfContacts = (ContactArray*(I./(n)))::Array{Float64,2} #Vectorised version

    # For each age class, calculate the derivative for each infection compartment
    for i=1:M
            #---------------------------------------------------------------
            ### Iterate through each set of exposure history and
            ### associated susceptible compartments
            for j = 1:ExpHistNum

                 #Population-level susceptibility against all strains for exposure hist. grp j,
                 #non-vaccinated, age group i
                 Sum_ModSusStrain_NotV = 0.
                 for kk = 1:NumOfStrains
                     Sum_ModSusStrain_NotV += AgeSuscep[i,kk]*beta[kk]*InfContacts[i,kk]*ExpHistArray[kk,j,i]
                 end
                 #Equivalent to sum(AgeSuscep[i,:].*beta.*InfContacts[i,:].*ExpHistArrayTranspose[j,:])

                #Non-vaccinated population, S compartments
                dPop[i+((j-1)*M)] = -Sum_ModSusStrain_NotV*S_NotV[i,j] - mu[i]*S_NotV[i,j] #Susceptible
                #dPop[i+((j-1)*M)] = -sum(AgeSuscep[i,:].*beta.*InfContacts[i,:].*ExpHistArrayTranspose[j,:])*S_NotV[i,j] - mu[i]*S_NotV[i,j] #Susceptible

                #Population-level susceptibility against all strains for exposure hist. grp i,
                #Vaccinated
                Sum_ModSusStrain_V = 0.
                for kk = 1:NumOfStrains
                    Sum_ModSusStrain_V += AgeSuscep[i,kk]*(1.0-alpha[i,kk])*beta[kk]*InfContacts[i,kk]*ExpHistArray[kk,j,i]
                end
                #Equivalent to sum(AgeSuscep[i,:].*beta.*(1.0-alpha[i,:]).*InfContacts[i,:].*ExpHistArrayTranspose[j,:])

                #Vaccinated population, S compartments
                dPop[(ExpHistNum*M)+((j-1)*M)+i] = -Sum_ModSusStrain_V*S_V[i,j] + mu[i]*S_NotV[i,j] #Susceptible
                #dPop[(ExpHistNum*M)+((j-1)*M)+i] = -sum(AgeSuscep[i,:].*beta.*(1.0-alpha[i,:]).*InfContacts[i,:].*ExpHistArrayTranspose[j,:])*S_V[i,j] + mu[i]*S_NotV[i,j] #Susceptible

            end

           #---------------------------------------------------------------
           ### ITERATE THROUGH EACH E,I,R COMPARTMENT
           for k = 1:NumOfStrains

            #---------------------------------------------------------------
            ### Non-vaccinated population, E, I, R compartments

            #Population-level susceptibility for age class i, strain k, non-vaccinated
            Sum_ModSus_NotV = 0.
            for jj = 1:ExpHistNum
                Sum_ModSus_NotV += ExpHistArray[k,jj,i]*S_NotV[i,jj] #Get modified susceptibility value for each exposure history grp
            end

            dPop[SusODENum+((k-1)*M)+i] = AgeSuscep[i,k]*beta[k]*InfContacts[i,k]*Sum_ModSus_NotV -sigma[k]*E_NotV[i,k] - mu[i]*E_NotV[i,k]  #Latent

            dPop[SusODENum+1*ODENumPerNonSusGrp+((k-1)*M)+i] = sigma[k]*E_NotV[i,k] - gamma[k]*I_NotV[i,k] - mu[i]*I_NotV[i,k]  #Infectious
            dPop[SusODENum+2*ODENumPerNonSusGrp+((k-1)*M)+i] = gamma[k]*I_NotV[i,k] - mu[i]*R_NotV[i,k] #Recovered/Immune

            #---------------------------------------------------------------
            ### Vaccinated population, E, I, R compartments

            #Population-level susceptibility for age class i, strain k, non-vaccinated
            Sum_ModSus_V = 0.
            for jj = 1:ExpHistNum
               Sum_ModSus_V += ExpHistArray[k,jj,i]*S_V[i,jj] #Get modified susceptibility value for each exposure history grp
            end

            dPop[SusODENum+3*ODENumPerNonSusGrp+((k-1)*M)+i] = AgeSuscep[i,k]*beta[k]*(1.0-alpha[i,k])*InfContacts[i,k]*Sum_ModSus_V -sigma[k]*E_V[i,k] + mu[i]*E_NotV[i,k]  #Latent

            dPop[SusODENum+4*ODENumPerNonSusGrp+((k-1)*M)+i] = sigma[k]*E_V[i,k] - gamma[k]*I_V[i,k] + mu[i]*I_NotV[i,k]  #Infectious
            dPop[SusODENum+5*ODENumPerNonSusGrp+((k-1)*M)+i] = gamma[k]*I_V[i,k] + mu[i]*R_NotV[i,k] #Recovered/Immune

            #---------------------------------------------------------------------------
            # ODE for cumulative number of cases
            dPop[SusODENum+6*ODENumPerNonSusGrp+((k-1)*M)+i] = sigma[k]*(E_NotV[i,k]+E_V[i,k])
          end
    end
end

function  PopulateContactArray(ContactArrayTemp)
#Inputs
    # ContactArrayTemp - Array of average number of adequate contacts. Ages 90+ not yet modified
#Outputs
    # ContactArrayTransformed - Array of average number of adequate contacts WITH Ages 90+ MODIFIED

    #Initialise array
    ContactArrayTransformed = ContactArrayTemp

    #Ages 0-90 interacting with ages 90+. Match value for contacts with 89-90
    ContactArrayTransformed[1:90,91:end] = ones(90,11).*ContactArrayTemp[1:90,90]

    #Ages 90+ interacting with ages 0-90. Match values of contacts for 89-90 aged individual with ages 0-90
    ContactArrayTransformed[91:end,1:90] = ones(11,90).*ContactArrayTemp[90,1:90]' #Transpose from column vector to row vector

    #Interatvions between 90+ age groups. Match 89-90 age bracket.
    ContactArrayTransformed[91:end,91:end] .= ContactArrayTemp[90,90]

return ContactArrayTransformed
end

function RunSeasonalFluModelFMAlt(ContactArray,MortalityFile,ONSPopnDistEngland,
    SimnParam,InfectionParam,AgeSuscep,AgeGrpSuscepTuple,
    VaccUptakeBySeason,LeakyVaccVarBySeason,LeakyTransFlag,
    ExpHistArrayFnInputs,ExpHistNum,ExpHistVaccType,MultiSeasonImmPropn,
    InfPropn_StartOfSeason,ICFromFile,SimnRunType)
#Inputs
# ContactArray - contact data for age class mixing
# MortalityFile, ONSPopnDistEngland - Import mortality rates and
# age distributions
# SimnParam - Specify month of year simulation will begin (1-Jan, 2-Feb,..., 12-Dec),
#             simn length and ode45 timestep value
# infection_param - Vector: beta - Transmission rate, sigma - rate of loss of latency, gamma - rate of loss of infectiousness
# AgeSuscep - Array: Age & strain specific susceptibilities (row per age class, column per strain)
# AgeGrpSuscepTuple - Cell: Three entries
#                   --> Cell 1 - AgeGrpSuscepParamNum
#                   --> Cell 2 - AgeSuscepLowerBounds
#                   --> Cell 3 - AgeSuscepUpperBounds
# BirthParam/DeathParam - Birth&death parameter
# VaccUptakeBySeason - rate of immunisation, daily uptake rates
# LeakyVaccVarBySeason -
#                   Two entry cell.
#                   --> Cell 1 - Vaccination efficacy (on susceptibility),
#                   --> Cell 2 - Infectiousness reduction for those vaccinated
#                  In each, multiple cells, one per season:
#                  Within each season cell, Two dimensional array. Row per age, column per strain.
# LeakyTransFlag - Flag variable to specify if infectiousness of infected
#                  vacc. groups is reduced (0 - Not reduced, 1 - Reduced)# ExpHistArray - Interaction between exposure history and susceptibility to
#                  the current season strain variant
# ExpHistArrayFnInputs - Two entry cell.
#                   --> Cell 1 - 2D array. Column per exposure history. Row per strain.
#                Interaction between exposure history and susceptibility to
#                  the current season strain variant.
#                  --> Cell 2 - Exposure history params used to populate
#                  ExpHistArray
# ExpHistNum - Number of exposure history classes in use
# ExpHistVaccType - Exposure history susceptibility modification form for vaccine-related states# InfPropn_StartOfSeason - % of population who are initialised as infected at
#                       beginning of each flu season
# MultiSeasonImmPropn - Constant: Proportion of population that are susceptible at end of season within
#                       a natural infection exposure history class that remain in that expsoure
#                       history category (rather than move to the Naive exposure history class)
# ICFromFile - Indicate whether initial conditions will be obtained from
#               external file and, if so, the file name
# SimnRunType - Type of run: 1 - exploratory; 2 - historical; 3 - alternative vacc. schemes

#Outputs
# Compartment counts at each timestep

#------------------------------------------------------------------------------
###  LOAD CONTACT DATA (NOW AN INPUT TO THIS FUNCTION!)

#------------------------------------------------------------------------------
### IMPORT MORTALITY RATES (NOW AN INPUT TO THIS FUNCTION!)

#------------------------------------------------------------------------------
### IMPORT INITIAL AGE DISTRIBUTION (NOW AN INPUT TO THIS FUNCTION!)

#------------------------------------------------------------------------------
### DEMOGRAPHY SIMN PARAMETERS

#Use first row of ONSPopnDistEngland data array
#Corresponds to 2010 absolute population values stratified by age
PopnDist = ONSPopnDistEngland[1,:]
AbsolutePopn = ONSPopnDistEngland[1,:]

#Assign settled age profile (at start of simulation year) to variable
current_n = PopnDist/sum(PopnDist)

#Number of age classes
M = length(current_n)::Int64
TotalPopn = 6e7
DemogParam = [M,TotalPopn]::Array
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

#------------------------------------------------------------------------------
### DISEASE DYNAMICS SIMN/FLAG PARAMETERS

#Disaggregate simulation paramter input vector
SimnStartDate = convert(Int64,SimnParam[1])
ODEBurnIn = SimnParam[2]
ODEStaticPopnTime = SimnParam[3]
ODEInferenceTime = SimnParam[4]
ODEForwardSimnTime = SimnParam[5]
MaxTime = SimnParam[6]
ODESampleTime = SimnParam[7]
timestep = SimnParam[8]
StoreFlag_PopnFOI = SimnParam[9]

#Values to pass to tspan in ODE fn
DaysPerMonth=[31 28 31 30 31 30 31 31 30 31 30 31]
CumulDaysCount=zeros(length(DaysPerMonth),1)
CumulDaysCount[2:end]=cumsum(DaysPerMonth[1:end-1])

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
###  TRANSMISSION RELATED PARAMETERS - CALCULATE BETA VALUES

#Disaggreagte InfectionParam
gamma = InfectionParam[3]::Array{Float64,1}
R_0 = InfectionParam[4]::Array{Float64,1}

#Transform ContactArray. Divide by population structure
#Divide row ii of ContactArray by current_n(ii)
ContactArrayTemp = ContactArray./current_n #Current_n column vector, orientation okay!

#Set values for 90+ age categories
#Match 89-90 single year age bracket values
ContactArrayTransformed = PopulateContactArray(ContactArrayTemp)

#Get spectral radius of transformed contact array.
spec_rad = maximum(abs.(eigvals(ContactArrayTransformed)))::Float64 #the spectral radius is the maximum of the absolute values of the eigenvalues

#Compute transmission rates to get desired basic reproduction number
#Assign to first tuple of InfectionParam
beta = gamma.*R_0/spec_rad
InfectionParam[1] = beta::Array{Float64,1}

#------------------------------------------------------------------------------
### EXPOSURE HISTORY PARAMETERS (NOW AN INPUT TO THIS FUNCTION!)

#------------------------------------------------------------------------------
###  DEMOGRAHPY PARAMETERS

#Import morality rates
MortalityProbByAge = readdlm(MortalityFile,',')

#Get number of strains to be accounted for
NumOfStrains = size(InfPropn_StartOfSeason,2)

#------------------------------------------------------------------------------
###  VACCINATION PARAMETERS (NOW AN INPUT TO THIS FUNCTION!)

#------------------------------------------------------------------------------
###  SET UP INITIAL CONDITIONS

# #Set up initial infected
# I0_NotV = zeros(M,NumOfStrains)
# for i = 1:NumOfStrains
#      I0_NotV[:,i] = current_n*InfPropn_StartOfSeason[i]
# end

#Set up initial infected, pick out specified row from SeasonStartPropnInf
TotalTime_H1N1Only = ODEBurnIn + ODEStaticPopnTime

#Row 1 for when only H1N1 infection allowed
#Row 2 when infection by any strain allowed
I0_NotV = zeros(M,NumOfStrains)
if TotalTime_H1N1Only > 0 #Initially, only H1N1 strain infection allowed
    for ii = 1:NumOfStrains
        I0_NotV[:,ii] = current_n*InfPropn_StartOfSeason[1,ii]
    end
else
    for ii = 1:NumOfStrains
         I0_NotV[:,ii] = current_n*InfPropn_StartOfSeason[2,ii]
    end
end

#Set initial conditions for ODEs
#Set up IC for susceptible groups, dependent on exposure history
ICFromFileFlag = ICFromFile[1]::Array{Int64,1} #Specify if initial conditions should be obtained from file
if ICFromFileFlag == 1
    vars = matread(ICFromFile[2]) #IF POSSIBLE, IMPORT WITH INITIAL INF PROPN ALREADY ACCOUNTED FOR!
    S0_NotV = vars["S0_NotV"]
#     RelMassPerExpHist=MassPerExpHist./sum(MassPerExpHist(:)) %Density in each age/exposure history Sus. compartment
#     j=1
#     while j<=ExpHistNum
#         #Modify susceptibles assigned to exposure history j
#         #Split initial infectious weighted by exposure history & age distribution
#         S0_NotV[j]=S0_NotV[j]-(RelMassPerExpHist[j,:].*I0_NotV)
#         j=j+1
#     end
else #All suscep. are put in naive, no exposure in previous year class
    S0_NotV = zeros(M,ExpHistNum)
    S0_NotV[:,1] = current_n - sum(I0_NotV,dims=2) #All susceptibles assigned to naive group
end

#All other compartments begin empty
S0_V = zeros(M,ExpHistNum)
E0_NotV=zeros(M,NumOfStrains); R0_NotV=zeros(M,NumOfStrains);
E0_V=zeros(M,NumOfStrains); I0_V=zeros(M,NumOfStrains); R0_V=zeros(M,NumOfStrains);
C0 = zeros(M,NumOfStrains)

#Initialise time and storage arrays for each set of compartments
T0=0.0
NumYrsRecorded = ceil(ODESampleTime/365)
CompartmentTraceSize = (ODESampleTime/timestep) + NumYrsRecorded #Add NumYrsRecorded to include initial conditions at start of each flu season
CompartmentTraceSize=convert(Int64,CompartmentTraceSize)

NumYrsSimTotal = ceil(MaxTime/365)
FOI_CompartmentTraceSize = (MaxTime/timestep) + 1 #Add one to include initial condition
FOI_CompartmentTraceSize = convert(Int64,FOI_CompartmentTraceSize)

Store_S_NotV= zeros(CompartmentTraceSize,M,ExpHistNum) #3rd dimension for exposure history
Store_S_V = zeros(CompartmentTraceSize,M,ExpHistNum) #3rd dimension for exposure history
Store_E_NotV = zeros(CompartmentTraceSize,M,NumOfStrains)
Store_I_NotV = zeros(CompartmentTraceSize,M,NumOfStrains)
Store_R_NotV = zeros(CompartmentTraceSize,M,NumOfStrains)
Store_E_V = zeros(CompartmentTraceSize,M,NumOfStrains)
Store_I_V = zeros(CompartmentTraceSize,M,NumOfStrains)
Store_R_V = zeros(CompartmentTraceSize,M,NumOfStrains)
Store_C = zeros(CompartmentTraceSize,M,NumOfStrains) #Cumulative number of cases

#Initialise storage array for population age distribution
Store_PopnDist = zeros(CompartmentTraceSize,M)

#Initialise storage array for time
Store_T = zeros(CompartmentTraceSize)

#If required, initialise storage array for population-level FOI
if StoreFlag_PopnFOI == 1
   Store_PopnFOI = zeros(FOI_CompartmentTraceSize,M,NumOfStrains)
end

#---------------------------------------------------------------------------
### INITIALISE VACCINE UPTAKE AND EFFICACY VALUES (BASED ON SIMN. RUN TYPE)

if SimnRunType == 1
    #SYNTHETIC DATA, no amendments to uptake/efficacy data needed
    vacc_per_day = VaccUptakeBySeason
    LeakyVaccVar = [LeakyVaccVarBySeason]::Array{Float64,2}
elseif SimnRunType == 2
    #BASED ON HISTORICAL DATA
    if TotalTime_H1N1Only > 0
        #Initially use 2009/2010 data (first row of arrays)
        RecordedSeasonIdx = 1
    else
        #If no H1N1 only period, use 2010/2011 data (second row of vacc. arrays)
        RecordedSeasonIdx = 2
    end

    vacc_per_day = VaccUptakeBySeason[RecordedSeasonIdx]::Array{Float64,2}

    if LeakyTransFlag == 0
        LeakyVaccVar = [LeakyVaccVarBySeason[RecordedSeasonIdx]]::Array
    elseif LeakyTransFlag == 1
        alpha = LeakyVaccVarBySeason[1][RecordedSeasonIdx]::Array{Float64,2}
        delta = LeakyVaccVarBySeason[2][RecordedSeasonIdx]::Array{Float64,2}
        LeakyVaccVar = [alpha,delta]::Array
    end
elseif SimnRunType == 3
    #RUN ALTERNATIVE VACC. SCHEMES
    if TotalTime_H1N1Only > 0
        #Initially use 2009/2010 data (first row of arrays)
        RecordedSeasonIdx = 1
    else
        RecordedSeasonIdx = 2
    end

    vacc_per_day = VaccUptakeBySeason[RecordedSeasonIdx]::Array{Float64,2}

    if LeakyTransFlag == 0
        LeakyVaccVar = [LeakyVaccVarBySeason[RecordedSeasonIdx]]::Array
    elseif LeakyTransFlag == 1
        alpha = LeakyVaccVarBySeason[1][RecordedSeasonIdx]::Array{Float64,2}
        delta = LeakyVaccVarBySeason[2][RecordedSeasonIdx]::Array{Float64,2}
        LeakyVaccVar = [alpha,delta]::Array
    end
else
    error("Unknown SimnRunType value, $SimnRunType")
end

#---------------------------------------------------------------------------
### CHECK CHOICE OF EXPOSURE HISTORY, UPDATE USING VACCINE EFFICACY IF NEEDED

#Pass to function
ExpHistArray = ExpHistArrayFnInputs[1]::Array{Float64,3}
ExpHistArrayParams = ExpHistArrayFnInputs[2]

ExpHistArray = ExpHistUpdate_FM(ExpHistVaccType,ExpHistArray,ExpHistArrayParams,LeakyVaccVar[1],M,
                                            ExpHistNum,NumOfStrains,AgeGrpSuscepTuple)

#------------------------------------------------------------------------------
### MAIN MODEL LOOP

#Get total number of ODEs in system
#ODE_EqnNumPerAge_Sus = ExpHistNum*2
ODE_EqnNum_Sus = 2*ExpHistNum*M
ODE_EqnNum_NonSus = 7*M*NumOfStrains
ODE_TotalEqnNum = (ODE_EqnNum_NonSus)+(ODE_EqnNum_Sus)

#Initialise initial condition vector
#OrigIC = zeros(Float64,ODE_TotalEqnNum)

#Specify current month
MonthIdx = SimnStartDate

#Initialise storage variables
InitialStoreFlag = 0 #Flag to specify if any values have been stored yet
StoreStartIdx = 1 #Initialise row indexing counter for storage arrays

InitialStoreFlag_PopnFOI = 0 #Flag to specify if any popnFOI values have been stored yet
StorePopnFOI_StartIdx = 1 #Initialise row indexing counter for popnFOI storage arrays]

#Initialise counter used to index row of ONS population distribution data array
PopnDistCounter = 1

while T0<MaxTime

    #tspan = (0.0,365.0)
    tspanForSol = (T0,T0+DaysPerMonth[MonthIdx]) #Specify time bounds to used in ODE solver
    tspan = T0:timestep:T0+DaysPerMonth[MonthIdx] #Timestep outputs

    #Concatenate suscepible class initial conditions
    IC =  [S0_NotV[:]; S0_V[:]; E0_NotV[:]; I0_NotV[:]; R0_NotV[:]; E0_V[:]; I0_V[:]; R0_V[:]; C0[:]]

    #Error check
    if any(isnan,IC)
        error("NaN found in initial conditions")
    end

    #Move negative values back to 0
    #IC[IC.<0.0] = 0.0

    #println("Sum OrigIC is:", sum(OrigIC[1:(end-M)]))
    #Normalise IC so it sums to one!
    IC[1:end-(M*NumOfStrains)] = (IC[1:end-(M*NumOfStrains)]/sum(IC[1:end-(M*NumOfStrains)]))::Array{Float64,1}

    #---------------------------------------------------------------------------
    ### RUN ODE SOLVER
    #Can now add parameters as input to ODE problem!
    p = [DemogParam, ContactArrayTransformed, InfectionParam, AgeSuscep,
            vacc_per_day, LeakyVaccVar, LeakyTransFlag,
            CumulDaysCount[SimnStartDate],ExpHistArray,NumOfStrains]

    #Code for checking type-related problems
    #du = similar(IC)
    #@code_warntype SeasonalFlu_ODEModelFMAlt_LeakyVacc(du,IC,p,10)

    #Define ODE problem and solve
    prob = ODEProblem(SeasonalFlu_ODEModelFMAlt_LeakyVacc,IC,tspanForSol,p)

    #can use alg_hints=[:stiff] if we have a stiff problem where we need high accuracy, but don't know the best stiff algorithm for this problem
    if InitialStoreFlag==0 || MonthIdx == SimnStartDate #If start of flu season, store initial value, save_start = true
        sol = solve(prob, Tsit5(), saveat=timestep, abstol = 1e-12, reltol = 1e-10,
       #callback = PositiveDomain(save=false),force_dtmin=true)
       isoutofdomain=(u,p,t) -> any(x -> x < 0, u))
    else
       sol = solve(prob, Tsit5(), saveat=timestep, save_start = false, abstol = 1e-12, reltol = 1e-10,
       #callback = PositiveDomain(save=false),force_dtmin=true)
       isoutofdomain=(u,p,t) -> any(x -> x < 0, u))
    end

    #---------------------------------------------------------------------------
    ### PERFORM STORAGE OF COMPARTMENT VALUES
    SusClassTotal = M*ExpHistNum::Int64 #number of equations for S^N or S^V
    NonSusClassTotal = M*NumOfStrains::Int64 #number of equations per E,I,R,C state

    # Amalgamate output in arrays for current time window
    T = tspan'

    # Collate into Array
   OutputNum=length(sol)::Int64
   pop = zeros(OutputNum,ODE_TotalEqnNum)
   for k = 1:OutputNum
       pop[k,:] = sol[k]::Array{Float64,1}
   end

    #Separate columns into susceptible and non-susceptible compartments
    SusNotVPop = reshape(pop[:,1:SusClassTotal]::Array{Float64,2},:,M,ExpHistNum)::Array{Float64,3}
    SusVPop = reshape(pop[:,SusClassTotal+1:SusClassTotal*2]::Array{Float64,2},:,M,ExpHistNum)::Array{Float64,3}
    ENotV_Pop = reshape(pop[:,SusClassTotal*2+1:(SusClassTotal*2) + NonSusClassTotal],:,M,NumOfStrains)::Array{Float64,3}
    INotV_Pop = reshape(pop[:,(SusClassTotal*2)+NonSusClassTotal+1:(SusClassTotal*2) + 2*NonSusClassTotal],:,M,NumOfStrains)::Array{Float64,3}
    RNotV_Pop = reshape(pop[:,(SusClassTotal*2)+(2*NonSusClassTotal)+1:(SusClassTotal*2) + 3*NonSusClassTotal],:,M,NumOfStrains)::Array{Float64,3}
    EV_Pop = reshape(pop[:,(SusClassTotal*2)+(3*NonSusClassTotal)+1:(SusClassTotal*2) + 4*NonSusClassTotal],:,M,NumOfStrains)::Array{Float64,3}
    IV_Pop = reshape(pop[:,(SusClassTotal*2)+(4*NonSusClassTotal)+1:(SusClassTotal*2) + 5*NonSusClassTotal],:,M,NumOfStrains)::Array{Float64,3}
    RV_Pop = reshape(pop[:,(SusClassTotal*2)+(5*NonSusClassTotal)+1:(SusClassTotal*2) + 6*NonSusClassTotal],:,M,NumOfStrains)::Array{Float64,3}
    C_Pop = reshape(pop[:,(SusClassTotal*2)+(6*NonSusClassTotal)+1:(SusClassTotal*2) + 7*NonSusClassTotal],:,M,NumOfStrains)::Array{Float64,3}

    #If beyond burn in, store compartment values
    if T[end]>ODEBurnIn

         #Get proportion of population in each age group at each timstep
        n_PopnDist = sum(SusNotVPop,dims=3) + sum(SusVPop,dims=3) + sum(ENotV_Pop,dims=3) +
                     sum(INotV_Pop,dims=3) + sum(RNotV_Pop,dims=3) + sum(EV_Pop,dims=3) +
                     sum(IV_Pop,dims=3) + sum(RV_Pop,dims=3)

         if InitialStoreFlag == 0
            StoreEndIdx = length(tspan)
            Store_T[StoreStartIdx:StoreEndIdx] = tspan'
            InitialStoreFlag = 1 #Amend flag value so will not pass through this loop again!
        elseif MonthIdx == SimnStartDate #Store initial condition if at start of flu season
            StoreEndIdx = StoreStartIdx + length(tspan) - 1
            #Need to realign indexing. Adding on length(tspan) makes index 1 too large

            Store_T[StoreStartIdx:StoreEndIdx] = tspan'
        else
            StoreEndIdx = StoreStartIdx + length(tspan) - 2
            #Subtract 2 as will not store initial condition,
            Store_T[StoreStartIdx:StoreEndIdx] = tspan[2:end]'
        end

            #Store susceptible compartment values
            Store_S_NotV[StoreStartIdx:StoreEndIdx,:,:] = SusNotVPop
            Store_S_V[StoreStartIdx:StoreEndIdx,:,:] =  SusVPop

            #Store non-susceptible compartment values
            Store_E_NotV[StoreStartIdx:StoreEndIdx,:,:] = ENotV_Pop
            Store_I_NotV[StoreStartIdx:StoreEndIdx,:,:] = INotV_Pop
            Store_R_NotV[StoreStartIdx:StoreEndIdx,:,:] = RNotV_Pop
            Store_E_V[StoreStartIdx:StoreEndIdx,:,:] = EV_Pop
            Store_I_V[StoreStartIdx:StoreEndIdx,:,:] = IV_Pop
            Store_R_V[StoreStartIdx:StoreEndIdx,:,:] = RV_Pop
            Store_C[StoreStartIdx:StoreEndIdx,:,:] = C_Pop

            #Store population distribution at each timestep
            Store_PopnDist[StoreStartIdx:StoreEndIdx,:] = n_PopnDist

            StoreStartIdx = StoreEndIdx + 1 #Update StartIdx value
    end

    #----------------------------------------------------------------------
    #Store total number of infected, can then be used to recover risk group
    #information!
    if StoreFlag_PopnFOI == 1
        if InitialStoreFlag_PopnFOI == 0
            StorePopnFOI_EndIdx = length(tspan)

            #Scale force of infection if appropriate
            if LeakyTransFlag == 0 #Unmodfiied infectiousness
                Store_PopnFOI[StorePopnFOI_StartIdx:StorePopnFOI_EndIdx,:,:] =
                    INotV_Pop + IV_Pop #Get total number of infecteds. Vector, entry per strain
            elseif LeakyTransFlag == 1 #Reduced infectiousness active
                delta = LeakyVaccVar[:,2]
                Store_PopnFOI[StorePopnFOI_StartIdx:StorePopnFOI_EndIdx,:,:] = INotV_Pop + (1-delta').*IV_Pop
                #Get scaled force of infection. Vector, entry per strain
            else
                error("Invalid value of $LeakyTransFlag for LeakyTransFlag, should be 0 or 1")
            end

            InitialStoreFlag_PopnFOI = 1 #Amend flag value so will not pass through this loop again!
         elseif MonthIdx == SimnStartDate #Store initial condition if at start of flu season
            StorePopnFOI_StartIdx = StorePopnFOI_StartIdx - 1 #Want to replace final reads from previus flu season with new IC
            StorePopnFOI_EndIdx = StorePopnFOI_StartIdx + length(tspan) - 1
            #Need to realign indexing. Adding on length(tspan) makes index 1 too large

            #Scale force of infection if appropriate
            if LeakyTransFlag == 0 #Unmodfiied infectiousness
                Store_PopnFOI[StorePopnFOI_StartIdx:StorePopnFOI_EndIdx,:,:] =
                    INotV_Pop + IV_Pop #Get total number of infecteds. Vector, entry per strain
            elseif LeakyTransFlag == 1 #Reduced infectiousness active
                delta = LeakyVaccVar[:,2]
                Store_PopnFOI[StorePopnFOI_StartIdx:StorePopnFOI_EndIdx,:,:] = INotV_Pop + (1-delta').*IV_Pop
                #Get scaled force of infection. Vector, entry per strain
            else
                error("Invalid value of $LeakyTransFlag for LeakyTransFlag, should be 0 or 1")
            end

        else
            StorePopnFOI_EndIdx = StorePopnFOI_StartIdx + length(tspan) - 2 #Subtract 2 as will not store initial condition,

            #Scale force of infection if appropriate
            if LeakyTransFlag == 0 #Unmodfiied infectiousness
                Store_PopnFOI[StorePopnFOI_StartIdx:StorePopnFOI_EndIdx,:,:] =
                    INotV_Pop + IV_Pop #Get total number of infecteds. Vector, entry per strain
            elseif LeakyTransFlag == 1 #Reduced infectiousness active
                delta = LeakyVaccVar[:,2]
                Store_PopnFOI[StorePopnFOI_StartIdx:StorePopnFOI_EndIdx,:,:] = INotV_Pop + (1-delta').*IV_Pop
                #Get scaled force of infection. Vector, entry per strain
            else
                error("Invalid value of $LeakyTransFlag for LeakyTransFlag, should be 0 or 1")
            end
        end
        StorePopnFOI_StartIdx = StorePopnFOI_EndIdx + 1 #Update StartIdx value
    end

    #---------------------------------------------------------------------------
    ### ASSIGN ODE MODEL OUTPUTS TO DESCRIPTIVE VARIABLE NAMES
    #Store susceptible compartment values
    S_NotV = reshape(SusNotVPop[end,:,:],:,ExpHistNum)::Array{Float64,2}
    S_V = reshape(SusVPop[end,:,:],:,ExpHistNum)::Array{Float64,2}

    #Store non-susceptible compartment columns to descriptive names
    E_NotV = reshape(ENotV_Pop[end,:,:],:,NumOfStrains)::Array{Float64,2}
    I_NotV = reshape(INotV_Pop[end,:,:],:,NumOfStrains)::Array{Float64,2}
    R_NotV = reshape(RNotV_Pop[end,:,:],:,NumOfStrains)::Array{Float64,2}
    E_V = reshape(EV_Pop[end,:,:],:,NumOfStrains)::Array{Float64,2}
    I_V = reshape(IV_Pop[end,:,:],:,NumOfStrains)::Array{Float64,2}
    R_V = reshape(RV_Pop[end,:,:],:,NumOfStrains)::Array{Float64,2}
    C = reshape(C_Pop[end,:,:],:,NumOfStrains)::Array{Float64,2}

    ###############################################################
    ### PERFORM END OF SEASON/AGE CLASS MOVEMENTS - BEGIN
    ###############################################################
    if MonthIdx == 8  #End of season, reinitialise infecteds
        #Get current age structure proportions (for risk group k)
        #First sum over the susceptible compartments
        S_NotV_Sum = sum(S_NotV,dims=2)
        S_V_Sum = sum(S_V,dims=2)
        S_Unvacc_MultiSeasonImmEligable = S_NotV[:,2:NumOfStrains+1] + S_NotV[:,NumOfStrains+3:end]
        S_Vacc_MultiSeasonImmEligable = S_V[:,2:NumOfStrains+1] + S_V[:,NumOfStrains+3:end]

        #---------------------------------------------------------------------------
        ### MAPPINGS TO EACH EXPOSURE HISTORY CLASS
        temp_S0_NotV = zeros(M,ExpHistNum) #Placeholder for Susceptible, unvaccinated initial condition

        #No exposure through infection/succesful vacc. (in effect, sum over S^N classes)
        temp_S0_NotV[:,1] = (S_NotV[:,1] + S_NotV[:,6]) + ((1-MultiSeasonImmPropn)*sum(S_Unvacc_MultiSeasonImmEligable,dims=2))

        #Natural infection by strain i, no vaccination
        temp_S0_NotV[:,2:5] = E_NotV + I_NotV + R_NotV + (MultiSeasonImmPropn*S_Unvacc_MultiSeasonImmEligable)

        #Vaccination, with no natural infection
        #Account for loss of "foldover" residual immunity
        #temp_S0_NotV[:,6] = S_V_Sum
        temp_S0_NotV[:,6] = (S_V[:,1] + S_V[:,6]) + ((1-MultiSeasonImmPropn)*sum(S_Vacc_MultiSeasonImmEligable,dims=2))

        #Natural infection by strain i, and vaccinated
        #Account for any "foldover" residual immunity
        temp_S0_NotV[:,7:10] = E_V + I_V + (MultiSeasonImmPropn*(S_V[:,2:5]+S_V[:,7:10]))

        S0_V = zeros(M,ExpHistNum) #Reset S^V class to zeros

        #-----------------------------------------------------------------------
        ### UPDATE POPULATION DISTRIBUTION (IF REQUIRED) - START
        #-----------------------------------------------------------------------

        #Add on non-susceptible compartments
        current_n = S_NotV_Sum + sum(E_NotV,dims=2) + sum(I_NotV,dims=2) +
            sum(R_NotV,dims=2) + S_V_Sum + sum(E_V,dims=2) + sum(I_V,dims=2) + sum(R_V,dims=2)

        #println("T0 is:", T0)

        if T[end] < ODEBurnIn + ODEStaticPopnTime #Using static population
            updated_n = current_n[:] #Ensure updated_n is vector, not nx1 array!
            S0_NotV = temp_S0_NotV

            #Youngest age class, all will be in "naive" susceptible class
            S0_NotV[1,1] = updated_n[1]
            S0_NotV[1,2:10] .= 0.0
        elseif T[end] < ODEBurnIn + ODEStaticPopnTime + ODEInferenceTime #Use historical popultion data, 2010-2016
            #Pick out specified row of ONS popn dist array
            PopnDistCounter = PopnDistCounter + 1

            #Get updated age partitioned population distribution
            AbsolutePopn = ONSPopnDistEngland[PopnDistCounter,:]
            updated_n = AbsolutePopn/sum(AbsolutePopn)

            #Find relative magnitude of next season population split per age compared to
            #current season
            scale_n = zeros(Float64,M-1)

            #Sclaing of next season 1-2 age group popn size to current season 0-1 age group popn size, and so on
            #scale_n[1:end-1] = updated_n[2:end-1]./current_n[1:end-2]
            for ii = 1:M-2
                 scale_n[ii] = updated_n[ii+1]/current_n[ii]
            end

            #Sclaing of next season 90+ age group popn size to current season combined 89-90 and 90+ group popn size
            scale_n[end] = updated_n[end]/sum(current_n[end-1:end])

            #Update susceptible and infectious class initial conditions
            S0_NotV[1,:] .= 0.0 #Ensure no mass in lowest age class for natural infection/vaccinated exposure histories
            S0_NotV[2:end-1,:] = temp_S0_NotV[1:end-2,:].*scale_n[1:end-1]

            #For 90+, eldest age class, take average across 89-90 class
            #and 90+ class
            S0_NotV[end,:] = (temp_S0_NotV[end-1,:] + temp_S0_NotV[end,:])*scale_n[end]

            S0_NotV[1,1] = updated_n[1] #Youngest age class, all will be in "naive" susceptible class
        else  #2017 onwards

            #Apply mortality probability to each age class
            #DeathAdjAbsolutePopn = AbsolutePopn.*(1-MortalityProbByAge)'
            DeathAdjAbsolutePopn = zeros(M)
            for ii = 1:M
                 DeathAdjAbsolutePopn[ii] = AbsolutePopn[ii]*(1-MortalityProbByAge[ii])
            end

            #Perform age class movements
            AbsolutePopn[2:end-1] = DeathAdjAbsolutePopn[1:end-2]
            AbsolutePopn[end] = DeathAdjAbsolutePopn[end-1] + DeathAdjAbsolutePopn[end]

            #Repopulate youngest age class
            AbsolutePopn[1] = 653467

            #Calculate updated age-stratified population distribution
            updated_n = AbsolutePopn/sum(AbsolutePopn)

            #Find relative magnitude of next season population split per age compared to
            #current season
            scale_n = zeros(Float64,M-1)

            #Sclaing of next season 1-2 age group popn size to current season 0-1 age group popn size, and so on
            #scale_n[1:end-1] = updated_n[2:end-1]./current_n[1:end-2]
            for ii = 1:M-2
                 scale_n[ii] = updated_n[ii+1]/current_n[ii]
            end

            #Sclaing of next season 90+ age group popn size to current season combined 89-90 and 90+ group popn size
            scale_n[end] = updated_n[end]/sum(current_n[end-1:end])

            #Update susceptible and infectious class initial conditions
            S0_NotV[1,:] .= 0.0 #Ensure no mass in lowest age class for natural infection/vaccinated exposure histories
            S0_NotV[2:end-1,:] = temp_S0_NotV[1:end-2,:].*scale_n[1:end-1]

            #For 90+, eldest age class, take average across 89-90 class
            #and 90+ class
            S0_NotV[end,:] = (temp_S0_NotV[end-1,:] + temp_S0_NotV[end,:])*scale_n[end]

            S0_NotV[1,1] = updated_n[1] #Youngest age class, all will be in "naive" susceptible class
        end

        #-----------------------------------------------------------------------
        ### UPDATE POPULATION DISTRIBUTION (IF REQUIRED) - END
        #-----------------------------------------------------------------------

        #-----------------------------------------------------------------------
        ### UPDATE TRANSMISSION RELATED PARAMETERS
        ### (USE REVISED POPN STRCUTURE) - BEGIN
        #-----------------------------------------------------------------------

        #gamma and R_0 values unchanged. Do not need to assign to variable again
        #gamma = InfectionParam[3]::Array{Float64,1}
        #R_0 = InfectionParam[4]::Array{Float64,1}

        #Transform ContactArray. Divide by population structure
        #Divide row ii of ContactArray by updated_n(ii)
        ContactArrayTemp = ContactArray./updated_n #updated_n column vector, orientation okay!

        #Set values for 90+ age categories
        #Match 89-90 single year age bracket values
        ContactArrayTransformed = PopulateContactArray(ContactArrayTemp)

        #Get spectral radius of transformed contact array.
        spec_rad = maximum(abs.(eigvals(ContactArrayTransformed)))::Float64 #the spectral radius is the maximum of the absolute values of the eigenvalues

        #Compute transmission rates to get desired basic reproduction number
        #Assign to first tuple of InfectionParam
        beta = gamma.*R_0/spec_rad
        InfectionParam[1] = beta::Array{Float64,1}

        #-----------------------------------------------------------------------
        ### UPDATE TRANSMISSION RELATED PARAMETERS
        ### (USE REVISED POPN STRCUTURE) - END
        #-----------------------------------------------------------------------

        #-----------------------------------------------------------------------
        ### Initialise infectious popn, update susceptible compartments
        # for i = 1:NumOfStrains
        #      I0_NotV[:,i] = updated_n*InfPropn_StartOfSeason[i]
        # end

        #Row 1 for when only H1N1 infection allowed
        #Row 2 when infection by any strain allowed
        if T[end] < TotalTime_H1N1Only  #Initially, only H1N1 strain infection allowed
            TotalSeasonStartPropnInf = sum(InfPropn_StartOfSeason[1,:]) #Proportion of popn initially infected
            for i = 1:NumOfStrains
                I0_NotV[:,i] = updated_n*InfPropn_StartOfSeason[1,i]
            end
        else
            TotalSeasonStartPropnInf = sum(InfPropn_StartOfSeason[2,:]) #Proportion of popn initially infected
            for i = 1:NumOfStrains
                 I0_NotV[:,i] = updated_n*InfPropn_StartOfSeason[2,i]
            end
        end

        RelMassPerExpHist = S0_NotV/sum(S0_NotV[:])

        #Modify susceptibles assigned to exposure history j
        #Split initial infectious weighted by exposure history & age distribution
        S0_NotV = S0_NotV - (RelMassPerExpHist*TotalSeasonStartPropnInf)

        #-----------------------------------------------------------------------
        ### Reset all IC for E and R non-vaccination compartments, and all vacc.compartments
        R0_NotV = zeros(M,NumOfStrains); E0_NotV = zeros(M,NumOfStrains)
        E0_V = zeros(M,NumOfStrains)
        I0_V = zeros(M,NumOfStrains)
        R0_V = zeros(M,NumOfStrains)

        #-----------------------------------------------------------------------
        ### Keep tracking cumulative number of infected, no reset required
        C0 = C

        #-----------------------------------------------------------------------
        ### IF BEYOND TotalTime_H1N1Only,
        #UPDATE EXPHISTARRAY, VACCINE UPTAKE AND EFFICACY VALUES (BASED ON SIMN. RUN TYPE)

        if TotalTime_H1N1Only == 0
            ExpArrayAmendTimeThreshold = ODEBurnIn
        else
            ExpArrayAmendTimeThreshold = TotalTime_H1N1Only
        end

        if T[end]>=ExpArrayAmendTimeThreshold && SimnRunType!=1
            #Note, SYNTHETIC DATA, no amendments to uptake/efficacy data needed

            #Get updated ExpHistArray
            ExpHistArray = ExpHistUpdate_FM(ExpHistVaccType,ExpHistArray,ExpHistArrayParams,LeakyVaccVar[1],M,
                                                        ExpHistNum,NumOfStrains,AgeGrpSuscepTuple)

            #Update Idx to pick out relevant season data from array
            if RecordedSeasonIdx < length(VaccUptakeBySeason)
                RecordedSeasonIdx = RecordedSeasonIdx + 1
            end

            if SimnRunType == 2 || SimnRunType == 3
                vacc_per_day = VaccUptakeBySeason[RecordedSeasonIdx]

                #For modified susceptibility/transmissibility due to leaky
                #vaccine, pick out column corresponding to RecordedSeasonIdx
                if LeakyTransFlag == 0
                    LeakyVaccVar = [LeakyVaccVarBySeason[RecordedSeasonIdx]]::Array
                elseif LeakyTransFlag == 1
                    alpha = LeakyVaccVarBySeason[1][RecordedSeasonIdx]::Array{Float64,2}
                    delta = LeakyVaccVarBySeason[2][RecordedSeasonIdx]::Array{Float64,2}
                    LeakyVaccVar = [alpha,delta]
                end
            else
                error("Unknown SimnRunType value, $SimnRunType")
            end
        end

    else #End of month, but not end of season. Perform age group movements only
        #-----------------------------------------------------------------------
        ### Non-susceptible compartment updates
        # Carry forward initial conditions from previous month for 2 years and greater
        E0_NotV = E_NotV
        I0_NotV = I_NotV
        R0_NotV = R_NotV
        E0_V = E_V
        I0_V = I_V
        R0_V = R_V

        #Carry forward cumulative infection counts for all age groups
        C0 = C

        #-----------------------------------------------------------------------
        ### Update susceptible compartments
        S0_NotV = S_NotV
        S0_V = S_V
    end
    ###############################################################
    ### PERFORM END OF SEASON/AGE CLASS MOVEMENTS - END
    ###############################################################

    T0 = T[end]

	#update month index
	MonthIdx = mod(MonthIdx + 1,12)
	if MonthIdx == 0
		MonthIdx = 12
	end
end

if StoreFlag_PopnFOI == 0
    return Store_T,Store_S_NotV,Store_S_V,Store_E_NotV,Store_E_V,Store_I_NotV,Store_I_V,
        Store_R_NotV,Store_R_V,Store_C,Store_PopnDist
else
    return Store_T,Store_S_NotV,Store_S_V,Store_E_NotV,Store_E_V,Store_I_NotV,Store_I_V,
        Store_R_NotV,Store_R_V,Store_C,Store_PopnDist,Store_PopnFOI
end

end
