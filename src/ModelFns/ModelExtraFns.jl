#Purpose:
#Contains functions to that may be used to
#perform Adaptive Population Monte Carlo Approximate Bayesian Computation
#-------------------------------------------------------------------------------


### GROUPING OF FUNCTIONS TO BE PASSED INTO THE APMC SCHEME ###
# (i) FUNCTION TO PERTURB SAMPLES WHEN GENERATING PROPOSED PARAMETER SETS
# (ii) FUNCTIONS TO DRAW SAMPLES FROM PRIOR
# (iii) SPECIFY PRIOR DENSITY
# (iv) SUMMARY STATISTIC CALCULATION FUNCTIONS
# (v) FUNCTIONS TO APPLY ASCERTAINMENT PROBABILITY TO INFECTION CASE LOAD
# (vi) FUNCTIONS TO CONSTRUCT SUSCEPTIBILITY ARRAY
# (vii) FUNCTION TO RUN MODEL SIMULATION & PRODUCE DESIRED OUTPUTS TO FEED INTO SUMMARY STATISTIC FUNCTION
# (viii) FUNCTIONS USED IN MODEL SIMULATION

#-------------------------------------------------------------------------------
# (i) FUNCTION TO PERTURB SAMPLES WHEN GENERATING PROPOSED PARAMETER SETS
#-------------------------------------------------------------------------------

function  FMAlt_OLCMPerturbVarFn(RetainedParams,RetainedWeights,SurvivorParticleIdx,GenT)
# Optimal local covariance matrix (OLCM)
# Uses a multivariate normal distribution with a covariance matrix based on a subset of the particles from the
# previous iteration, whose distances are smaller than the threshold of the current iteration

#Inputs:
#   RetainedParams,RetainedWeights - Current set of samples with associated weights
# 	SurvivorParticleIdx - Boolean true/false states defining the subset of the particles from the previous iteration, whose distances are smaller than the threshold of the current iteratio
#	  GenT - Generation number, used to determine if should sample using a global covariance
#Outputs:
#   C - Variance-covariance matrix. Three dimensional, slice per particle

#Define number of parameters and particles
RetainedParamsDims = size(RetainedParams)
ParamNum = RetainedParamsDims[2]  #Number of parameters matches number of columns in parameter array
ParticleNum = RetainedParamsDims[1] #Number of particles in generation matches number of rows in parameter array

#Update covariance arrays
if GenT == 0 #First generation,
	C_SingleSlice = 2.0*cov(RetainedParams,AnalyticWeights(RetainedWeights[:]),corrected=true)::Array{Float64}

	#Check if C is non-singular. Modify if singular.
	while rank(C_SingleSlice) != ParamNum #Check
		println("Before C eigvals: $(eigvals(C_SingleSlice))")
		C_SingleSlice = C_SingleSlice + 1e-12I
		println("After C eigvals: $(eigvals(C_SingleSlice))")
	end

	#Repeat covariance matrix, once per parameter set
	ParticleNum = size(RetainedParams,1) #Number of particles in generation matches number of rows in parameter array
	C = repeat(C_SingleSlice, outer=(1,1,ParticleNum))
else

	#Get parameter sets from the previous iteration, whose distances are smaller than the threshold of the current iteration
	# i.e. "those that survived", row ID of RetainedParams array dennoted by SurvivorParticleIdx
	SPP_thetas = RetainedParams[vec(SurvivorParticleIdx),:] 	#SPP, subset of previous particles

	#Pick out weights for subset of the particles from the
	# previous iteration, whose distances are smaller than the threshold of the current iteration
	# Normalise those weights
	SPP_Weights = RetainedWeights[vec(SurvivorParticleIdx)]
	SPP_NormWeights = SPP_Weights/sum(SPP_Weights)

	#Compute mean of the survivor particle population
    m = sum(SPP_thetas.*SPP_NormWeights,dims=1)
	#println(SPP_thetas)
	#println(m)

	#Initialise variance-covariance array
	C = zeros(ParamNum, ParamNum,ParticleNum)

	#Loop through each particle and compute covariance
	for kk = 1:1:ParticleNum
		Current_C = zeros(ParamNum, ParamNum)
		for jj = 1:1:ParamNum
			for ii = jj:1:ParamNum
				Current_C[ii, jj] = sum(SPP_NormWeights.*(SPP_thetas[:,ii] .- m[ii]).*(SPP_thetas[:, jj] .- m[jj])) + (m[ii] - RetainedParams[kk,ii])*(m[jj] - RetainedParams[kk,jj])

				#Covariance matrix is symmetric.
				#Assign transposed off-diagonal value
				if ii != jj
					Current_C[jj, ii] = Current_C[ii, jj]
				end
      end
		end


		#Check if Current_C is non-singular. Blow it up if singular.
		while rank(Current_C) != ParamNum #Check
			#println(abs(minimum(real(eigvals(Current_C)))))
			Current_C = Current_C + 1e-12I
		end

		#Assign revised covariance to C (collection of covariance arrays)
		C[:,:,kk] = Current_C
	end

	return C
end
end

function  FMAlt_ComponentWisePerturbVarFn(RetainedParams,RetainedWeights)
#Inputs:
#   RetainedParams,RetainedWeights - Current set of samples with associated weights

#Outputs:
#   C - Variance-covariance matrix

# if ndims(RetainedParams) > 1
#     C = PerturbVarScale*diag(cov(RetainedParams,Weights(RetainedWeights[:])))::Vector{Float64}
# else
#     C = PerturbVarScale*var(RetainedParams,Weights(RetainedWeights[:]))::Vector{Float64}
# end

    #Get variance-covariance array
	C = 0.1*diag(cov(RetainedParams,AnalyticWeights(RetainedWeights[:]),corrected=true))::Vector{Float64}

	return C
end

function  FMAlt_CovarPerturbVarFn(RetainedParams,RetainedWeights)
#Inputs:
#   RetainedParams,RetainedWeights - Current set of samples with associated weights

#Outputs:
#   C - Variance-covariance matrix

    #Get variance-covariance array
	C = 0.1*cov(RetainedParams,AnalyticWeights(RetainedWeights[:]),corrected=true)::Array{Float64}

	return C
end

#--------------------------------------------------------------------------
# (ii) FUNCTIONS TO DRAW SAMPLES FROM PRIOR
#--------------------------------------------------------------------------
function  FMAlt_SampleFirstGenFn_FromFileExample(N::Int64)
#Inputs:
#   N - Number of samples from LHC parameter set to take forward as first generation
    ParticlesFromPrior = readdlm("FMAlt_APMCOutputFilesJULIA/APMCsamples_FMAlt_JulRun#101_RetainedParams.txt",'\t')

    #Assign weights to particle
    Particle_Weights = readdlm("FMAlt_APMCOutputFilesJULIA/APMCsamples_FMAlt_JulRun#101_RetainedWeights.txt",'\t')

    #Get SurvivorParticleIdx
    SurvivorParticleIdx = readdlm("FMAlt_APMCOutputFilesJULIA/APMCsamples_FMAlt_JulRun#101_SurvivorParticleIdx.txt",'\t')

	  #Get summary statistic measure for each paramter set
    Particle_SummStat = readdlm("FMAlt_APMCOutputFilesJULIA/APMCsamples_FMAlt_JulRun#101_RetainedSummStat.txt",'\t')

    return ParticlesFromPrior, Particle_Weights, SurvivorParticleIdx, Particle_SummStat
end



#Six season fit, with susceptibility based on age
#Non-uniform priors in use for R_0, strain susceptibility modifier & ascertainment probabilities
#Three susceptibility age groups:  (i) 0-17; (ii) 18-64; (iii) 65+.
#Linear piecewise age-dependent ascertainment prob. per season. AND influenza type specific.
#Separate curves for each of the four influenza viruses
#Age 100 has max ascertainment. Five additional scaling factors, knots at ages: (i) 0; (ii) 2; (iii) 18; (iv) 65; (v) 85.
function  FMAlt_SampleFirstGenFn_SixSeasons_AgeSuscepTypeFour_AgeAndStrainAscertainLinearPiecewise(N::Int64)
#Inputs:
#   N - Number of samples from LHC parameter set to take forward as first generation

	#Distributions to be sampled from for each variable
  d = [Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125), #R_0
			Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
			Uniform(0,1),Uniform(0,1),Uniform(0,1), #Susceptibility per age band
	    Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs for A(H1N1)
		  Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1), #Ascertainment scaling for A(H1N1) (as proportion of value at age 100+)
		  Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs for A(H3N2)
		  Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1), #Ascertainment scaling for A(H3N2) (as proportion of value at age 100+)
		  Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs for B/Yam
		  Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1), #Ascertainment scaling for B/Yam (as proportion of value at age 100+)
		  Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs for B/Vic
		  Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1)] #Ascertainment scaling for B/Vic (as proportion of value at age 100+)

    ParamToFitNum = length(d) #Total number of parameters to fit

	#Sample and assign to array
    ParticlesFromPrior = zeros(N,ParamToFitNum)
    for ii = 1:ParamToFitNum
		  ParticlesFromPrior[:,ii] = rand(d[ii],N)
    end

  #---------------------------------------------------------------------------------------
  ### Enforce constraint on susceptibility, 18-64 category must be lower than the rest ###
  #---------------------------------------------------------------------------------------

  #Per particle, find minimum susceptibility in all classes except 18-64
  NonAdultSuscepIdx = [8,10]
  SuscepLimit = minimum(ParticlesFromPrior[:,NonAdultSuscepIdx],dims=2) #Particle per row, so take minimum in each row

  #For each particle, draw value from uniform distribution
  for ii = 1:N
    d_AdultSuscep = Uniform(0,SuscepLimit[ii])
    ParticlesFromPrior[ii,9] = rand(d_AdultSuscep,1)[1]
  end

  #Error check, are any adult suscep vals are greater than suscep values for other age groups
	SuscepValConstraintCheck = sum(ParticlesFromPrior[:,9].>ParticlesFromPrior[:,NonAdultSuscepIdx],dims=2)

	#If parameter values non-compatible, assign zero prior probability.
	if sum(SuscepValConstraintCheck) > 0
				error("Age group susceptibility value error, fails constraint check.")
	end

  #---------------------------------------------------------------------------------------

	#Assign weights to particle
    Particle_Weights = ones(N)

	return ParticlesFromPrior, Particle_Weights
end

#Six season fit, with susceptibility based on age
#Non-uniform priors in use for R_0, strain susceptibility modifier & ascertainment probabilities
#Three susceptibility age groups:  (i) 0-17; (ii) 18-64; (iii) 65+.
#Linear piecewise age-dependent ascertainment prob. per season. AND influenza type specific.
#Separate curves for Type A & Type B
#Age 100 has max ascertainment. Five additional scaling factors, knots at ages: (i) 0; (ii) 2; (iii) 18; (iv) 65; (v) 85.
function  FMAlt_SampleFirstGenFn_SixSeasons_AgeSuscepTypeFour_AgeAndTypeAscertainLinearPiecewise(N::Int64)
#Inputs:
#   N - Number of samples from LHC parameter set to take forward as first generation

	#Distributions to be sampled from for each variable
	d = [Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125), #R_0
		Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
		Uniform(0,1),Uniform(0,1),Uniform(0,1), #Susceptibility per age band
		Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs for type A
		Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1), #Ascertainment scaling for type A (as proportion of value at age 100+)
		Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs for type B
		Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1)] #Ascertainment scaling for type B (as proportion of value at age 100+)

    ParamToFitNum = length(d) #Total number of parameters to fit

	#Sample and assign to array
    ParticlesFromPrior = zeros(N,ParamToFitNum)
    for ii = 1:ParamToFitNum
		  ParticlesFromPrior[:,ii] = rand(d[ii],N)
    end

  #---------------------------------------------------------------------------------------
  ### Enforce constraint on susceptibility, 18-64 category must be lower than the rest ###
  #---------------------------------------------------------------------------------------

  #Per particle, find minimum susceptibility in all classes except 18-64
  NonAdultSuscepIdx = [8,10]
  SuscepLimit = minimum(ParticlesFromPrior[:,NonAdultSuscepIdx],dims=2) #Particle per row, so take minimum in each row

  #For each particle, draw value from uniform distribution
  for ii = 1:N
    d_AdultSuscep = Uniform(0,SuscepLimit[ii])
    ParticlesFromPrior[ii,9] = rand(d_AdultSuscep,1)[1]
  end

  #Error check, are any adult suscep vals are greater than suscep values for other age groups
	SuscepValConstraintCheck = sum(ParticlesFromPrior[:,9].>ParticlesFromPrior[:,NonAdultSuscepIdx],dims=2)

	#If parameter values non-compatible, assign zero prior probability.
	if sum(SuscepValConstraintCheck) > 0
				error("Age group susceptibility value error, fails constraint check.")
	end

  #---------------------------------------------------------------------------------------

	#Assign weights to particle
    Particle_Weights = ones(N)

	return ParticlesFromPrior, Particle_Weights
end

#Six season fit, with susceptibility based on age
#Non-uniform priors in use for R_0, strain susceptibility modifier & ascertainment probabilities
#Three susceptibility age groups:  (i) 0-17; (ii) 18-64; (iii) 65+.
#Linear piecewise age-dependent ascertainment prob. per season (& still season dependent!)
#Age 100 has max ascertainment. Five additional scaling factors, knots at ages: (i) 0; (ii) 2; (iii) 18; (iv) 65; (v) 85.
function  FMAlt_SampleFirstGenFn_SixSeasons_AgeSuscepTypeFour_AgeAscertainLinearPiecewise(N::Int64)
#Inputs:
#   N - Number of samples from LHC parameter set to take forward as first generation

	#Distributions to be sampled from for each variable
	d = [Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125), #R_0
		Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
		Uniform(0,1),Uniform(0,1),Uniform(0,1), #Susceptibility per age band
		Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs
		Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1)] #Ascertainment scaling (as proportion of value at age 100+)

    ParamToFitNum = length(d) #Total number of parameters to fit

	#Sample and assign to array
    ParticlesFromPrior = zeros(N,ParamToFitNum)
    for ii = 1:ParamToFitNum
		  ParticlesFromPrior[:,ii] = rand(d[ii],N)
    end

  #---------------------------------------------------------------------------------------
  ### Enforce constraint on susceptibility, 18-64 category must be lower than the rest ###
  #---------------------------------------------------------------------------------------

  #Per particle, find minimum susceptibility in all classes except 18-64
  NonAdultSuscepIdx = [8 10]
  SuscepLimit = minimum(ParticlesFromPrior[:,NonAdultSuscepIdx], dims=2) #Particle per row, so take minimum in each row

  #For each particle, draw value from uniform distribution
  for ii = 1:N
    d_AdultSuscep = Uniform(0,SuscepLimit[ii])
    ParticlesFromPrior[ii,9] = rand(d_AdultSuscep,1)[1]
  end
  #---------------------------------------------------------------------------------------

	#Assign weights to particle
    Particle_Weights = ones(N)

	return ParticlesFromPrior, Particle_Weights
end

#Six season fit, with susceptibility based on age
#Non-uniform priors in use for R_0, strain susceptibility modifier & ascertainment probabilities
#Age-dependent ascertainment prob. per season (& still season dependent!)
#Four age groups:  (i) 0-17; (ii) 18-64; (iii) 65-84; (iv) 85+.
function  FMAlt_SampleFirstGenFn_SixSeasons_AgeSuscepTypeFour_AgeAscertainTypeThree(N::Int64)
#Inputs:
#   N - Number of samples from LHC parameter set to take forward as first generation

	#Distributions to be sampled from for each variable
	d = [Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125), #R_0
		Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
		Uniform(0,1),Uniform(0,1),Uniform(0,1), #Susceptibility per age band
		Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs
		Gamma(4,1),Gamma(4,1)] #Ascertainment scale by age

    ParamToFitNum = length(d) #Total number of parameters to fit

	#Sample and assign to array
    ParticlesFromPrior = zeros(N,ParamToFitNum)
    for ii = 1:ParamToFitNum
		  ParticlesFromPrior[:,ii] = rand(d[ii],N)
    end

  #---------------------------------------------------------------------------------------
  ### Enforce constraint on susceptibility, 18-64 category must be lower than the rest ###
  #---------------------------------------------------------------------------------------

  #Per particle, find minimum susceptibility in all classes except 18-64
  NonAdultSuscepIdx = [8 10]
  SuscepLimit = minimum(ParticlesFromPrior[:,NonAdultSuscepIdx],2) #Particle per row, so take minimum in each row

  #For each particle, draw value from uniform distribution
  for ii = 1:N
    d_AdultSuscep = Uniform(0,SuscepLimit[ii])
    ParticlesFromPrior[ii,9] = rand(d_AdultSuscep,1)[1]
  end
  #---------------------------------------------------------------------------------------

	#Assign weights to particle
    Particle_Weights = ones(N)

	return ParticlesFromPrior, Particle_Weights
end



#-------------------------------------------------------------------------------
#  (iii) SPECIFY PRIOR DENSITY
#-------------------------------------------------------------------------------


#Six season fit, with susceptibility based on age only
#Non-uniform priors used for R_0, strain modifier and ascertainment probability
#Three susceptibility age groups: (i) 0-17; (ii) 18-64; (iii) 65+.
#Linear piecewise age-dependent ascertainment prob. per season. AND influenza type specific
#Separate curves for each of the four influenza viruses
#Age 100 has max ascertainment. Five additional scaling factors, knots at ages: (i) 0; (ii) 2; (iii) 18; (iv) 65; (v) 85.
function APMC_FMAltPrior_SixSeasons_AgeSuscepTypeFour_AgeAndStrainAscertainLinearPiecewise(x)

#Inputs:
#   x - proposed particle value to be tested

#Outputs:
#   PriorProb - Prior likelihood value of current particle (parameter set)
#   InBoundFlag - Denotes whether proposed particle values are within prior bounds

	#Get number of dimensions of x
	Particle_nDims = ndims(x)

	#--------------------------------------------------------------------------
	#COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
	#use)

	#Distributions to be sampled from for each variable
	d = [Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125), #R_0
			Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
			Uniform(0,1),Uniform(0,1),Uniform(0,1), #Susceptibility per age band
	    Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs for A(H1N1)
		  Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1), #Ascertainment scaling for A(H1N1) (as proportion of value at age 100+)
		  Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs for A(H3N2)
		  Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1), #Ascertainment scaling for A(H3N2) (as proportion of value at age 100+)
		  Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs for B/Yam
		  Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1), #Ascertainment scaling for B/Yam (as proportion of value at age 100+)
		  Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs for B/Vic
		  Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1)] #Ascertainment scaling for B/Vic (as proportion of value at age 100+)

	if Particle_nDims == 1
		ParticleNum = 1
		ParamNum = length(x)

		PriorProb = pdf.(d[1], x[1])

	    for ii = 2:ParamNum
	        PriorProb = PriorProb.*pdf.(d[ii], x[ii])
	    end

    #-----------------------------------------------------------------------------------
    ### Perform susceptibility parameter constraint check
    ### Enforce adult age groups have lower susceptibility (18-64 age group less than others)
    #-----------------------------------------------------------------------------------

			#Assign particle param values to distinct variable names
			AdultSuscepVal = x[9]

			#Specify values of other suscep param elements
			OtherSuscepValsIdx = [8,10] #Get relevant indices
			OtherSuscepVals = x[OtherSuscepValsIdx]

			#Check if any adult suscep vals are greater than suscep values for other age groups
			SuscepValConstraintCheck = sum(AdultSuscepVal.>OtherSuscepVals)

			#If parameter values non-compatible, assign zero prior probability.
			if SuscepValConstraintCheck > 0
				PriorProb = 0
			end
   #-----------------------------------------------------------------------------------

		#Alter flag value for implausible parameter set to 0
		if PriorProb==0
			InBoundFlag = 0::Int64
		else
			InBoundFlag = 1::Int64
		end
	else
		ParticleNum = size(x,1)
		ParamNum = size(x,2)

		PriorProb = pdf.(d[1], x[:,1])
	    for ii = 2:ParamNum
			  PriorProb = PriorProb.*pdf.(d[ii], x[:,ii]);
	    end

    #-----------------------------------------------------------------------------------
    ### Perform susceptibility parameter constraint check
    ### Enforce adult age groups have lower susceptibility (18-64 age group less than others)
    #-----------------------------------------------------------------------------------

		#Assign particle param values to distinct variable names
		AdultSuscepVals = x[:,9]

		#Specify values of other suscep param elements
		OtherSuscepValsIdx = [8,10] #Get relevant indices
		OtherSuscepVals = x[:,OtherSuscepValsIdx]

		#Check if any adult suscep vals are greater than suscep values for other age groups
		SuscepValConstraintCheck = sum(AdultSuscepVals.>OtherSuscepVals,dims=2)

		#If parameter values non-compatible, assign zero prior probability.
      	PriorProb[vec(SuscepValConstraintCheck).>0] .= 0

   #-----------------------------------------------------------------------------------

		#Alter flag value for implausible parameter set to 0
		InBoundFlag = ones(Int64,ParticleNum)
		InBoundFlag[PriorProb.==0] .= 0
	end

	return PriorProb,InBoundFlag
end

#Six season fit, with susceptibility based on age only
#Non-uniform priors used for R_0, strain modifier and ascertainment probability
#Three susceptibility age groups: (i) 0-17; (ii) 18-64; (iii) 65+.
#Linear piecewise age-dependent ascertainment prob. per season. AND influenza type specific
#Separate curves for Type A & Type B
#Age 100 has max ascertainment. Five additional scaling factors, knots at ages: (i) 0; (ii) 2; (iii) 18; (iv) 65; (v) 85.
function APMC_FMAltPrior_SixSeasons_AgeSuscepTypeFour_AgeAndTypeAscertainLinearPiecewise(x)

#Inputs:
#   x - proposed particle value to be tested

#Outputs:
#   PriorProb - Prior likelihood value of current particle (parameter set)
#   InBoundFlag - Denotes whether proposed particle values are within prior bounds

	#Get number of dimensions of x
	Particle_nDims = ndims(x)

	#--------------------------------------------------------------------------
	#COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
	#use)

	#Distributions to be sampled from for each variable
	d = [Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125), #R_0
			Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
			Uniform(0,1),Uniform(0,1),Uniform(0,1), #Susceptibility per age band
		  Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs for type A
			Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1), #Ascertainment scaling values at knot ages (type A)
			Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.005),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs for type B
  			Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1)] #Ascertainment scaling values at knot ages (type B)

	if Particle_nDims == 1
		ParticleNum = 1
		ParamNum = length(x)

		PriorProb = pdf.(d[1], x[1])

	    for ii = 2:ParamNum
	        PriorProb = PriorProb.*pdf.(d[ii], x[ii])
	    end

    #-----------------------------------------------------------------------------------
    ### Perform susceptibility parameter constraint check
    ### Enforce adult age groups have lower susceptibility (18-64 age group less than others)
    #-----------------------------------------------------------------------------------

			#Assign particle param values to distinct variable names
			AdultSuscepVal = x[9]

			#Specify values of other suscep param elements
			OtherSuscepValsIdx = [8,10] #Get relevant indices
			OtherSuscepVals = x[OtherSuscepValsIdx]

			#Check if any adult suscep vals are greater than suscep values for other age groups
			SuscepValConstraintCheck = sum(AdultSuscepVal.>OtherSuscepVals)

			#If parameter values non-compatible, assign zero prior probability.
			if SuscepValConstraintCheck > 0
				PriorProb = 0
			end
   #-----------------------------------------------------------------------------------

		#Alter flag value for implausible parameter set to 0
		if PriorProb==0
			InBoundFlag = 0::Int64
		else
			InBoundFlag = 1::Int64
		end
	else
		ParticleNum = size(x,1)
		ParamNum = size(x,2)

		PriorProb = pdf.(d[1], x[:,1])
	    for ii = 2:ParamNum
			  PriorProb = PriorProb.*pdf.(d[ii], x[:,ii]);
	    end

    #-----------------------------------------------------------------------------------
    ### Perform susceptibility parameter constraint check
    ### Enforce adult age groups have lower susceptibility (18-64 age group less than others)
    #-----------------------------------------------------------------------------------

		#Assign particle param values to distinct variable names
		AdultSuscepVals = x[:,9]

		#Specify values of other suscep param elements
		OtherSuscepValsIdx = [8,10] #Get relevant indices
		OtherSuscepVals = x[:,OtherSuscepValsIdx]

		#Check if any adult suscep vals are greater than suscep values for other age groups
		SuscepValConstraintCheck = sum(AdultSuscepVals.>OtherSuscepVals,dims=2)

		#If parameter values non-compatible, assign zero prior probability.
      	PriorProb[vec(SuscepValConstraintCheck).>0] .= 0

   #-----------------------------------------------------------------------------------

		#Alter flag value for implausible parameter set to 0
		InBoundFlag = ones(Int64,ParticleNum)
		InBoundFlag[PriorProb.==0] .= 0
	end

	return PriorProb,InBoundFlag
end

#Six season fit, with susceptibility based on age only
#Non-uniform priors used for R_0, strain modifier and ascertainment probability
#Three susceptibility age groups: (i) 0-17; (ii) 18-64; (iii) 65+.
#Linear piecewise age-dependent ascertainment prob. per season (& still season dependent!)
#Age 100 has max ascertainment. Five additional scaling factors, knots at ages: (i) 0; (ii) 2; (iii) 18; (iv) 65; (v) 85.
function APMC_FMAltPrior_SixSeasons_AgeSuscepTypeFour_AgeAscertainLinearPiecewise(x)

#Inputs:
#   x - proposed particle value to be tested

#Outputs:
#   PriorProb - Prior likelihood value of current particle (parameter set)
#   InBoundFlag - Denotes whether proposed particle values are within prior bounds

	#Get number of dimensions of x
	Particle_nDims = ndims(x)

	#--------------------------------------------------------------------------
	#COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
	#use)

	#Distributions to be sampled from for each variable
	d = [Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125), #R_0
			Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
			Uniform(0,1),Uniform(0,1),Uniform(0,1), #Susceptibility per age band
		  Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs
			Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1),Uniform(0,1)]

	if Particle_nDims == 1
		ParticleNum = 1
		ParamNum = length(x)

		PriorProb = pdf.(d[1], x[1])

	    for ii = 2:ParamNum
	        PriorProb = PriorProb.*pdf.(d[ii], x[ii])
	    end

    #-----------------------------------------------------------------------------------
    ### Perform susceptibility parameter constraint check
    ### Enforce adult age groups have lower susceptibility (18-64 age group less than others)
    #-----------------------------------------------------------------------------------

			#Assign particle param values to distinct variable names
			AdultSuscepVal = x[9]

			#Specify values of other suscep param elements
			OtherSuscepValsIdx = [8,10] #Get relevant indices
			OtherSuscepVals = x[OtherSuscepValsIdx]

			#Check if any adult suscep vals are greater than suscep values for other age groups
			SuscepValConstraintCheck = sum(AdultSuscepVal.>OtherSuscepVals)

			#If parameter values non-compatible, assign zero prior probability.
			if SuscepValConstraintCheck > 0
				PriorProb = 0
			end
   #-----------------------------------------------------------------------------------

		#Alter flag value for implausible parameter set to 0
		if PriorProb==0
			InBoundFlag = 0::Int64
		else
			InBoundFlag = 1::Int64
		end
	else
		ParticleNum = size(x,1)
		ParamNum = size(x,2)

		PriorProb = pdf.(d[1], x[:,1])
	    for ii = 2:ParamNum
			PriorProb = PriorProb.*pdf.(d[ii], x[:,ii]);
	    end

    #-----------------------------------------------------------------------------------
    ### Perform susceptibility parameter constraint check
    ### Enforce adult age groups have lower susceptibility (18-64 age group less than others)
    #-----------------------------------------------------------------------------------

			#Assign particle param values to distinct variable names
			AdultSuscepVals = x[:,9]

			#Specify values of other suscep param elements
			OtherSuscepValsIdx = [8,10] #Get relevant indices
			OtherSuscepVals = x[:,OtherSuscepValsIdx]

			#Check if any adult suscep vals are greater than suscep values for other age groups
			SuscepValConstraintCheck = sum(AdultSuscepVals.>OtherSuscepVals,dims=2)

			#If parameter values non-compatible, assign zero prior probability.
			SuscepValConstraintIdx = convert(BitArray,SuscepValConstraintCheck)

      		#PriorProb[vec(SuscepValConstraintIdx)] = 0
			PriorProb[vec(SuscepValConstraintCheck).>0] .= 0

   #-----------------------------------------------------------------------------------

		#Alter flag value for implausible parameter set to 0
		InBoundFlag = ones(Int64,ParticleNum)
		InBoundFlag[PriorProb.==0] .= 0
	end

	return PriorProb,InBoundFlag
end

#Six season fit, with susceptibility based on age only
#Non-uniform priors used for R_0, strain modifier and ascertainment probability
#Age-dependent ascertainment probability (& still season dependent)
#Three age groups: (i) 0-17; (ii) 18-64; (iii) 65+.
function APMC_FMAltPrior_SixSeasons_AgeSuscepTypeFour_AgeAscertainTypeThree(x)

#Inputs:
#   x - proposed particle value to be tested

#Outputs:
#   PriorProb - Prior likelihood value of current particle (parameter set)
#   InBoundFlag - Denotes whether proposed particle values are within prior bounds

	#Get number of dimensions of x
	Particle_nDims = ndims(x)

	#--------------------------------------------------------------------------
	#COMPUTE PRIOR LIKELIHOOD (syntax dependent upon number of partciles in
	#use)

	#Distributions to be sampled from for each variable
	d = [Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.125),Gamma(10,0.15), #R_0
			Uniform(0,1),Uniform(0,1),Uniform(0,1), #Immunity propagation params
			Uniform(0,1),Uniform(0,1),Uniform(0,1), #Susceptibility per age band
		  Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02),Gamma(2,0.02), #Ascertainment probs
			Gamma(4,1),Gamma(4,1)]

	if Particle_nDims == 1
		ParticleNum = 1
		ParamNum = length(x)

		PriorProb = pdf.(d[1], x[1])

	    for ii = 2:ParamNum
	        PriorProb = PriorProb.*pdf.(d[ii], x[ii])
	    end

    #-----------------------------------------------------------------------------------
    ### Perform susceptibility parameter constraint check
    ### Enforce adult age groups have lower susceptibility (18-64 age group less than others)
    #-----------------------------------------------------------------------------------

			#Assign particle param values to distinct variable names
			AdultSuscepVal = x[9]

			#Specify values of other suscep param elements
			OtherSuscepValsIdx = [8,10] #Get relevant indices
			OtherSuscepVals = x[OtherSuscepValsIdx]

			#Check if any adult suscep vals are greater than suscep values for other age groups
			SuscepValConstraintCheck = sum(AdultSuscepVal.>OtherSuscepVals)

			#If parameter values non-compatible, assign zero prior probability.
			if SuscepValConstraintCheck > 0
				PriorProb = 0
			end
   #-----------------------------------------------------------------------------------

   #-----------------------------------------------------------------------------------
    ### Perform ascertainment parameter constraint check
    ### Ensure scaled age probability is below one in all age classes
    #-----------------------------------------------------------------------------------

			#Assign particle param values to distinct variable names
			AscertainSeasonVals = x[end-7:end-2]
      AscertainAgeScaleVals = x[end-1:end]

		  #Take product of age scaling values with base values
      ProductAscertainVals = AscertainSeasonVals.*AscertainAgeScaleVals'

		  #Check if any exceed 1. If so, constraint check is failed
		  AscertainValConstraintCheck = sum(ProductAscertainVals.>1)

			#If parameter values non-compatible, assign zero prior probability.
			if AscertainValConstraintCheck > 0
				PriorProb = 0
			end
   #-----------------------------------------------------------------------------------

		#Alter flag value for implausible parameter set to 0
		if PriorProb==0
			InBoundFlag = 0::Int64
		else
			InBoundFlag = 1::Int64
		end
	else
		ParticleNum = size(x,1)
		ParamNum = size(x,2)

		PriorProb = pdf.(d[1], x[:,1])
	    for ii = 2:ParamNum
			PriorProb = PriorProb.*pdf.(d[ii], x[:,ii]);
	    end

    #-----------------------------------------------------------------------------------
    ### Perform susceptibility parameter constraint check
    ### Enforce adult age groups have lower susceptibility (18-64 age group less than others)
    #-----------------------------------------------------------------------------------

			#Assign particle param values to distinct variable names
			AdultSuscepVals = x[:,9]

			#Specify values of other suscep param elements
			OtherSuscepValsIdx = [8,10] #Get relevant indices
			OtherSuscepVals = x[:,OtherSuscepValsIdx]

			#Check if any adult suscep vals are greater than suscep values for other age groups
			SuscepValConstraintCheck = sum(AdultSuscepVals.>OtherSuscepVals,2)

			#If parameter values non-compatible, assign zero prior probability.
      		PriorProb[vec(SuscepValConstraintCheck).>0] = 0
   #-----------------------------------------------------------------------------------

    #-----------------------------------------------------------------------------------
    ### Perform ascertainment parameter constraint check
    ### Ensure scaled age probability is below one in all age classes
    #-----------------------------------------------------------------------------------

			#Assign particle param values to distinct variable names
			AscertainSeasonVals = x[:,end-7:end-2]
      		AscertainAgeScaleVals = x[:,end-1:end]

		  #Take product of age scaling values with base values
        AscertainValConstraintCheck = zeros(ParticleNum)
        for ii = 1:ParticleNum
          ProductAscertainVals = AscertainSeasonVals[ii,:].*AscertainAgeScaleVals[ii,:]'

          #Check if any exceed 1. If so, constraint check is failed
          AscertainValConstraintCheck[ii] = sum(ProductAscertainVals.>1)
        end

			#If parameter values non-compatible, assign zero prior probability.
      PriorProb[vec(AscertainValConstraintCheck.>0)] = 0

   #-----------------------------------------------------------------------------------

		#Alter flag value for implausible parameter set to 0
		InBoundFlag = ones(Int64,ParticleNum)
		InBoundFlag[PriorProb.==0] .= 0
	end

	return PriorProb,InBoundFlag

end

#--------------------------------------------------------------------------
# (iv) SUMMARY STATISTIC CALCULATION FUNCTIONS
#--------------------------------------------------------------------------

function  FMAlt_SummStatFun_PoissDev(ObservedData,SimnData)

    #Disaggregate simulation data
    x = SimnData[1]
    Total_I = SimnData[2]

	#Compute poisson deviance
	SummStatIterNum = length(ObservedData)
	TempSum = 0.0
    for ii = 1:SummStatIterNum
        if ObservedData[ii].!=0
            TempSum = TempSum + (ObservedData[ii]*log(ObservedData[ii]/x[ii])) - (ObservedData[ii] - x[ii])
        end
    end
    SummStatVal_Overall = 2*(TempSum + sum(x[ObservedData.==0]))

    #Perform temporal check
    NumOfSeasons = convert(Int64,size(Total_I,1)/366) #%Number of seasons obtained by dividing number of daily records by days in yr (+1 to account for day 0 recording1)
	AgeGrpNum = convert(Int64,size(Total_I,2)) #Number of age groups, number of columns of Total_I
	Total_I_Array = zeros(NumOfSeasons,AgeGrpNum,366)
	Total_I_Transpose = Total_I'
    for jj = 1:NumOfSeasons
       StartIdx = ((jj-1)*366) + 1
       EndIdx = jj*366

	   #Allocate seasons worth of records to new array, slice per season
       #Take transpose of Total_I so array dimensions agree
       Total_I_Array[jj,:,:] = Total_I_Transpose[:,StartIdx:EndIdx]
    end

	#Find day of seasonal year in which infection is at peak
	MaxInfValBySeason,inds = findmax(Total_I_Array[4:end,:,:],dims=3)
	MaxInfIdxBySeason = map(x->x[3], inds)


    #If peak outside Sep-Feb, set TemporalFlag to 0 and put amended
    #SummStatVal as infinity
    if sum(MaxInfIdxBySeason.>182) == 0 #181 days September-February. Plus account for initial value
        AmendedSummStatVal = SummStatVal_Overall
    else
        AmendedSummStatVal = Inf
    end

    return AmendedSummStatVal

end

#-------------------------------------------------------------------------------
# (v) FUNCTIONS TO APPLY ASCERTAINMENT PROBABILITY TO INFECTION CASE LOAD
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


# Season-dependent, Age agnostic
function  SeasonDepOnlyAscertainmentFn(AscertainProb,SeasonRatePerAgeStrain,
										AgeBandedCaseRateTotalSum_FitCheckSeasons,
										FitAggAgeFlag,AgeBandNum,
										AgeBandBounds,M)

	#Error check on parameter vector length
	RetainedSeasonNum = size(SeasonRatePerAgeStrain,1)
	if length(AscertainProb) != RetainedSeasonNum
		error("Number of parameters allocated to AscertainProb is $(length(AscertainProb)). Should be $RetainedSeasonNum.")
	end

	FluToGPconversion = 100000*AscertainProb #Scale to rate per 100,000 popualation
	for jj = 1:length(AscertainProb)
		SeasonRatePerAgeStrain[jj,:,:] = AgeBandedCaseRateTotalSum_FitCheckSeasons[jj,:,:].*FluToGPconversion[jj]
	end

	return SeasonRatePerAgeStrain
end

# Age-dependent, season agnostic
function  AgeDepOnlyAscertainmentFn(AscertainProb,SeasonRatePerAgeStrain,
										AgeBandedCaseRateTotalSum_FitCheckSeasons,
										FitAggAgeFlag,AgeBandNum,
										AgeBandBounds,M)

	#Error check on parameter vector length
	if length(AscertainProb) != AgeBandNum
		error("Number of parameters allocated to AscertainProb is $(length(AscertainProb)). Should be $AgeBandNum.")
	end

	if FitAggAgeFlag == 0  #Use single year age classes

		#Disaggregate AgeBandBounds
		AgeBandLowerBounds = AgeBandBounds[1,:]::Array{Int64,1}
		AgeBandUpperBounds = AgeBandBounds[2,:]::Array{Int64,1}

		#Initialise array - Ascertain prob. by yr of age
		AscertainProbByAge = zeros(M)

		#Populate AscertainProbByAge, iterate through each age band
		#Assign AscertainProb to relevant single year of age classess
		for ii = 1:AgeBandNum
			AscertainAgeStartIdx = AgeBandLowerBounds[ii] + 1 #Ages begin at 0. Add 1 to align with array indexing
			AscertainAgeEndIdx = AgeBandUpperBounds[ii] + 1
			if ii == AgeBandNum #If eldest age category, sum to final column.
				AscertainProbByAge[AscertainAgeStartIdx:end] .= AscertainProb[ii]
			else
				AscertainProbByAge[AscertainAgeStartIdx:AscertainAgeEndIdx] .= AscertainProb[ii]
			end
		end

		for jj = 1:length(AscertainProbByAge)
			SeasonRatePerAgeStrain[:,jj,:] = AgeBandedCaseRateTotalSum_FitCheckSeasons[:,jj,:].*AscertainProbByAge[jj]*100000
		end
	elseif FitAggAgeFlag == 1 #Use intervals/age buckets designated by AgeBandBounds
		FluToGPconversion = 100000*AscertainProb #Scale to rate per 100,000 popualation
		for jj = 1:length(AscertainProb)
			SeasonRatePerAgeStrain[:,jj,:] = AgeBandedCaseRateTotalSum_FitCheckSeasons[:,jj,:].*FluToGPconversion[jj]
		end
	end

	return SeasonRatePerAgeStrain
end

# Seasons specific baseline, with age group scaling. Step function.
function  SeasonWithAgeScaleStepAscertainmentFn(AscertainProb,SeasonRatePerAgeStrain,
										AgeBandedCaseRateTotalSum_FitCheckSeasons,
										FitAggAgeFlag,AgeBandNum,
										AgeBandBounds,M)


	#Error check on parameter vector length
	RetainedSeasonNum = size(SeasonRatePerAgeStrain,1)
	if length(AscertainProb) != (RetainedSeasonNum+AgeBandNum-1)
		error("Number of parameters allocated to AscertainProb is $(length(AscertainProb)). Should be $(RetainedSeasonNum+AgeBandNum-1).")
	end

	#First, scale values by seasonal ascertainment probability
	AscertainProbBySeason = AscertainProb[1:RetainedSeasonNum]
	FluToGPconversion = 100000*AscertainProbBySeason #Scale to rate per 100,000 popualation
	for jj = 1:RetainedSeasonNum
		SeasonRatePerAgeStrain[jj,:,:] = AgeBandedCaseRateTotalSum_FitCheckSeasons[jj,:,:].*FluToGPconversion[jj]
	end

	#Second, scale age groups with age-specific ascertainment modification
	AscertainProbAgeScale = AscertainProb[RetainedSeasonNum+1:end]

	if FitAggAgeFlag == 0  #Use single year age classes

		#Disaggregate AgeBandBounds
		AgeBandLowerBounds = AgeBandBounds[1,:]::Array{Int64,1}
		AgeBandUpperBounds = AgeBandBounds[2,:]::Array{Int64,1}

		#Initialise array - Ascertain prob. by yr of age
		AscertainProbScaleByAge = zeros(M)

		#Populate AscertainProbByAge, iterate through each age band
		#Assign AscertainProb to relevant single year of age classess
		for ii = 1:AgeBandNum
			AscertainAgeStartIdx = AgeBandLowerBounds[ii] + 1 #Ages begin at 0. Add 1 to align with array indexing
			AscertainAgeEndIdx = AgeBandUpperBounds[ii] + 1
			if ii == AgeBandNum #If eldest age category, sum to final column.
				AscertainProbScaleByAge[AscertainAgeStartIdx:end] .= AscertainProbAgeScale[ii-1]
			elseif ii == 1 #Youngest age category, unmodified ascertainment prob.
				AscertainProbScaleByAge[AscertainAgeStartIdx:AscertainAgeEndIdx] .= 1
			else #All other age categories, assign across Idx interval age scale modifier
				AscertainProbScaleByAge[AscertainAgeStartIdx:AscertainAgeEndIdx] .= AscertainProbAgeScale[ii-1]
			end
		end

		for jj = 1:length(AscertainProbScaleByAge)
			SeasonRatePerAgeStrain[:,jj,:] = SeasonRatePerAgeStrain[:,jj,:].*AscertainProbScaleByAge[jj]
		end
	elseif FitAggAgeFlag == 1 #Use intervals/age buckets designated by AgeBandBounds

		#Leave youngest age category unscaled (in effect has strain ascertainment scale value of 1)
		#i.e. First column unaltered
		#SeasonRatePerAgeStrain[:,2:end,:] = SeasonRatePerAgeStrain[:,2:end,:].*AscertainProbAgeScale
		for jj = 1:length(AscertainProbAgeScale)
			SeasonRatePerAgeStrain[:,jj+1,:] = SeasonRatePerAgeStrain[:,jj+1,:].*AscertainProbAgeScale[jj]
		end
	end

	return SeasonRatePerAgeStrain
end

# Seasons specific baseline, with age group scaling. Piecewise linear function.
function  SeasonWithAgeScaleLinearAscertainmentFn(AscertainProb,SeasonRatePerAgeStrain,
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
		FluToGPconversionBySeason = 100000*MaxAscertainProbBySeason
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


# Seasons specific baseline, with age group scaling. Piecewise linear function.
function  SeasonWithAgeScaleLinearAscertainmentFn_TypeSpecific(AscertainProb,SeasonRatePerAgeStrain,
										AgeBandedCaseRateTotalSum_FitCheckSeasons,
										FitAggAgeFlag,AgeBandNum,
										AgeBandBounds,M)

	#Error check on parameter vector length
	RetainedSeasonNum = size(SeasonRatePerAgeStrain,1)
	if length(AscertainProb) != (2*(RetainedSeasonNum+AgeBandNum))
		error("Number of parameters allocated to AscertainProb is $(length(AscertainProb)). Should be $(2*(RetainedSeasonNum+AgeBandNum)).")
	end

	#Note, only invoked if fitting to single year age classes
	#Fitting to binned data, throws error message & ends programme
	if FitAggAgeFlag == 0

		#######################################################
		### Type A & B computations performed concurrently! ###
		#######################################################

		### First, scale values by max seasonal ascertainment probability (attained by those aged 100yrs) ###
		MaxAscertainProbBySeason_TypeA = AscertainProb[1:RetainedSeasonNum]
		MaxAscertainProbBySeason_TypeB = AscertainProb[(RetainedSeasonNum+AgeBandNum)+1:(RetainedSeasonNum+AgeBandNum)+RetainedSeasonNum]

		FluToGPconversionBySeason_TypeA = 100000*MaxAscertainProbBySeason_TypeA
		FluToGPconversionBySeason_TypeB = 100000*MaxAscertainProbBySeason_TypeB
		for jj = 1:RetainedSeasonNum
			SeasonRatePerAgeStrain[jj,:,1:2] = AgeBandedCaseRateTotalSum_FitCheckSeasons[jj,:,1:2].*FluToGPconversionBySeason_TypeA[jj] #Type A scaling
			SeasonRatePerAgeStrain[jj,:,3:4] = AgeBandedCaseRateTotalSum_FitCheckSeasons[jj,:,3:4].*FluToGPconversionBySeason_TypeB[jj] #Type B scaling

		end

		### Second, scale age groups with age-specific ascertainment modification ###

		#Assign ascertainment age band end value scaling factors to variable
		AscertainProbAgeScale_TypeA = AscertainProb[RetainedSeasonNum+1:RetainedSeasonNum+AgeBandNum]
		AscertainProbAgeScale_TypeB = AscertainProb[(RetainedSeasonNum+AgeBandNum)+RetainedSeasonNum+1:end]

		#End endpoint value of 1 (max age has highest ascertainment, unscaled)
		AscertainProbAgeScaleIncEndPoint_TypeA = [AscertainProbAgeScale_TypeA;1]
		AscertainProbAgeScaleIncEndPoint_TypeB = [AscertainProbAgeScale_TypeB;1]

		#Concatenate influenza ascertainment prob by type into a single array
		AscertainProbAgeScaleIncEndPoint = hcat(AscertainProbAgeScaleIncEndPoint_TypeA,AscertainProbAgeScaleIncEndPoint_TypeB)

		#Disaggregate AgeBandBounds
		AgeBandLowerBounds = AgeBandBounds[1,:]::Array{Int64,1}
		AgeBandUpperBounds = AgeBandBounds[2,:]::Array{Int64,1}

		#Initialise array - Ascertain prob. by yr of age and influenza type
		InfluenzaTypeNum = 2
		AscertainProbScaleByAge = zeros(M,InfluenzaTypeNum)

		#Iterate over each year of age ii
		for ii = 1:M

			#Convert loop index to year of age
			YrOfAgeVal = ii - 1

		  	#Check upper bound, first <= than.
		  	AgeBandVectorIdx = findfirst(YrOfAgeVal .<= AgeBandUpperBounds)

			for InfluenzaTypeIdx = 1:InfluenzaTypeNum
			  	#Identify values at edge of bins, Entry idx and idx+1.
			  	LeftEdgeVal = AscertainProbAgeScaleIncEndPoint[AgeBandVectorIdx,InfluenzaTypeIdx]
			  	RightEdgeVal = AscertainProbAgeScaleIncEndPoint[AgeBandVectorIdx+1,InfluenzaTypeIdx]

			  	#For age ii, ascertain scale value is
			  	#AscertainScale[idx] + (IdxInBin/NumYrsInBin)*(AscertainScale[idx+1] - AscertainScale[idx])
				if AgeBandVectorIdx == length(AgeBandLowerBounds) #FOr last age category, get age years up to and including 100
					AgeBandSpan = 100 - AgeBandLowerBounds[AgeBandVectorIdx]
				else
				  AgeBandSpan = AgeBandLowerBounds[AgeBandVectorIdx+1] - AgeBandLowerBounds[AgeBandVectorIdx]
				end

			  	IdxInAgeBand = (ii-1) - AgeBandLowerBounds[AgeBandVectorIdx] #Subtract 1 as ages begin from zero!
			  	AscertainProbScaleByAge[ii,InfluenzaTypeIdx] = LeftEdgeVal + ((IdxInAgeBand/AgeBandSpan)*(RightEdgeVal - LeftEdgeVal))

				#Error check
				if AgeBandSpan == 0
					error("An age band contains only a single year of age. Incompatible with the selected piecewise linear ascertainment function.")
				end
			end

		end

		for jj = 1:M
			SeasonRatePerAgeStrain[:,jj,1:2] = SeasonRatePerAgeStrain[:,jj,1:2].*AscertainProbScaleByAge[jj,1] #Apply type A age ascertainment profile
			SeasonRatePerAgeStrain[:,jj,3:4] = SeasonRatePerAgeStrain[:,jj,3:4].*AscertainProbScaleByAge[jj,2] #Apply type B age ascertainment profile
		end

	else
		error("Banded age data (FitAggAgeFlag set to 1) not compatible with choice of ascertainment function.")
	end

	return SeasonRatePerAgeStrain
end


# Seasons specific baseline, with age group scaling. Piecewise linear function.
function  SeasonWithAgeScaleLinearAscertainmentFn_StrainSpecific(AscertainProb,SeasonRatePerAgeStrain,
										AgeBandedCaseRateTotalSum_FitCheckSeasons,
										FitAggAgeFlag,AgeBandNum,
										AgeBandBounds,M)

	#Error check on parameter vector length
	RetainedSeasonNum = size(SeasonRatePerAgeStrain,1)
	if length(AscertainProb) != (4*(RetainedSeasonNum+AgeBandNum))
		error("Number of parameters allocated to AscertainProb is $(length(AscertainProb)). Should be $(4*(RetainedSeasonNum+AgeBandNum)).")
	end

	#Note, only invoked if fitting to single year age classes
	#Fitting to binned data, throws error message & ends programme
	if FitAggAgeFlag == 0

		############################################################
		### Strain-specific computations performed concurrently! ###
		############################################################

		### First, scale values by max seasonal ascertainment probability (attained by those aged 100yrs) ###
		MaxAscertainProbBySeason_H1 = AscertainProb[1:RetainedSeasonNum]
		MaxAscertainProbBySeason_H3 = AscertainProb[(RetainedSeasonNum+AgeBandNum)+1:(RetainedSeasonNum+AgeBandNum)+RetainedSeasonNum]
		MaxAscertainProbBySeason_BYam = AscertainProb[(2*(RetainedSeasonNum+AgeBandNum))+1:(2*(RetainedSeasonNum+AgeBandNum))+RetainedSeasonNum]
   	MaxAscertainProbBySeason_BVic = AscertainProb[(3*(RetainedSeasonNum+AgeBandNum))+1:(3*(RetainedSeasonNum+AgeBandNum))+RetainedSeasonNum]

		FluToGPconversionBySeason_H1 = 100000*MaxAscertainProbBySeason_H1
		FluToGPconversionBySeason_H3 = 100000*MaxAscertainProbBySeason_H3
   	FluToGPconversionBySeason_BYam = 100000*MaxAscertainProbBySeason_BYam
    FluToGPconversionBySeason_BVic = 100000*MaxAscertainProbBySeason_BVic
		for jj = 1:RetainedSeasonNum
			SeasonRatePerAgeStrain[jj,:,1] = AgeBandedCaseRateTotalSum_FitCheckSeasons[jj,:,1].*FluToGPconversionBySeason_H1[jj] #A(H1N1) scaling
      SeasonRatePerAgeStrain[jj,:,2] = AgeBandedCaseRateTotalSum_FitCheckSeasons[jj,:,2].*FluToGPconversionBySeason_H3[jj] #A(H3N2) scaling
			SeasonRatePerAgeStrain[jj,:,3] = AgeBandedCaseRateTotalSum_FitCheckSeasons[jj,:,3].*FluToGPconversionBySeason_BYam[jj] #B/Yam scaling
			SeasonRatePerAgeStrain[jj,:,4] = AgeBandedCaseRateTotalSum_FitCheckSeasons[jj,:,4].*FluToGPconversionBySeason_BVic[jj] #B/Vicscaling
		end

		### Second, scale age groups with age-specific ascertainment modification ###

		#Assign ascertainment age band end value scaling factors to variable
		AscertainProbAgeScale_H1 = AscertainProb[RetainedSeasonNum+1:RetainedSeasonNum+AgeBandNum]
		AscertainProbAgeScale_H3 = AscertainProb[(RetainedSeasonNum+AgeBandNum)+RetainedSeasonNum+1:2*(RetainedSeasonNum+AgeBandNum)]
    AscertainProbAgeScale_BYam = AscertainProb[(2*(RetainedSeasonNum+AgeBandNum))+RetainedSeasonNum+1:3*(RetainedSeasonNum+AgeBandNum)]
    AscertainProbAgeScale_BVic = AscertainProb[(3*(RetainedSeasonNum+AgeBandNum))+RetainedSeasonNum+1:4*(RetainedSeasonNum+AgeBandNum)]

		#End endpoint value of 1 (max age has highest ascertainment, unscaled)
		AscertainProbAgeScaleIncEndPoint_H1 = [AscertainProbAgeScale_H1;1]
		AscertainProbAgeScaleIncEndPoint_H3 = [AscertainProbAgeScale_H3;1]
    AscertainProbAgeScaleIncEndPoint_BYam = [AscertainProbAgeScale_BYam;1]
    AscertainProbAgeScaleIncEndPoint_BVic = [AscertainProbAgeScale_BVic;1]

		#Concatenate influenza ascertainment prob by type into a single array
		AscertainProbAgeScaleIncEndPoint = hcat(AscertainProbAgeScaleIncEndPoint_H1,AscertainProbAgeScaleIncEndPoint_H3,AscertainProbAgeScaleIncEndPoint_BYam,AscertainProbAgeScaleIncEndPoint_BVic)

		#Disaggregate AgeBandBounds
		AgeBandLowerBounds = AgeBandBounds[1,:]::Array{Int64,1}
		AgeBandUpperBounds = AgeBandBounds[2,:]::Array{Int64,1}

		#Initialise array - Ascertain prob. by yr of age and influenza type
		InfluenzaTypeNum = 4
		AscertainProbScaleByAge= zeros(M,InfluenzaTypeNum)

		#Iterate over each year of age ii
		for ii = 1:M

			#Convert loop index to year of age
			YrOfAgeVal = ii - 1

		  	#Check upper bound, first <= than.
		  	AgeBandVectorIdx = findfirst(YrOfAgeVal .<= AgeBandUpperBounds)

			for InfluenzaTypeIdx = 1:InfluenzaTypeNum
			  	#Identify values at edge of bins, Entry idx and idx+1.
			  	LeftEdgeVal = AscertainProbAgeScaleIncEndPoint[AgeBandVectorIdx,InfluenzaTypeIdx]
			  	RightEdgeVal = AscertainProbAgeScaleIncEndPoint[AgeBandVectorIdx+1,InfluenzaTypeIdx]

			  	#For age ii, ascertain scale value is
			  	#AscertainScale[idx] + (IdxInBin/NumYrsInBin)*(AscertainScale[idx+1] - AscertainScale[idx])
				if AgeBandVectorIdx == length(AgeBandLowerBounds) #FOr last age category, get age years up to and including 100
					AgeBandSpan = 100 - AgeBandLowerBounds[AgeBandVectorIdx]
				else
				  AgeBandSpan = AgeBandLowerBounds[AgeBandVectorIdx+1] - AgeBandLowerBounds[AgeBandVectorIdx]
				end

			  	IdxInAgeBand = (ii-1) - AgeBandLowerBounds[AgeBandVectorIdx] #Subtract 1 as ages begin from zero!
			  	AscertainProbScaleByAge[ii,InfluenzaTypeIdx] = LeftEdgeVal + ((IdxInAgeBand/AgeBandSpan)*(RightEdgeVal - LeftEdgeVal))

				#Error check
				if AgeBandSpan == 0
					error("An age band contains only a single year of age. Incompatible with the selected piecewise linear ascertainment function.")
				end
			end
		end

    for InfluenzaTypeIdx = 1:InfluenzaTypeNum
		  for jj = 1:M
			  SeasonRatePerAgeStrain[:,jj,InfluenzaTypeIdx] = SeasonRatePerAgeStrain[:,jj,InfluenzaTypeIdx].*AscertainProbScaleByAge[jj,InfluenzaTypeIdx] #Apply age ascertainment profile for given strain
		  end
    end

	else
		error("Banded age data (FitAggAgeFlag set to 1) not compatible with choice of ascertainment function.")
	end

	return SeasonRatePerAgeStrain
end

#-------------------------------------------------------------------------------
# (vi) FUNCTIONS TO CONSTRUCT SUSCEPTIBILITY ARRAY
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

#Age & strain specific suscep.
function AgeAndStrainDepSuscepFn(x,AgeGrpSuscepParamNum,TransmissExpHistParamNum,NumOfStrains)

	AgeSuscepParamNum = (AgeGrpSuscepParamNum*NumOfStrains)::Int64
	AgeSuscepParams = x[TransmissExpHistParamNum+1:TransmissExpHistParamNum+AgeSuscepParamNum]::Array{Float64,1}

	#Put AgeSuscepParams vector into array form
	#Row per susceptibility grouping, column per strain
	AgeSuscepByGroup = reshape(AgeSuscepParams,(AgeGrpSuscepParamNum,NumOfStrains))

	return AgeSuscepByGroup, AgeSuscepParamNum
end

#Age grp specific, with strain modifier
function AgeDepWithStrainScalingSuscepFn(x,AgeGrpSuscepParamNum,TransmissExpHistParamNum,NumOfStrains)
	AgeSuscepParamNum = (AgeGrpSuscepParamNum + (NumOfStrains - 1))::Int64
	AgeSuscepParamsOnly = x[TransmissExpHistParamNum+1:TransmissExpHistParamNum+AgeGrpSuscepParamNum]::Array{Float64,1}
	StrainModifierSuscepParamsOnly = x[TransmissExpHistParamNum+AgeGrpSuscepParamNum+1:TransmissExpHistParamNum+AgeSuscepParamNum]::Array{Float64,1}

	#Construct AgeSuscepByGroup from AgeSuscepParamsOnly &
	# StrainModifierSuscepParamsOnly vector
	#Row per susceptibility grouping, column per strain
	AgeSuscepByGroup_StrainUnmod = repeat(AgeSuscepParamsOnly,outer =(1,NumOfStrains)) #Populate each column with age-specific suscep.
	AgeSuscepByGroup = [AgeSuscepByGroup_StrainUnmod[:,1] AgeSuscepByGroup_StrainUnmod[:,2:end].*StrainModifierSuscepParamsOnly'];

	return AgeSuscepByGroup, AgeSuscepParamNum
end

#Strain specific, with age group modifier
function StrainDepWithAgeScalingSuscepFn(x,AgeGrpSuscepParamNum,TransmissExpHistParamNum,NumOfStrains)
	AgeSuscepParamNum = (NumOfStrains + (AgeGrpSuscepParamNum - 1))::Int64
	StrainSuscepParamsOnly = x[TransmissExpHistParamNum+1:TransmissExpHistParamNum+NumOfStrains]::Array{Float64,1}
	AgeModifierSuscepParamsOnly = x[TransmissExpHistParamNum+NumOfStrains+1:TransmissExpHistParamNum+AgeSuscepParamNum]::Array{Float64,1}

	#Construct AgeSuscepByGroup from StrainSuscepParamsOnly &
	# AgeModifierSuscepParamsOnly vector
	#Row per susceptibility grouping, column per strain
	AgeSuscepByGroup_AgeUnmod = repeat(StrainSuscepParamsOnly,outer = (1,AgeGrpSuscepParamNum)) #Populate each column with age-specific suscep.

	#Row 1, youngest age band, has modifier of 1. Unaltered.
	#Update suscep. of other age bands using age modification factor
	#(transpose AgeModifierSuscepParamsOnly so multiply row by matching column element in AgeModifierSuscepParamsOnly')
	AgeSuscepByGroupTrans = [AgeSuscepByGroup_AgeUnmod[:,1] AgeSuscepByGroup_AgeUnmod[:,2:end].*AgeModifierSuscepParamsOnly']
	AgeSuscepByGroup = AgeSuscepByGroupTrans'

	return AgeSuscepByGroup, AgeSuscepParamNum
end

#Age group speicifc only (independent of strain)
function AgeOnlySuscepFn(x,AgeGrpSuscepParamNum,TransmissExpHistParamNum,NumOfStrains)

	AgeSuscepParamNum = AgeGrpSuscepParamNum::Int64
	AgeSuscepParamsOnly = x[TransmissExpHistParamNum+1:TransmissExpHistParamNum+AgeGrpSuscepParamNum]::Array{Float64,1}

	#Construct AgeSuscepByGroup from AgeSuscepParamsOnly
	#Row per susceptibility grouping, column per strain
	AgeSuscepByGroup = repeat(AgeSuscepParamsOnly,outer = (1,NumOfStrains)) #Populate each column with age-specific suscep.

	return AgeSuscepByGroup, AgeSuscepParamNum
end

#Age group speicifc only (independent of strain)
#Piecewise linear function (rather than a step function)
function PiecewiseLinearAgeSuscepFn(x,AgeGrpSuscepParamNum,TransmissExpHistParamNum,NumOfStrains)

	AgeSuscepParamNum = (AgeGrpSuscepParamNum + 1)::Int64

	#Assign susceptibility values (at knot ages) to variable
	AgeSuscepParamsVals = x[TransmissExpHistParamNum+1:TransmissExpHistParamNum+AgeSuscepParamNum]::Array{Float64,1}

	#Construct AgeSuscepByGroup from AgeSuscepParamsOnly
	#Row per susceptibility grouping, column per strain
	AgeSuscepByGroup = repeat(AgeSuscepParamsVals,outer = (1,NumOfStrains)) #Populate each column with age-specific suscep.

	return AgeSuscepByGroup, AgeSuscepParamNum
end

#-------------------------------------------------------------------------------
# (vii) FUNCTION TO RUN MODEL SIMULATION & PRODUCE DESIRED OUTPUTS TO FEED INTO
# SUMMARY STATISTIC FUNCTION
#-------------------------------------------------------------------------------

# Replicate FMAlt_NonCalibVersion.
function  RunModelFMAlt(FixedModelParams,x)
#Inputs:
#   FixedModelParams - Parameters with consistent values across all runs
#   x - Values of parameters that are being inferred

#Outputs:
#   SimnData - Model output to be fed into summary statistic function

#--------------------------------------------------------------------------
### Outline of steps
#--------------------------------------------------------------------------
# (i) Assign FixedModelParams to variables
# (ii) For parameters being inferred, update associated parameters
# (iii) Run ODE model and get end of seasons case counts
# (iv) For each age band, sum across relevant ages
# (v) Modify end of season counts by ascertainment probabilities
# (vi) Declare output variable names

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
	AscertainProbFn = FixedModelParams[19]
	FitAggAgeTuple = FixedModelParams[20]
	AgeGrpSuscepTuple = FixedModelParams[21]

	#----------------------------------------------------------------------
	### (ii) For parameters being inferred, update associated parameters
	#----------------------------------------------------------------------

	### Update R0 and then beta! ###
    R_0 = x[1:4]::Array{Float64,1}
	#gamma = InfectionParam[3]::Array{Float64,1} #rate of loss of infectiousness
    #beta = (gamma.*R_0)./spec_rad ##ecover strain transmission rates from contact structure & R_0
    InfectionParam[4] = R_0 #Update R_0 row in Infection Param (fourth row)

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
										ExpHistArrayParams,AgeGrpSuscepTuple);
	ExpHistArrayFnInputs = [ExpHistArray,ExpHistArrayParams]

	### Populate AgeSuscep array ###

	#Get age-binned susceptibility array
	AgeSuscepByGroup, AgeSuscepParamNum = AgeSuscepTypeFn(x,AgeGrpSuscepParamNum,TransmissExpHistParamNum,NumOfStrains)

    #Age Suscep by yr of age to be passed to flu model run function
    #Row per age band, column per strain
    AgeSuscep = zeros(M,NumOfStrains)

	if AgeSuscepTypeFn == PiecewiseLinearAgeSuscepFn  #Piecewise linear function

		#Get suscep values at knot ages
		SuscepValsKnotAges = AgeSuscepByGroup[:,1]

		#Iterate over each year of age ii
		for ii = 1:M

			#Convert loop index to year of age
			YrOfAgeVal = ii - 1

			#Check upper bound, first <= than.
			SuscepBandVectorIdx = findfirst(YrOfAgeVal .<= AgeSuscepUpperBounds)

			#Identify values at edge of bins, Entry idx and idx+1.
			LeftEdgeVal = SuscepValsKnotAges[SuscepBandVectorIdx]
			RightEdgeVal = SuscepValsKnotAges[SuscepBandVectorIdx+1]

			#For age ii, value is
			#Scale[idx] + (IdxInBin/NumYrsInBin)*(Scale[idx+1] - Scale[idx])
			if SuscepBandVectorIdx == length(AgeSuscepLowerBounds) #For last age category, get age years up to and including 100
				 SuscepBandSpan = 100 - AgeSuscepLowerBounds[SuscepBandVectorIdx]
			else
			  	SuscepBandSpan = AgeSuscepLowerBounds[SuscepBandVectorIdx+1] - AgeSuscepLowerBounds[SuscepBandVectorIdx]
			end

			IdxInSuscepBand = (ii-1) - AgeSuscepLowerBounds[SuscepBandVectorIdx] #Subtract 1 as ages begin from zero!
			SuscepValAgeAdjus = LeftEdgeVal + ((IdxInSuscepBand/SuscepBandSpan)*(RightEdgeVal - LeftEdgeVal))
			AgeSuscep[ii,:] .= SuscepValAgeAdjus

			#Error check
			if SuscepBandSpan == 0
				error("An age band contains only a single year of age. Incompatible with the selected piecewise linear susceptibility function.")
			end
		end
	else #Step function
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
	end
	#----------------------------------------------------------------------
    ### (iii) Run ODE model and get end of seasons case counts
    #----------------------------------------------------------------------

	#Run the model!
	StoreArrays = RunSeasonalFluModelFMAlt(ContactArray,MortalityFile,ONSPopnDistEngland,
											SimnParam,InfectionParam,AgeSuscep,AgeGrpSuscepTuple,
											VaccUptakeBySeason,LeakyVaccVarBySeason,LeakyTransFlag,
											ExpHistArrayFnInputs,ExpHistNum,ExpHistVaccType,MultiSeasonImmPropn,
											InfPropn_StartOfSeason,ICFromFile,SimnRunType)

	#Get timestep value & indicator value StoreFlag_PopnFOI
	timestep = SimnParam[8]
	StoreFlag_PopnFOI = SimnParam[9]

	#Disaggregate StoreArrays
	T = StoreArrays[1]::Array{Float64,1}
	C = StoreArrays[10]::Array{Float64,3}
	E_NotV = StoreArrays[4]::Array{Float64,3}
	E_V = StoreArrays[5]::Array{Float64,3}
	I_NotV = StoreArrays[6]::Array{Float64,3}
	I_V = StoreArrays[7]::Array{Float64,3}
	PopnDist = StoreArrays[11]::Array{Float64,2}

	if StoreFlag_PopnFOI == 1
		PopnFOI = StoreArrays[12]::Array{Float64,3}
	end

	#Visits in week (or month) t, denoted c_{t}, difference in the cumulative proportion
	#of consultations over the previous timestep
	#i.e. c_{t} = p(C(t)-C(t-1)),
	SliceVarIdx = convert(Int64,(365/timestep)+1) #Increment to access arrays at to obtain end of influenza season values
	CumulCaseCount = cat(C[1:SliceVarIdx:end,:,:],C[end:end,:,:]; dims = 1) #Get cumul. case count at specified intervals
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

				AgeBandInfAsPropnOvPopn = sum(UnNormCaseCountPerSeason_YrAgeGrps[:,StartIdx:end,:],dims=2)
				AgeBandPopnAsPropnOvPopn = sum(PopnDistBySeason[:,StartIdx:end,:],dims=2)


				#AgeBandedCaseRateTotalSum_FitCheckSeasons[:,ii,:] =
	             #   sum(StrainSeasonRateTotalSum_YrAgeGrps_FitCheckSeasons[:,StartIdx:end,:],2)
	        else #If not eldest age category, sum up to EndIdx.
				AgeBandInfAsPropnOvPopn = sum(UnNormCaseCountPerSeason_YrAgeGrps[:,StartIdx:EndIdx,:],dims=2)
				AgeBandPopnAsPropnOvPopn = sum(PopnDistBySeason[:,StartIdx:EndIdx,:],dims=2)
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
	SeasonRatePerAgeStrain = AscertainProbFn(AscertainProb,SeasonRatePerAgeStrain,
											AgeBandedCaseRateTotalSum_FitCheckSeasons,
											FitAggAgeFlag,AgeBandNum,
											AgeBandBounds,M)


	#----------------------------------------------------------------------
	### (vi) Declare output variable names
	#----------------------------------------------------------------------
	#Compute infected temporal profile
    #Stratified by age group, but aggregated over strains
	Total_I_AllDimn = sum(I_NotV,dims=3) + sum(I_V,dims=3) + sum(E_NotV,dims=3) + sum(E_V,dims=3) #Note, returns a three dimensional array, n x n x 1!

	#Remove singleton dimensions from Total_I_AllDimn
	Total_Infected = dropdims(Total_I_AllDimn; dims=3)

	#Assign outputs to be used in Summary Statisitc function to array
	#Vary amount of outputs. If tracked, output population level force of infection (FOI)
	if StoreFlag_PopnFOI == 1
		SeasonRatePerAgeStrainPerPerson = SeasonRatePerAgeStrain./100000 #Convert "per 100,000" rate to "per person" rate
		SimnData = [SeasonRatePerAgeStrainPerPerson,AgeBandedCaseRateTotalSum_FitCheckSeasons,PopnDist,PopnFOI]
	else
		SimnData = [SeasonRatePerAgeStrain,Total_Infected,PopnDist]
	end
	return SimnData::Array{Array{Float64},1}
end

#-------------------------------------------------------------------------------
# (viii) FUNCTIONS USED IN MODEL SIMULATION
#-------------------------------------------------------------------------------

function BuildExpHistArray(ExpHistVaccType,NumOfStrains,M,ExpHistArrayParams,AgeGrpSuscepTuple)
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
# MAIN LOOP - RUN APMC USING FMAlt FRAMEWORK
#-------------------------------------------------------------------------------
function APMC_FMAlt(RunID,ObvsData,SeasonsToSimulate,AscertainProbFn,ExpHistVaccType,
	FitAggAgeTuple,AgeGrpSuscepTuple,
    N_alpha,alpha,MinAcceptRate,MaxGen,PerturbVarScale,
    PriorFn,SummStatFn,SampleFromPriorFn,ModelSimnFn,FirstGenFromFileFlag)

	#-------------------------------------------------------------------------------
	# SET UP APMC SCHEME RELATED PARAMETERS
	#-------------------------------------------------------------------------------

	#Calculate number of samples before retention phase
	ScalingFactor = 1/alpha
	N = convert(Int64,round(N_alpha*ScalingFactor))

	#-------------------------------------------------------------------------------
	# DEFINE AND GROUP MODEL SIMULATION FIXED PARAMETERS
	#-------------------------------------------------------------------------------

	#-------------------------------------------------------------------------------
	# SPECIFY TYPE OF RUN THROUGH FLAG VARIABLE
	#-------------------------------------------------------------------------------
	# (INFLUENCES VACCINE UPTAKE/EFFICACY, & USE OF ODEBurnIn, ODEH1N1OnlyTime, ODEAlleStrainTime)
	#1 - exploratory; 2 - historical; 3 - alternative vacc. scheme
	SimnRunType = 2

	#------------------------------------------------------------------------------
	###  LOAD CONTACT DATA
	ContactArray = readdlm("../../../../Data/ContactData/UKAdeqContact_Over100Vers.csv",',')

	#------------------------------------------------------------------------------
	### IMPORT MORTALITY RATES
	#MortalityFile = "C:/Users/edwar/Documents/GitHub/DoH-FluVaccine/Data/DemographicData/MortalityRateByAge_ONS90plus.txt"
	MortalityFile = "../../../../Data/DemographicData/ModifiedData/MortalityProbPerAge0to100_EH.txt"

	#------------------------------------------------------------------------------
	### IMPORT INITIAL AGE DISTRIBUTION

	#Set proportion of individuals initially in each age
	#ONSPopnDistEngland20102017 = readdlm("../../../../Data/DemographicData/ONSPopnDistEngland20102017_0to100.txt",',')
	ONSPopnDistEngland20102018 = readdlm("../../../../Data/DemographicData/ONSPopnDistEngland20102018_0to100.txt",',')

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
	    ODEBurnIn = 0*365
	    ODEStaticPopnTime = 1*365
		ODEInferenceTime = (SeasonsToSimulate-1)*365;
	    ODEForwardSimnTime = 0*365
	elseif SimnRunType == 3

	    #Specify number of seasons that will use historical data
	    HistoricalSeasonNum = 7

	    ODEBurnIn = 0*365
	    ODEStaticPopnTime = 20*365
	    ODEInferenceTime = HistoricalSeasonNum*365
	    ODEForwardSimnTime = 4*365
	else
	    error("Incorrect RunType entered")
	end

	MaxTime = ODEBurnIn + ODEStaticPopnTime + ODEInferenceTime + ODEForwardSimnTime

	#Compute total time post burn-in
	ODESampleTime = ODEStaticPopnTime + ODEInferenceTime + ODEForwardSimnTime

	#Specify timestep to use in ode45 scheme
	timestep = 1.

	#Concatenate simulation variables
	SimnParam=[SimnStartDate,ODEBurnIn,ODEStaticPopnTime,ODEInferenceTime,ODEForwardSimnTime,
				MaxTime,ODESampleTime,timestep,StoreFlag_PopnFOI]
	#------------------------------------------------------------------------------

	#------------------------------------------------------------------------------
	###  TRANSMISSION RELATED PARAMETERS
	NumOfStrains = 4 #Specify number of strains in the system (A/H1N1, A/H3N2, two B lineages)
	gamma = 1/3.8*ones(NumOfStrains) #recovery rate
	sigma = [1/1.4,1/1.4,1/0.6,1/0.6] #latent rate

	#Calculate beta, group infectious status parameters into vector
	R_0 = [1.8,1.2,1.5,1.6]
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
	if SimnRunType == 1
	    #SYNTHETIC DATA

	    #Set vaccine uptake rate
	    #Entry (ii,kk)- Age class ii, day kk
	    NumVaccDays = 91
	    TargetCovPerSeason = ones(M,NumVaccDays)*((91/365)/3)
	    VaccUptakeBySeason = zeros(M,365) #Proportion given the vaccine during each day
	    VaccUptakeBySeason[:,244:334] .= TargetCovPerSeason/NumVaccDays

	elseif SimnRunType == 2
	    #BASED ON HISTORICAL DATA

	    #Cell per season
	    #Per cell, 2D array
	    #rows for age (0 to 90+), cols for calendar day of year
	    VaccUptakeBySeason = Vector{Array}(undef,SeasonsToSimulate)

	    # Import the data - Pandemic flu vacc (2009/2010 season)
		PandemicFluVaccUptake = XLSX.readdata("C:/Users/edwar/Documents/GitHub/DoH-FluVaccine/Data/VaccUptake/DailyWeeklyUptakeRates_HistoricalFluSeasons/EstimatedDailyUptakeData_Shifted/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc_EMH_May2019.xlsx","Both","C3:NC103")
	    #PandemicFluVaccUptake = XLSX.readdata("../../../../Data/VaccUptake/DailyWeeklyUptakeRates_HistoricalFluSeasons/EstimatedDailyUptakeData_Shifted/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc_EMH_May2019.xlsx","Both","C3:NC103")

	    # Collate into Array, Assign to storage cell
	    VaccUptakeBySeason[1] = Array{Float64, 2}(PandemicFluVaccUptake)

	    #Import the data - Sesonal flu vacc
	    SheetNames = ["2010_2011","2011_2012","2012_2013","2013_2014",
	        "2014_2015","2015_2016","2016_2017","2017_2018"]

	    for ii = 1:length(SheetNames)
			HistoricalSeasonalFluVaccUptake = XLSX.readdata("C:/Users/edwar/Documents/GitHub/DoH-FluVaccine/Data/VaccUptake/DailyWeeklyUptakeRates_HistoricalFluSeasons/EstimatedDailyUptakeData_Shifted/AgeStrucModel_ByYrOfAge_DailyVaccUptakeBySeasonCalYr_All_EMH_May2019.xlsx","$(SheetNames[ii])","C3:NC103")
	        #HistoricalSeasonalFluVaccUptake = XLSX.readdata("../../../../Data/VaccUptake/DailyWeeklyUptakeRates_HistoricalFluSeasons/EstimatedDailyUptakeData_Shifted/AgeStrucModel_ByYrOfAge_DailyVaccUptakeBySeasonCalYr_All_EMH_May2019","$(SheetNames[ii])","C3:NC103")

	        # Collate into Array, Assign to storage cell
	        VaccUptakeBySeason[ii+1] = Array{Float64, 2}(HistoricalSeasonalFluVaccUptake)
	            #Add 1 to idx as first entry is for 2009/2010 season
	    end

	elseif SimnRunType == 3
	    #RUN ALTERNATIVE VACC. SCHEMES

	    #Cell per season
	    #Per cell, 2D array
	    #rows for age (0 to 90+), cols for calendar day of year
		VaccUptakeBySeason = Vector{Array}(undef,SeasonsToSimulate)

	    # Import the data - Pandemic flu vacc (2009/2010 season)
		PandemicFluVaccUptake = XLSX.readdata("C:/Users/edwar/Documents/GitHub/DoH-FluVaccine/Data/VaccUptake/DailyWeeklyUptakeRates_HistoricalFluSeasons/EstimatedDailyUptakeData_Shifted/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc_EMH_May2019.xlsx","Both","C3:NC103")
	    #PandemicFluVaccUptake = XLSX.readdata("../../../../Data/VaccUptake/DailyWeeklyUptakeRates_HistoricalFluSeasons/EstimatedDailyUptakeData_Shifted/AgeStrucModel_ByYrOfAge_DailyVaccUptakeCalYr_PandemicFluVacc_EMH_May2019.xlsx","Both","C3:NC103")

	    # Collate into Array, Assign to storage cell
	    VaccUptakeBySeason[1] = Array{Float64, 2}(PandemicFluVaccUptake)

	    #Import the data - Sesonal flu vacc
	    SheetNames = ["2010_2011","2011_2012","2012_2013","2013_2014",
	        "2014_2015","2015_2016","2016_2017","2017_2018"]

	    for ii = 1:length(SheetNames)
			HistoricalSeasonalFluVaccUptake = XLSX.readdata("C:/Users/edwar/Documents/GitHub/DoH-FluVaccine/Data/VaccUptake/DailyWeeklyUptakeRates_HistoricalFluSeasons/EstimatedDailyUptakeData_Shifted/AgeStrucModel_ByYrOfAge_DailyVaccUptakeBySeasonCalYr_All_EMH_May2019.xlsx","$(SheetNames[ii])","C3:NC103")
	        #HistoricalSeasonalFluVaccUptake = XLSX.readdata("../../../../Data/VaccUptake/DailyWeeklyUptakeRates_HistoricalFluSeasons/EstimatedDailyUptakeData_Shifted/AgeStrucModel_ByYrOfAge_DailyVaccUptakeBySeasonCalYr_All_EMH_May2019","$(SheetNames[ii])","C3:NC103")

	        # Collate into Array, Assign to storage cell
	        VaccUptakeBySeason[ii+1] = Array{Float64, 2}(HistoricalSeasonalFluVaccUptake)
	            #Add 1 to idx as first entry is for 2009/2010 season
	    end
	else
	    error("Incorrect RunType entered")
	end

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
	        VaccEffTemp = XLSX.readdata("../../../../Data/VaccEfficacy/VaccEfficacy_ByYrOfAge.xlsx","$(SheetNames[ii])","C3:F103")

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
	        VaccEffTemp = XLSX.readdata("../../../../Data/VaccEfficacy/VaccEfficacy_ByYrOfAge.xlsx","$(SheetNames[ii])","C3:F103")

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
	RetainedSeasonNum = size(ObvsData,1)

	#--------------------------------------------------------------------------
 	### AGGREGATE FIXED PARAMETERS
	#--------------------------------------------------------------------------
	FixedModelParams = [ContactArray,MortalityFile,ONSPopnDistEngland20102018,
	                    SimnRunType,ExpHistVaccType,StoreFlag_PopnFOI,
	                    SimnParam,NumOfStrains,ExpHistNum,M,InfectionParam,
	                    MultiSeasonImmPropn,VaccUptakeBySeason,LeakyTransFlag,
	                    LeakyVaccVarBySeason,InfPropn_StartOfSeason,ICFromFile,
	                    RetainedSeasonNum,AscertainProbFn,
	                    FitAggAgeTuple,AgeGrpSuscepTuple]

	#--------------------------------------------------------------------------
	# SET END-OF-GENERATION OUTPUT TEXT FILE INFO
	#--------------------------------------------------------------------------
	OutputFileName = "FMAlt_APMCOutputFilesJULIA/EndOfGenAPMC_FMAlt"

	#--------------------------------------------------------------------------
	# RUN APMC SCHEME
	#--------------------------------------------------------------------------
	RetainedParams,RetainedWeights,RetainedSummStat,GenT = APMC_loop(ObvsData,FirstGenFromFileFlag,N,alpha,N_alpha,MinAcceptRate,MaxGen,PerturbVarScale,
	    PriorFn,SummStatFn,SampleFromPriorFn,FixedModelParams,ModelSimnFn,RunID,OutputFileName)

	#-------------------------------------------------------------------------------
	# SAVE TO FILE
	#-------------------------------------------------------------------------------

	#Specify filename to write sample values
	FName_GenT = "FMAlt_APMCOutputFilesJULIA/APMCsamples_FMAlt_JulRun#$(RunID)_GenT.txt"
	FName_RetainedParams = "FMAlt_APMCOutputFilesJULIA/APMCsamples_FMAlt_JulRun#$(RunID)_RetainedParams.txt"
	FName_RetainedSummStat = "FMAlt_APMCOutputFilesJULIA/APMCsamples_FMAlt_JulRun#$(RunID)_RetainedSummStat.txt"
	FName_RetainedWeights = "FMAlt_APMCOutputFilesJULIA/APMCsamples_FMAlt_JulRun#$(RunID)_RetainedWeights.txt"

	#--------------------------------------------------------------------------
	### SAVE SUMMARY STATISTICS TO FILE
	#--------------------------------------------------------------------------
	writedlm(FName_GenT,GenT)
	writedlm(FName_RetainedParams,RetainedParams)
	writedlm(FName_RetainedSummStat,RetainedSummStat)
	writedlm(FName_RetainedWeights,RetainedWeights)

end
