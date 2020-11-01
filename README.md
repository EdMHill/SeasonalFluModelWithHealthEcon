# SeasonalFluModelWithHealthEcon
This repository houses code developed for the analysis presented in the scientific publication "Optimising age coverage of seasonal influenza vaccination in England: A mathematical and health economic evaluation" by Edward M. Hill, Stavros Petrou, Henry Forster, Simon de Lusignan, Ivelina Yonova, and Matt J. Keeling.  

Publication details: Hill et al. (2020) Optimising age coverage of seasonal influenza vaccination in England: A mathematical and health economic evaluation. 16(10): e1008278. doi: 10.1371/journal.pcbi.1008278. URL: https://doi.org/10.1371/journal.pcbi.1008278.

Please find below an explainer of the directory structure within this repository.  

# Data

The GP consultation data and Hospital Episode Statistics (HES) data contain confidential information, with public data deposition non-permissible for socioeconomic reasons. The GP consultation data resides with the RCGP Research and Surveillance Centre and are available via the RCGP RSC website (www.rcgp.org.uk/rsc). The HES database resides with NHS Digital and are available via the HES webpage (https://digital.nhs.uk/data-and-information/data-tools-and-services/data-services/hospital-episode-statistics).

## ContactData
Contact patterns by age group. Sourced from the study "Inferring the Structure of Social Contacts from Demographic Data in the Analysis of Infectious Diseases Spread" by L. Fumanelli et al. PLoS Comput Biol 8(9): e1002673. doi: 10.1371/journal.pcbi.1002673.

## DemographicData
Population size and mortality data (from ONS)

## HealthEconParamData
Summary data on the cost and relative rate of occurrence (compared to GP consultations) for mortality, inpatient and outpatient attendances.

## ILIData
Using RCGP data, compute cumulative ILI cases (sum of weekly ILI rate per 100,000) for 2009/10 - 2017/2018 seasons.  
Non-age & age-structured model variants.  
Values are scaled using weekly influenza positivity data.

## RCGPSamplePositivity
Construct arrays storing weekly influenza positivity values.     
Monitored through the RCGP sentinel swabbing scheme in England.

## RiskGrpPropnsData
By single year age groups, the proportion of the population of the given age that were in a seasonal influenza "at-risk' group. Season-by-season values spanning the 2000/01 - 2017/2018 influenza seasons.

## VaccEfficacy
Population-level vaccine efficacy estimates.

## VaccUptake
Daily vaccine uptake proportions for at-risk, low risk, entire population groups under differing vaccination scenarios.

## WHOFluNet
Data on circulating influenza virus composition from World Health Oragnisation FluNet database

# Results

## EvaluatingModelFit
Script to produce bar plots containing outputs from 100 model simulation replicates.

## FMAlt_OptimFitOutput
Contains scripts to analyse outputs from the parameter optimisation scheme.

## HealthEconEval
Threshold vaccine dose prices for simulated vaccination strategies.

# src

## HealthEconEval
Scripts to carry out health economic evaluation.

## ParamOptim
Fit the transmission model to the data using particle swarm optimisation.

## PosteriorPredictiveModelSimns
Run simulations of the model using parameter sets representing samples from the posterior distribution.

## ModelFns
Functions to numerically solve ODEs and update exposure history array
