#### Program information ######################################################
# Package: MenA epi dynamics                                                  #
# Source file name: MenA_VaccSims.R                                           #
# Version Date 01/21/2020                                                     #
#_____________________________________________________________________________#
# Purpose: Top-level simulation program. Runs multiple iterations of the      #
# MenA simulations under a specified set of inputs.                           #
#_____________________________________________________________________________#
# Inputs: Script builds some inputs from source files on population size and  #
# posterior parameter estimates.                                              #
#_____________________________________________________________________________#
# Outputs:                                                                    #
#_____________________________________________________________________________#
# Steps in this program:                                                      #
# (1) Set up libraries                                                        #
# (2) User-set parameters                                                     #
# (3) Import and format data/functions                                        #
# (4) Run simulations                                                         #
#_____________________________________________________________________________#
# Author: Mike Jackson                                                        #
# This code is adapted from the VIMC MenA vaccination models, so much of the  #
# structure and many of the inputs are derived from that project.             #
#_____________________________________________________________________________#

## TO DO:
# 1. Create function to get desired outputs
# 2. Create set of inputs for Emergent Pathogen models
# 3. Adapt top-level program to use these inputs and function


### (1) Set up libraries ######################################################

library(lubridate)
library(doParallel)
library(dplyr)
library(data.table)
library(reshape2)
library(tidyr)

### (2) User-set parameters ###################################################
# Edit the code in this section to specify country and scenario.              #
# Update 2019.12.20 - include option to pseudo-automate. If set to true, the  #
# program will find the first country/scenario combination that has not yet   #
# been completed and will run it, from the scenario_tracker.csv file. If set  #
# to false, manually enter the country, program option, and region type.      #

begin <- Sys.time()
start <- as.Date("1951-01-01") # Use 1/1/1951 to give 50 years burn-in
end <- as.Date("2050-12-31")   # Use 50 years of simulation time
PSA <- TRUE
seed <- 6543  # Seed for random sto, use same for all scenarios
nSims <- 200  # Update: 100 takes around 12 minutes if using 4 cores.
phi <- 0.2 # Set stochastic parameter


# Directory containing population inputs
input.dir  <- "C:/Users/O992928/Documents/MenA epi dynamics/Data/Input data"
# Directory for simulation outputs
output.dir <- "C:/Users/O992928/Documents/MenA epi dynamics/Data/Sim output"
# Directory containing R scripts
script.dir <- "C:/Users/O992928/documents/MenA epi dynamics/MenA_Epi_Code"

mycountry <- "BFA"
myregion <- "hyper"    
# Vaccination program is legacy of original code use in VIMC
# For this project, we have no vaccination, so set to "none"
vacc_program <- "none" 
vacc_subprogram <- "default"

### (3) Import and format data/functions ######################################
# Import scripts used in the simulations, country-specific parameters, and    #
# vaccination program details.                                                #

# (A) Import scripts
if (dir.exists(script.dir)) {
  if (file.exists(paste0(script.dir, "/","MenA_paramCheck.R"))==FALSE) {
    msg<(paste0("MenA_paramcheck.R not found in ", script.dir))
    stop(msg)
  }
} else {
  script.dir<- getSrcDirectory(function(dummy) {dummy})
  if (file.exists(paste0(script.dir, "/","MenA_paramCheck.R"))==FALSE) {
    print(paste0("MenA_paramCheck.R not found in ", script.dir))
    stop("This script requires 6 other scripts; please put in same directory as this one or specify script directory.")
  }
}
setwd(script.dir)
source("MenA_paramCheck.R")
source("ModelInputUtilities.R")
source("MenA_OneSim.R")
source("MenA_helper_functions.R")
source("MenA_summarization_functions.R")


# (B) Check parameters set above
setparams <- as.list(c(mycountry, as.character(start), as.character(end), myregion, PSA, vacc_program, phi, seed, nSims, input.dir, output.dir))
names(setparams) <- c("mycountry", "start", "end", "myregion", "PSA", "vacc_program", "phi", "seed", "nSims", "input.dir", "output.dir")
if (CheckSetParameters(setparams)==FALSE) {
    stop(spmessage)
} else {
  if (length(spmessage)>1) { print(spmessage) }
}

# (C) Import country-specific parameters
myparams.full <- GetDemographicParameters(path=input.dir,  mycountry=mycountry, start=start, end=end)
if (CheckDemogParameters(myparams.full)==FALSE) {
  stop(dpmessage)
} else {
  if (length(dpmessage)>1) { print(dpmessage) }
}

# (D) Import vaccination program details
if (vacc_program!="none") {
  myvacc <- GetVaccScenario(mycountry=mycountry, scenario=vacc_program, sub.scenario=vacc_subprogram, directory=input.dir)
  if (is.data.frame(myvacc)==FALSE) { stop(vaccmsg)}  #check for output
  #make as vector of years where nothing happens (empty except for campaign only) for efficiency
  if (vacc_program=="campaign") {
    nodoses <- as.vector(myvacc[is.na(myvacc$DosesCampaign) | myvacc$DosesCampaign==0,"year"])
  }
}

# (E) Country-specific life expectancy
my.lifex <- GetLifeExp(path=input.dir, mycountry.s=mycountry)

# (F) Read in parameters calculated in ABC, or a row of parameters to be used by ABC.  
paramfixed <- GetModelParams(path=script.dir, region.val=myregion)

# (G) Initialize full population.
startSize <- myparams.full[myparams.full$year==year(start)-1, "totalpop"]
initpop.full <- InitializePopulation(scriptdir=script.dir, inputdir=input.dir, start=start, end=end, country=mycountry, region=myregion, startSize=startSize)
#check for errors
if (!(is.numeric(initpop.full))) {
  if (disterr!="") { print(disterr) } 
  if (dxerr!="") { print(dxerr) } 
  stop(initmsg)
}

# (H) Scale down the full population to the modeled population
# This means changing both the starting population size and
# the annual number of births.
pct.modeled <- GetModelPct(path=input.dir, mycountry.s=mycountry)
initpop <- initpop.full
initpop[,,1] <- initpop[,,1] * pct.modeled

myparams <- myparams.full
myparams$births <- myparams$births * pct.modeled

### (4) Run simulations #######################################################
# Set up random seed vector and use parallel processing to run simulations.   #

summarizeme <- 1
# Set the random number seed based on seed, then create a vector of random number seeds (consistent within seed values)
set.seed(seed, kind = NULL, normal.kind = NULL)
seed.vec <- unique(floor(runif(nSims*2, 0, 1000000)))[1:nSims]

# Begin simulations
cl <- makeCluster(4)  #scale this upwards if you're on a workstation with >16gb memory
registerDoParallel(cl)
my_data <- foreach(n=1:nSims, .verbose=TRUE, .packages = c("lubridate", "dplyr", "data.table", "reshape2")) %dopar% {
  set.seed(seed.vec[n])
  paramfixed.row <- paramfixed[n,]  
  finalpop<-MenASimulation(startdt=start, enddt=end, fp=paramfixed.row, initpop=initpop, vacc_program=vacc_program,
                           countryparams=myparams, region=myregion, country=mycountry)
  if (summarizeme > 0) {
    # Using PSA option for CFR
    cfr <- as.numeric(paramfixed.row[, c("cfr1", "cfr2", "cfr3", "cfr4", "cfr5", "cfr6")])
    summarizeOneSim(finalpop, n, cfr)
    } #end of conditional summarization
} #end of foreach loop
stopCluster(cl)


# Also run one time to get the cohort size
# VIMC wants the size of the full population, not just the vaccine-
# targetted population, so run this without scaling down
onerun <- MenASimulation(startdt=start, enddt=end, fp=paramfixed[4,], initpop=initpop.full, vacc_program=vacc_program,
                         countryparams=myparams.full, region=myregion, country=mycountry)
cohortSize <- getCohortSize(onerun)
totalPop <- cohortSize %>% 
  group_by(year) %>% summarize(tot=sum(cohortsize))


### (5) Create outputs ########################################################
# Calculate overall annual incidence, identify frequency of major epidemics,  #
# and get inter-epidemic incidence.                                           #

results.df <- data.frame(sim=1:nSims, epi.freq=numeric(nSims), 
                         nonepi.inc=numeric(nSims))
 
for (s in 1:nSims){
  # Get annual incidence
  cases.yr <- aggregate(Cases ~ year, data=my_data[[s]], FUN=sum)[,2] 
  Inc.yr <- 100000*cases.yr/totalPop$tot
  
  # Classify years as major epidemic or not
  epi.yr <- Inc.yr >= 100
  # Drop years that are part of multi-year epidemic and not the first year
  epi.shift <- c(0, epi.yr[1:(length(epi.yr)-1)])
  epi.clean <- ifelse(epi.yr==1 & epi.shift==1, NA, epi.yr)
  
  # Count epidemics, convert to average inter-epidemic period
  epi.count <- sum(epi.clean, na.rm=TRUE)
  results.df$epi.freq[s] <- length(epi.clean[is.na(epi.clean)==FALSE]) / epi.count
  
  # Average incidence in inter-epidemic periods
  results.df$nonepi.inc[s] <- mean(Inc.yr[Inc.yr < 100])
  
}

write.csv(results.df, file="C:/Users/O992928/Desktop/Sim_results.csv")