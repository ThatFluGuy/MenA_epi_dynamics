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

useNewStrains.top <- TRUE
newStrainYears.top <- c(2001, 2022, 2034)
newStrainDrop.top <- 1
newStrainDur.top <- 3
# Maximum duration of immunity (wh4) is ~36 years
# Try keeping newStrainDur <= 3 so immunity never exceeds 100 years


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
nSims <- 100  # Update: 100 takes around 12 minutes if using 4 cores.
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

# (I) If using the Emergent Pathogen model, change the duration 
# of Hs/Hc immunity per the described scaling factor

if (useNewStrains.top==TRUE){
  paramfixed$wh1 <- paramfixed$wh1*newStrainDur.top
  paramfixed$wh2 <- paramfixed$wh2*newStrainDur.top
  paramfixed$wh3 <- paramfixed$wh3*newStrainDur.top
  paramfixed$wh4 <- paramfixed$wh4*newStrainDur.top

}

### (4) Run simulations #######################################################
# Set up random seed vector and use parallel processing to run simulations.   #
  
summarizeme <- 1
# Set the random number seed based on seed, then create a vector of random number seeds (consistent within seed values)
set.seed(seed, kind = NULL, normal.kind = NULL)
seed.vec <- unique(floor(runif(nSims*2, 0, 1000000)))[1:nSims]

# Begin simulations
cl <- makeCluster(2)  #scale this upwards if you're on a workstation with >16gb memory
registerDoParallel(cl)
my_data <- foreach(n=1:nSims, .verbose=TRUE, .packages = c("lubridate", "dplyr", "data.table", "reshape2")) %dopar% {
  set.seed(seed.vec[n])
  paramfixed.row <- paramfixed[n,]  
  finalpop<-MenASimulation(startdt=start, enddt=end, fp=paramfixed.row, initpop=initpop, vacc_program=vacc_program,
                           countryparams=myparams, region=myregion, country=mycountry,
                           useNewStrains=useNewStrains.top, newStrainYears=newStrainYears.top,
                           newStrainDrop=newStrainDrop.top)
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
# Need to add calculation of age-specific prevalence                          #

results.inc <- data.frame(sim=1:nSims, epi.freq=numeric(nSims), 
                         nonepi.inc=numeric(nSims))

prev.compile <- data.frame(AgeGroup=numeric(0), epi.yr=logical(0),
                           Prev=numeric(0))

incid.compile <- data.frame(epi.yr=character(0), AgeGroup=numeric(0), 
                            Incid=numeric(0), stringsAsFactors = FALSE)

 
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
  results.inc$epi.freq[s] <- length(epi.clean[is.na(epi.clean)==FALSE]) / epi.count
  
  # Average incidence in inter-epidemic periods
  results.inc$nonepi.inc[s] <- mean(Inc.yr[Inc.yr < 100])
  
  # Classify years
  cases.yr <- aggregate(Cases ~ year, data=my_data[[s]], FUN=sum) 
  cases.yr$Inc.yr <- 100000*cases.yr$Cases / totalPop$tot
  cases.yr$epi.yr <- case_when(
    cases.yr$Inc.yr < 20 ~ "minor",
    cases.yr$Inc.yr < 100 ~ "middle",
    cases.yr$Inc.yr >= 100 ~ "major"
  )
  
  # Incidence by year and age group
  inc.df <- left_join(my_data[[s]], cohortSize, by=c("year", "AgeInYears"))
  
  inc.df$AgeGroup <- case_when(
    inc.df$AgeInYears < 1 ~ 0,
    inc.df$AgeInYears < 5 ~ 1,
    inc.df$AgeInYears < 10 ~ 5,
    inc.df$AgeInYears < 15 ~ 10,
    inc.df$AgeInYears < 20 ~ 15,
    inc.df$AgeInYears < 25 ~ 20,
    inc.df$AgeInYears < 30 ~ 25,
    TRUE ~ 30
  )
  
  inc.df2 <- inc.df %>% group_by(year, AgeGroup) %>%
    summarize(Cases=sum(Cases), cohortsize=sum(cohortsize))
  inc.df2$Incid <- 100000 * inc.df2$Cases / inc.df2$cohortsize
  
  inc.df3 <- left_join(inc.df2[, c("year", "AgeGroup", "Incid")], 
                       cases.yr[, c("year", "epi.yr")], by="year")
  
  incid.compile <- rbind(incid.compile, as.data.frame(inc.df3))
  
  
  # Get age-specific prevalence by type of year (major epidemic vs not)
  # First add in cohort sizes
  prev.df <- left_join(my_data[[s]], cohortSize, by=c("year", "AgeInYears"))
  
  prev.df$AgeGroup <- case_when(
    prev.df$AgeInYears < 1 ~ 0,
    prev.df$AgeInYears < 5 ~ 1,
    prev.df$AgeInYears < 10 ~ 5,
    prev.df$AgeInYears < 15 ~ 10,
    prev.df$AgeInYears < 20 ~ 15,
    prev.df$AgeInYears < 25 ~ 20,
    prev.df$AgeInYears < 30 ~ 25,
    TRUE ~ 30
  )
  
  # Sum carriers and cohort sizes by age group
  prev.df2 <- prev.df %>% group_by(year, AgeGroup) %>%
    summarize(Carriers=sum(Carriers), cohortsize=sum(cohortsize))
  
  # Calculate year- and agegroup-specific prevalence
  prev.df2$Prev <- prev.df2$Carriers / prev.df2$cohortsize
  
  # Assign years as major epidemic vs not
  prev.df3 <- left_join(prev.df2, data.frame(year=2000:2050, epi.yr), by="year")

  prev.compile <- rbind(prev.compile, 
                        as.data.frame(prev.df3[, c("AgeGroup", "epi.yr", "Prev")]))
  
}

write.csv(results.inc, file="C:/Users/O992928/Desktop/Sim_incidence.csv",
          row.names=FALSE)

# Incidence by age and epidemic type
inc.df4 <- incid.compile %>% group_by(epi.yr, AgeGroup) %>%
  summarize(meanInc=mean(Incid), stdInc=sd(Incid))

write.csv(inc.df4, file="C:/Users/O992928/Desktop/Sim_incid.csv", row.names=FALSE)

# Get mean and std prevalence by age group and epi.yr
prev.df4 <- prev.compile %>% group_by(AgeGroup, epi.yr) %>%
  summarize(meanPrev=mean(Prev), stdPrev=sd(Prev), N=length(Prev)) %>%
  arrange(epi.yr, AgeGroup)

write.csv(prev.df4, file="C:/Users/O992928/Desktop/Sim_prevalence.csv", row.names=FALSE)