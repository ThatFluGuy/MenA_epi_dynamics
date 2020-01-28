# MenA_epi_dynamics

## Background
In our simulation models of the impact of vaccination on incidence of serogroup A *Neisseria meningitidis* (NmA), we assume the pathogen is unchanged over time. Under this assumption, under a set of estimated waning parameters and assumed stochasticity of infection, we can reproduce the dynamics of NmA epidemics in the African Meningitis Belt. In particular, we see major epidemics (incidence >100 cases per 100,000 population) every 7-10 years, with minor epidemics in between.

However, one feature of meningococcal population dynamics is the decade-scale replacement of previously circulating NmA strains with newly emergent strains. For example, the ST-7 complex emerged in the early 1990s and displaced the previously-circulating ST-5 strains,[Caugant *PLoS One* 2012; Caugant *Vaccine* 2007; Nicolas *J Clin Microbio* 2005]. The ST-5 complex in turn had emerged in China in the 1960s and migrated to Africa in by the 1990s.[Nicolas *J Clin Microbio* 2000] This suggests that changes to the pathogen may be important drivers of NmA epidemiology in the pre-vaccine era. Specifically, new strains with novel sub-capsular antigents could evade pre-existing host immunity and launch new epidemics.

This difference is consequential even with widespread use of serogroup A meningococcal conjugate vaccines (MenA) in Africa. Mathematical models predicting the impact and cost-benefit of MenA vaccination program assume (counterfactual) dynamics in the absence of vaccination that are driven by host immunity. If major epidemics are driven by pathogen changes, these models may not correctly estimate the impact of MenA vaccination. 

This project explores ranges of parameters for immune waning and other simulation model parameters to recreate NmA dynamics assuming the periodic emergence of antigenically novel NmA strains.

## Methods
Unless otherwise specified, all modeling details are as described in Jackson et al *PLoS One* 2018. Informally, that was "version 2.0" of the KPWA MenA simulation model, which used approximate Bayesian computation (ABC) to estimate posterior distributions for all model parameters.

### Population structure
In this project we use UN data to define various population inputs, including population size, age structure, birth rates, and death rates. For simplicity, we use the versions of these data as compiled by VIMC in the 2019 modeling round.

### Immune dynamics model
The "immune dynamics" model assumes that NmA epidemics are driven solely by fluctations in population-level immunity to NmA colonization and disease. This is the model used in Jackson *PLoS One* 2018. In Jackson *PLoS One* 2018, we identified 113 simulated parameter sets that met the model fit threshold and were included in the posterior distribution. Here, following our approach for the 2019 VIMC deliverables, we expand this to 200 parameter sets by randomly sampling an additional 87 sets (with replacement) from the 113. 

### Emergent pathogen models
The "emergent pathogen" models assume that antigenically novel NmA strains emerge periodically, which escape pre-existing human host immunity. Based on observed data that suggests three major waves of hyperinvasive NmA strains from ~1960 through 2000, we assume stochastic emergence of new strains, with a mean time-to-emergence of 16.67 years. As a method of implementing this, we first assume that immunity to NmA colonization is longer lasting than in the Immune Dynamics model. We range this from lifelong down to the duration used in the Immune Dynamics model. Second, we (simplistically) assume that the proportion of the population in the low-immunity and high-immunity states drops by factor of x (0 < x < 1) at the time of emergence.

We will run multiple emergent pathogen models with different values of x and different durations of protection.

### Model evaluation
We run 200 iterations of each model and calculate the following quantities on each iteration:
a) Number of major epidemics (>100 cases per 100,000 population);
b) Mean time between major epidemics;
c) Mean incidence during years without major epidemics;
d) Mean incidence during years with major epidemics.

We first verify that the immune dynamics model recreates expected epidemic size/frequency.
From among the various emergent pathogen models, we assess which of them (if any) recreate expected epidemic size/frequency. We then describe parameter ranges that reproduce this behavior in the emergent pathogen models.
