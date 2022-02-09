# NetworkCOVID19

<!-- badges: start -->

<!-- badges: end -->

The goal of NetworkCOVID19 is to perform contact network based simulation of COVID-19 outbreaks. The orginal code are designed for studies described in "**Modelling the impact of rapid tests, tracing and distancing in lower-income countries suggest optimal policies varies with rural-urban settings**": https://www.medrxiv.org/content/10.1101/2021.03.17.21253853v1

## Installation

You can install the development version of NetworkCOVID19 from [GitHub](https://github.com/) with:

```r
# install.packages("devtools")
devtools::install_github("Xilin-Jiang/NetworkCOVID19")
```

## Simple breakdown example

The simulation is breaking into three steps for easy nevigation: first initialise the parameters; 
second generate and save networks, which is the synthetic populations containing the contact information;
third simulate the outbreaks. 

The first step is initialising all the parameters. Setting = 2 will specify an urban setting, 
social_distancing_flg = 1 refers to no social distancing. You could specify your own contact number
which will change the R0. 
```r
library(NetworkCOVID19)
para <- init_para(setting = 2, country_id = 1, social_distancing_flg = 1, contact_number = NA)
```

After initialising the parameter, you could constumise your own simulation setting by modifying para. 
A full explanation of parameter setting is in **Customise your simulation settings**. In order to specify
your own parameter, you should input a list of parameters into the init_para function. 

```
para <- list()
# use a smaller R0, for a population that are with high immunity level
para$R0 <- 1.3
# use a high number of daily contact number, to reflect densely populationed area
para$num_cc <- 20
# generate ERGM based contact network
nw <- network_generate(para)
# simulate transmission using the network generated
sim_rslt <- simulate_transmission(NW_SIM = nw)
```

## Full scale simulation
Following part require larger computational power. We show the wrapper function to simulate many epidemics under different
containment strategies.
```
results <- network_covid_simulate(rep_num = 1, network_num = 1, output = "example", para = NA)
plt <- plot_epidemic_curves(results, title_fig = "")
```


## Customise your simulation settings
Below are the code for customise your only simulation parameters; please use this command then followed by the init_para function. 
```r
  para <- list()
  para$setting <- setting # 1=rural 2=affluent 3=slum
  # Environment setup
  para$pop_sz <- 1000 # 5000
  para$community.pop_sz <- 200 # use this small network to construct large network
  para$Non_HH_CC_rate <- 1 # .8,.6,.4,.2 -- if set to 0.2, means only 20% of non-household contact is happening due to controling method

  ##########################################
  # loading country specific variables:
  # age distribution & household
  ##########################################
  para$age_dist <- c(0.481, 0.203, 0.316) # Uganda
  para$HH_dist <- c(.11, .22, .27, .4) # Uganda

  ##########################################
  # processing demographic information;
  # for details: refer to the supplementary methods
  ##########################################
  para$num_cc <- 7 # set daily close contact to be 7
  para$family_sz <- 5 # average household size 5
  # set the percentage of HH_cc
  para$percent_HH_cc <- .5

  # load the contact structure
  para$age_mix <- contact_all[[country_id]]
  para$Home_age_mix <- contact_home[[country_id]]
  # adjust the age matrix to represent the specified household contact rate
  para$age_mix <- para$Home_age_mix + (para$age_mix - para$Home_age_mix) *
    sum(para$Home_age_mix)/sum(para$age_mix - para$Home_age_mix) *  (1-para$percent_HH_cc)/para$percent_HH_cc

  # adjust R0 for 1) young people susceptibility 2) subclinical cases
  # adjust for the social distancing
  para$num_cc_scdst <- para$num_cc * ((1-para$percent_HH_cc)*para$Non_HH_CC_rate + para$percent_HH_cc) # reduce the number of cc
  para$age_mix_scdst <- para$Home_age_mix + (para$age_mix - para$Home_age_mix) * para$Non_HH_CC_rate
  para$percent_HH_cc_scdst <- para$percent_HH_cc/((1-para$percent_HH_cc)*para$Non_HH_CC_rate + para$percent_HH_cc)

  ###############################
  # parameters relevant to transmission
  ###############################

  para$E_0 <- 2 # number of initial importation
  para$infect_sz <- (-1)*para$pop_sz/1000 # intervention only starts after X per 1000 individuals are already infected

  R_UK <- 2.7
  num_cc_UK <- 13
  R0_baseline <- R_UK/num_cc_UK # the R0 is measured in affluent region
  para$R0 <- R0_baseline * para$num_cc # R naught
  
  # Transmission parameter & tool availability
  para$symptom_report <- 0.7 # precentage infected report symptom
  para$theta <- 0.1 # quarantine effect: precentage of remaining transmission quarantined individual have
  para$pcr_available <- 1000*para$pop_sz/1000 # -1 # daily maximum PCR tests per 1000 population.
  para$ppe_coef <- 1 # if people wear ppe, then their transmissibility will be different outside their family


  # subclinical rate
  para$sub_clini_rate <- 0.3
  para$asym_rate <- 0.2 # transmission rate of asymptomatic patient
  para$death_rate <- c(0.002,0.002,0.01) # death rate for three age groups
  para$death_delay <- 10 # delay time for a death to be recorded

  # Parameters about Infected pts
  para$ab_test_rate <- 0.7 # % accepted RDT testing among detected (symptomatic) patients
  para$pcr_test_rate <- 0.8 # % accepted pcr testing among detected (symptomatic) patients

  para$onsetiso <- 0.2 # Isolation compliance rate at onset based on symptom
  para$abiso <-0.9 # Isolation compliance rate among RDT testing positive
  para$pcriso <-0.9 # Isolation compliance rate among pcr testing positive

  para$delay_symptom <- 1 # days of delay after onset to detect symptomatic patients
  para$delay_ab <- 8 # -- not so related with the current strategy days of delay after onset to receive antibody test and obtain result; 
  para$delay_pcr <- 4 # days of delay after onset to report pcr test result (note this is the interval between symptom to results, in practice, if we assume sample collection takes 1-2 days, then with 2-3 days for results to come, it would be 4-5 days)

  # Parameters about tracing contects
  para$tracing_cc_onset <- 3 # set how many days we trace close contact back after symptom-based patient detection
  para$tracing_cc_ab <- para$delay_ab # set how many days we trace close contact back after a positive ab_test
  para$tracing_cc_pcr <- para$delay_pcr # set how many days we trace close contact back after a positive ab_test

  para$cc_success_symptom <- 0.85 # precentage close contact successfully traced after symptom-based patient detection
  para$cc_success_ab <- 0.75 # precentage close contact successfully traced after positive ab test
  para$cc_success_pcr <- 0.80 # precentage close contact successfully traced after positive pcr test

  para$qrate_symptom <- 0.5 # CC quarantine compliance rate based on symptom
  para$qrate_ab <- 0.7 # CC quarantine compliance rate based on positive RDT test
  para$qrate_pcr <- 0.7 # CC quarantine compliance rate based on positive pcr test

  # Parameters about testing tools
  para$ab_rate <- function(x, t_onset) 1/(1 + exp(7+t_onset-x)) # seroconversion rate from infection day, based on the clinical paper from Yumei Wen
  para$sensitivity_ab <- 0.9 # antibody test sensitivity
  para$sensitivity_pcr <- 0.999 # pcr test sensitivity
  para$samplefailure_pcr <- 0.3 # pcr sampling failure
  para$sensitivity_antig <- 0.8 # antigen sensitivity is 0.8


para <- init_para(para = para)
```


