###############################################
# Preparation functions
###############################################

# incubation time: shape = 7.92, scale = 0.69
incub <- function(x) floor(pmax(pmin(rgamma(x, shape = 7.92, scale = 0.69), 14),1)) # maximum incubation time is 14 days

# here we use the shape of weibull from the science paper.
# The onset should be the mode of weibull, so we compute the scale from the onset .
# weibull has the property of keeping the shape, and mode change with scale
# physical property of weibull also make sense: single event happen with an rate of t^k

g_transmissibility <- function(t_onset, para) para$R0_adj/para$num_cc * dweibull(1:(para$len_sim+1), shape = 2.826, scale = t_onset/0.8568)
# plot(g_transmissibility(20))
# weibull distribution: shape: 2.8 scale 5.7
# mean 5.665*(1-1/2.826)^(1/2.826)
# normalized transmissibility per contact para$R0_adj/para$num_cc

###############################################
# external functions (all function that could actually be used by other users)
###############################################
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @import stats
NULL

#' @import ergm
NULL

#' @import ergm.ego
NULL

#' initiate default parameters for simulation
#'
#' @param setting integer. rural simulation; 2 for urban, 3 for slum
#' @param country_id 1 stands for Uganda demographic
#' @param social_distancing_flg 1 stands for no social distancing
#'
#' @return a list contains all parameters at default setting for each scenatrio
#' @export
#'
#' @examples # initialise a normal para <- init_para()
init_para <- function(setting = 2, country_id = 1, social_distancing_flg = 1){
  para <- list()

  ###############################
  # parameters relevant to network
  ###############################
  para$setting <- setting # 1=rural 2=affluent 3=slum
  # Environment setup
  para$pop_sz <- 1000 # 5000
  para$ego.pop_sz <- 1000 # using 200 to fit the egocentric network
  para$Non_HH_CC_rate <- c(1,.8,.6,.4,.2)[social_distancing_flg]
  community_setting <- c("Rural", "Non-slum urban", "Slum")[setting]
  print(paste0("Simulate ",para$pop_sz," individuals in ",community_setting," setting"))
  ##########################################
  # loading country specific variables:
  # age distribution & household
  ##########################################
  if(country_id == 1){
    para$age_dist <- c(0.481, 0.203, 0.316) # Uganda
    para$HH_dist <- c(.11, .22, .27, .4) # Uganda
  }else if(country_id == 2){
    para$age_dist <- c(0.292, 0.193, 0.515) # South africa
    para$HH_dist <- c(.27, .35, .23, .15) # South Africa
  }else if(country_id == 3){
    para$age_dist <- c(0.419, 0.195, 0.386) # kenya
    para$HH_dist <- c(.19, .28, .3, .23) # kenya
  }else if(country_id == 4){
    para$age_dist <- c(0.440, 0.190, 0.370) # nigeria
    para$HH_dist <- c(.16, .26, .26, .32) # nigeria
  }

  ##########################################
  # processing demographic information;
  # for details: refer to the supplementary methods
  ##########################################
  # para$HH_affluent_dist <- c(.31, .5, .18, .02) # UK

  if (para$setting==1) {
    para$num_cc <- 7 # set daily close contact to be 7
    para$family_sz <- 5 # average household size 5
    # set the percentage of HH_cc
    para$percent_HH_cc <- .5
  }else if (para$setting==2) {
    para$num_cc <- 13
    para$family_sz <- 5
    para$percent_HH_cc <- .23
  }else if (para$setting==3) {
    para$num_cc <- 14
    para$family_sz <- 15
    para$percent_HH_cc <- .5
    para$HH_dist <- c(.00, .06, .17, .77) # afganistan
  }else print ("Parameter setting error")
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
  # Transmission parameter & tool availability
  if (para$setting==1) {
    para$symptom_report <- 0.7 # precentage infected report symptom
    para$R0 <- R0_baseline * para$num_cc # R naught
    para$theta <- 0.1 # quarantine effect: precentage of remaining transmission quarantined individual have
    para$pcr_available <- 1000*para$pop_sz/1000 # -1 # daily maximum PCR tests per 1000 population.
    para$ppe_coef <- 1 # if people wear ppe, then their transmissibility will be different outside their family
  }else if (para$setting==2) {
    para$symptom_report <- 0.7
    para$R0 <- R0_baseline * para$num_cc
    para$theta <- 0.1
    para$pcr_available <-1000*para$pop_sz/1000
    para$ppe_coef <- 1
  }else if (para$setting==3) {
    para$symptom_report <- 0.6
    para$R0 <- R0_baseline * para$num_cc
    para$theta <- 0.2
    para$pcr_available <- 1000*para$pop_sz/1000 # 2*para$pop_sz/1000
    para$ppe_coef <- 1
  }else print ("Parameter setting error")

  # subclinical rate
  para$sub_clini_rate <- 0.3
  para$asym_rate <- 0.2 # transmission rate of asymptomatic patient
  para$death_rate <- c(0.002,0.002,0.01)
  para$death_delay <- 10

  # Parameters about Infected pts
  para$ab_test_rate <- 0.7 # % accepted ab testing among detected (symptomatic) patients
  para$pcr_test_rate <- 0.8 # % accepted pcr testing among detected (symptomatic) patients

  para$onsetiso <- 0.2 # Isolation compliance rate at onset based on symptom
  para$abiso <-0.9 # Isolation compliance rate among ab testing positive
  para$pcriso <-0.9 # Isolation compliance rate among pcr testing positive

  para$delay_symptom <- 1 # days of delay after onset to detect symptomatic patients
  para$delay_ab <- 8 # days of delay after onset to receive ab test and obtain result
  para$delay_pcr <- 5 # days of delay after onset to report pcr test result

  # Parameters about tracing contects
  para$tracing_cc_onset <- 3 # set how many days we trace close contact back after symptom-based patient detection
  para$tracing_cc_ab <- para$delay_ab # set how many days we trace close contact back after a positive ab_test
  para$tracing_cc_pcr <- para$delay_pcr # set how many days we trace close contact back after a positive ab_test

  para$cc_success_symptom <- 0.85 # precentage close contact successfully traced after symptom-based patient detection
  para$cc_success_ab <- 0.75 # precentage close contact successfully traced after positive ab test
  para$cc_success_pcr <- 0.80 # precentage close contact successfully traced after positive pcr test

  para$qrate_symptom <- 0.5 # CC quarantine compliance rate based on symptom
  para$qrate_ab <- 0.7 # CC quarantine compliance rate based on positive ab test
  para$qrate_pcr <- 0.7 # CC quarantine compliance rate based on positive pcr test

  # Parameters about testing tools
  para$ab_rate <- function(x, t_onset) 1/(1 + exp(7+t_onset-x)) # seroconversion rate from infection day, based on the clinical paper from Yumei Wen
  para$sensitivity_ab <- 0.9 # ab test sensitivity
  para$sensitivity_pcr <- 0.999 # pcr test sensitivity
  para$samplefailure_pcr <- 0.3 # pcr sampling failure
  para$sensitivity_antig <- 0.8 # antigen sensitivity is 0.8

  # adjust R0 for next generation matrix
  para <- R0_adjust(para)

  return(para)
}

###############################################
# generate the contact network
###############################################
#' Simulate synthetic population estimate ERGM contact networks.
#'
#' @param para list contains all parameter, which should be from the function init_para; the paramter settings could be modified at this stage
#' @param output character string for the function to save the output into. The function will create and save results into "Networks" directory.
#' @param searched_clustering_number a parameter determin the strength of clustering within the population. Should use default setting unless
#' for very specific testing session (e.g. simulating lockdown scenario where there is no contact except designated shopping time)
#'
#' @return a list of two object: first is an ERGM object which could be used to simulate contact networks; second is a para object which was updated
#' with the synthetic population age and family information.
#' @export
#'
#' @examples
network_generate <- function(para, output = "example", searched_clustering_number = 4){
  para_ego <- para
  para_ego$pop_sz <- para$ego.pop_sz
  nw.ego <- initiate_nw(para_ego)[[1]]
  # save the network properties, including AGE, family_lable, geographical location
  para <- initiate_nw(para)[[2]]

  est1 <- NA
  grid_id <- 1
  #################################################################################
  # step 1: fit a small ERGM network to get the coefficients
  #################################################################################
  # perform grid search to get the local clustering coefficient
  while(is.na(est1)){
    clustering_effect <- (para_ego$pop_sz/2*para_ego$num_cc_scdst)*(1 - para_ego$percent_HH_cc_scdst) * searched_clustering_number #  14 (25) 23 (20)
    target.stats <- c(para_ego$pop_sz/2 * para_ego$num_cc_scdst, para_ego$pop_sz/2 * para_ego$num_cc_scdst * para_ego$percent_HH_cc_scdst,
                      (para_ego$age_mix_scdst/sum(para_ego$age_mix_scdst) * para_ego$pop_sz/2 * para_ego$num_cc_scdst)[1:5], clustering_effect, clustering_effect)

    try({
      suppressMessages(  est1 <- ergm(nw.ego ~ edges  + nodematch("family") + mm("age", levels2 = -6) + absdiff("clustering_x", pow=2) + absdiff("clustering_y", pow=2),
                                      target.stats = target.stats,control = control.ergm(MCMLE.maxit = 400, SAN.maxit = 200)) )
      ego.sim100 <- simulate(est1,nsim = 100)
      sim.stats <- attr(ego.sim100,"stats")
      trgt <- rbind(colMeans(sim.stats), est1$target.stats)
      deviation_target_statistics <- mean(abs(trgt[1,] - trgt[2,])/trgt[2,])
      if(deviation_target_statistics > 0.05){
        est1 <- NA
      }
    })

    searched_clustering_number <- searched_clustering_number + 1 # reduce the searched number if it did not converge; 6 is when there is no geographical clustering effect
    print(paste0("search ", grid_id, " for geographical clustering"))
    grid_id <- grid_id + 1
  }
  para$clustering_effect <- clustering_effect

  ########################################################
  # Step 2: Transform the parameters to an egocentric model
  ########################################################
  est.ego <-  as.egodata(simulate(est1))
  suppressMessages(
    ego.net.fitting <- ergm.ego(est.ego ~ edges  + nodematch("family") + mm("age", levels2 = -6) + absdiff("clustering_x", pow=2) + absdiff("clustering_y", pow=2),
                                control = control.ergm.ego(ergm = control.ergm(MCMLE.maxit = 400, SAN.maxit = 200)))
  )
  # save the data
  dir.create("Networks")
  NW_SIM <- list(ego.net.fitting,para,searched_clustering_number, deviation_target_statistics)
  save(NW_SIM, file = paste0("Networks/",output,".network.Rdata"))

  return(NW_SIM)
}


initiate_nw <- function(para){
  para$AGE <- unlist(sapply(1:length(para$age_dist), function(x) rep(x,round(para$age_dist[x] * para$pop_sz))))
  stopifnot(length(para$AGE) == para$pop_sz)
  nw <- network::network.initialize(n = para$pop_sz , directed = FALSE)
  # specify house hold
  adult_idx <- which(para$AGE == 3)
  family_core <- sample(adult_idx, para$pop_sz/para$family_sz, replace =F)# a family should have at least one adult
  # assign other adults randomly
  non_core_adult_idx <- adult_idx[! adult_idx %in% family_core]
  # family_adults <- split(non_core_adult_idx, sample(para$pop_sz/para$family_sz, length(non_core_adult_idx), replace = TRUE) )

  non_core_idx <- sample(c(which(para$AGE <= 2), non_core_adult_idx))
  # comupte the split point for family size
  family_sz_label <- cumsum(para$HH_dist)

  # assigning family labels
  family_lable <- rep(NA, para$pop_sz)
  family_lable[family_core] <- 1:(para$pop_sz/para$family_sz)

  # using a rolling rotation to create family based on the household size distribution
  for(sz_gp in 1:2){
    # for each group 2-3/4-5, I assign half of family group 2-3 to be 2 and the other half to be 3
    idx1 <- floor(family_sz_label[sz_gp] * para$pop_sz/para$family_sz):(para$pop_sz/para$family_sz)
    family_lable[non_core_idx[1:(idx1[length(idx1)] - idx1[1] + 1)]] <-  idx1
    non_core_idx <- non_core_idx[(1+length(idx1)):length(non_core_idx)]
    # pick half of thsese families to assiagn another member: the result will be in the family of size (2-3, that's the data we have unfortunately),
    # half of then are havign 2 members and the other half have 3 members
    idx2 <- floor((family_sz_label[sz_gp] + family_sz_label[sz_gp+1])/2 * para$pop_sz/para$family_sz):(para$pop_sz/para$family_sz)
    family_lable[non_core_idx[1:(idx2[length(idx2)] - idx2[1] + 1)]] <-  idx2
    non_core_idx <- non_core_idx[(1+length(idx2)):length(non_core_idx)]
  }
  # assign the remaining people randomly to those family with more than 5 people
  family_lable[non_core_idx] <- sample(floor(family_sz_label[3] * para$pop_sz/para$family_sz + 1):(para$pop_sz/para$family_sz),length(non_core_idx), replace = T)

  # use a spiral function to create the location grid: each household will be roughly equal distance
  phi <- c(1)
  for(i in 2:(para$pop_sz/para$family_sz)){
    phi[i] <- phi[i-1] +1/phi[i-1]
  }
  clustering_x <- phi * cos(phi *pi) + rnorm(length(phi))
  clustering_x <- clustering_x[family_lable]
  clustering_y <- phi * sin(phi *pi) + rnorm(length(phi))
  clustering_y <- clustering_y[family_lable]

  para$family_lable <- family_lable
  para$clustering_x <- clustering_x
  para$clustering_y <- clustering_y


  # this family and age label will be add to the population
  nw <- network::set.vertex.attribute(nw, "family", family_lable)
  nw <- network::set.vertex.attribute(nw, "age", para$AGE)
  # nw <- network::set.vertex.attribute(nw, "clustering", c(clustering_x,clustering_y))
  nw <- network::set.vertex.attribute(nw, "clustering_x", clustering_x)
  nw <- network::set.vertex.attribute(nw, "clustering_y", clustering_y)
  return(list(nw, para))
}

R0_adjust <- function(para){
  #########################################################################
  # adjust R0 using next generation matrix for 1) young people susceptibility 2) subclinical cases
  #########################################################################
  # ajust for R0
  # norm1 <- (sum(para$age_mix) - para$age_mix[c(1)]/2- sum(para$age_mix[c(2,4)])/4)/sum(para$age_mix)
  # using next generation matrix to compute norm1; only kept terms connected with young people susceptibility
  Cyy <- para$age_mix[c(1)] # number of comtact for young <--> young
  Coy <- sum(para$age_mix[c(2,4)]) # number of comtact for young <--> old
  Coo <- sum(para$age_mix[c(3,5,6)]) # number of comtact for old <--> old
  Ny <- para$age_dist[1] # number of young people, Cyy/Ny is the average number of young contact for a yong person
  No <- sum(para$age_dist[c(2,3)]) # number of young people, Cyy/Ny is the average number of young contact for a yong person
  y_sus <- 0.5 # susceptability of young person
  NGM <- matrix(c(y_sus * Cyy/Ny,  Coy/Ny, y_sus * Coy/No, Coo/No) , nrow = 2)
  trNGM <- sum(diag(NGM))
  detNGM <- det(NGM)
  Spectral_radius_half <- trNGM + (trNGM^2 - 4*detNGM )^0.5

  y_sus <- 1 # susceptability of young person
  NGM <- matrix(c(y_sus * Cyy/Ny,  Coy/Ny, y_sus * Coy/No, Coo/No) , nrow = 2)
  trNGM <- sum(diag(NGM))
  detNGM <- det(NGM)
  Spectral_radius_1 <- trNGM + (trNGM^2 - 4*detNGM )^0.5

  norm1 <- Spectral_radius_half/Spectral_radius_1
  # norm2 account for the redution of the low transmission rate of asymptomatic cases
  norm2 <- 1 - para$sub_clini_rate * (1 - para$asym_rate)
  para$R0_adj <- para$R0/(norm1 * norm2)
  ##################################
  # compute Re
  norm_age_mix_scdst <- para$age_mix_scdst/sum(para$age_mix_scdst)
  para$Re <- para$R0_adj/para$num_cc *
    ((1- norm_age_mix_scdst[c(1)]/2- sum(norm_age_mix_scdst[c(2,4)])/4) * para$num_cc_scdst) *
    norm2 # approximate Re (haven't taken age-susceptibility into account)
  return(para)
}


###############################################
# wrapper for transmission simulation
###############################################
#' Title
#'
#' @param NW_SIM A list of two object, which is the output of network_generate function.
#' @param input_location A file location to read the NW_SIM object.
#' @param PCR logical True means PCR test is deployed for the containment
#' @param RDT logical True means rapid diagnostic tests are deployed for the containment
#' @param len_sim numeric number of days to simulate for the outbreak, default 100.
#'
#' @return a dataframe contains the simulated results.
#' @export
#'
#' @examples
simulate_transmission <- function(NW_SIM = NA, input_location = "Networks/example.network.Rdata", PCR = F, RDT = F, len_sim = 100){
  if(is.na(NW_SIM)){
    if(stringr::str_detect(input_location, ".network.Rdata")){
      load(input_location)
    }else{
      load(paste0(input_location, ".network.Rdata"))
    }
  }
  ego.net.fitting <- NW_SIM[[1]]
  para <- NW_SIM[[2]]
  para$len_sim <- len_sim

  if(RDT){
    para$ab_test_rate  <- para$pcr_test_rate # consent for antigen
    para$ab_rate <- function(x, t_onset) (1-para$samplefailure_pcr) # failure rate of sample for antigen
    # para$abiso # same as ab
    para$sensitivity_ab <- para$sensitivity_antig # sensitivity is 0.8
    para$tracing_cc_ab <- para$tracing_cc_onset
    para$delay_ab <- para$delay_symptom # we do the test as soon as the symptom is discovered
  }else if(!RDT && !PCR){
    para$onsetiso <- 0
  }

  # C is the contact matrix, which is traced for 7 days. If i,j element is 7, it means the latest contact of i,j is today,
  # if it is 0, then there is no contact between i,j in the past 7 days,
  C <- matrix(0, para$pop_sz,para$pop_sz)

  # I is a vector indicating when an individual has been infected
  # NA means susceptible, 0 means infected date
  I <- matrix(NA, para$pop_sz,1)
  trace_inf_n <- matrix(0, para$pop_sz,1)
  # Z is a vector indicating if the infected is detectable
  Z <- matrix(F, para$pop_sz,1)

  # onset just save the incubation time
  O <- matrix(NA, para$pop_sz,1)

  # Q is the vector indicating when an individual has been quarantine, NA means not quarantine, 1 means first day
  Q <- matrix(NA, para$pop_sz,1)
  RDT_n <- 0 # total number of RDT used
  PCR_n <- 0 # total number of PCR used
  rdt <- matrix(0, len_sim)
  pcr <- matrix(0, len_sim)

  C_lst <- list()
  O_lst <- list()
  I_lst <- list()
  Q_lst <- list()

  # initial case
  # para$E_0 <- floor(runif(1,min = 1, max = 3))
  init_idx <- sample(1:para$pop_sz, para$E_0)
  I[init_idx] <- 0
  O[init_idx] <- incub(para$E_0)
  Z[init_idx] <- F
  for(t in 1:len_sim){
    C_lst[[t]] <- C
    C <- C_update(C, Q, ego.net.fitting, para)
    lst <- I_O_update(I, Q, C, O, Z, trace_inf_n, para)
    I <- lst[[1]]
    O <- lst[[2]]
    Z <- lst[[3]]
    trace_inf_n <- lst[[4]]
    I_lst[[t]] <- I
    O_lst[[t]] <- O
    lst1 <- Q_update(I, Q, C, O, Z, RDT_n, PCR_n, para, flg_ab = RDT, flg_multi_ab = FALSE, flg_PCR=PCR, PCR_need_yesterday)
    Q <- lst1[[1]]
    RDT_n <-lst1[[2]]
    PCR_n <-lst1[[3]]
    PCR_need_yesterday <- lst1[[4]]
    Q_lst[[t]] <- Q
    rdt[t] <- RDT_n
    pcr[t] <- PCR_n
  }
  Rt <- rep(NA, len_sim)
  death <- rep(0, len_sim)
  for(t in 1:len_sim){
    i <- I_lst[[t]]
    idx <- which(i == 2) # from the second day they are infectious
    if(length(idx)) Rt[t] <- mean(trace_inf_n[idx])
    # compute death
    death_case <- which(i == para$death_delay)
    if(length(death_case)) death[t] <- sum(runif(length(death_case)) < para$death_rate[para$AGE[death_case]])
  }
  # plot the daily new case
  ab_new_daily <- rep(NA, len_sim)
  cap <- rep(NA, len_sim) # store the Q numbers
  ab_new_daily[1] <- para$E_0
  cap[1] <- sum(!is.na(Q_lst[[1]]) & Q_lst[[t]] < 14)
  for(t in 2:len_sim){
    ab_new_daily[t] <- sum(is.na(I_lst[[t-1]]))  - sum(is.na(I_lst[[t]]))
    cap[t] <- sum(!is.na(Q_lst[[t]]) & Q_lst[[t]] < 14)
  }
  return(data.frame(new_daily_case = ab_new_daily, quarantine_daily = cap, Re_daily = Rt, RDT_used = rdt, PCR_used = pcr, death_daily = death))

}

###############################################
# Transmission update functions
###############################################
# update C
# the function below provide a way to update daily contact
C_update <- function(C, Q, est_nw, para){
  # one day pass; C is 12 when just infected, then C decrease by 1 everyday pass
  C <- pmax(C-1, 0)
  # then update close contact of this day. we assume the close contact to happen between the people around
  ##########################
  # simulated an actual net from the egocentric estimation
  population.ego.sim <- data.frame(family = para$family_lable, age = para$AGE, clustering_x = para$clustering_x, clustering_y = para$clustering_y)
  edge_lst <- network::as.edgelist(simulate(est_nw, popsize =population.ego.sim))
  # cc store the distance from the close contact to index case for each close contact of a particular day

  prob_cc <- ifelse(!is.na(Q[edge_lst[,1]]), para$theta, 1) * ifelse(!is.na(Q[edge_lst[,2]]),para$theta, 1)
  true_cc <- (runif(dim(edge_lst)[1]) < prob_cc)
  C[cbind(edge_lst[true_cc,1], edge_lst[true_cc,2])] <- 12
  return(C)
}
# update the Infected individual & keep an record of the incubation time
I_O_update <- function(I, Q, C, O, Z, trace_inf_n, para){
  I <- I + 1
  for(i in 1:para$pop_sz){
    # i is the infection source, idx is the infected
    # cc_idx find the infected individual and he has contact someone that day
    # (our model is almost certain to have contact, unless he is quarantined)
    # C[-i,i] "-i" is to avoid count the index case himself
    if( !is.na(I[i]) & (length(which(C[-i,i] == 12)) | length(which(C[i,-i] == 12)) ) ){
      cc_idx <- union(which(C[-i,i] == 12), which(C[i,-i] == 12))
      u <- runif(length(cc_idx))
      # get the individual transmissibility, it need to be normalized by the number of individuals
      transmissibility <- g_transmissibility(O[i], para) * ifelse(Z[i] == 1, para$asym_rate, 1)
      # pick the index of individual 1) who is not yet infected & 2) transmited
      # the expected transimisstion rate need to be normalized by daily contact 2*num_cc
      # reduce the sesceptability for young people (AGE == 1) by half
      susceptibility <- ifelse(para$AGE[cc_idx] == 1, .5, 1)
      u1 <- runif(length(cc_idx))
      idx <- cc_idx[(u < transmissibility[I[i]]) & is.na(I[cc_idx]) & u1 < susceptibility]
      I[idx] <- 0 # consider the incubation time
      u1 <- runif(length(idx))
      # if Z=2, the infected case is detectable, if Z=1, subclinical, if Z=F, not reported and infected
      Z[idx] <- sapply(u1, function(x) if(x<para$symptom_report*(1-para$sub_clini_rate)) 2 else if(x<(1-para$sub_clini_rate)) F else 1)
      O[idx] <- incub(length(idx)) # incubation duration
      # trace the Rt
      trace_inf_n[i] <- trace_inf_n[i] + length(idx)
    }
  }
  return(list(I, O, Z,trace_inf_n))
}
########################################################
# Quarantine functions
########################################################
# this function will quarantine after ab tests
# flg_ab <- FALSE
# flg_multi_ab <- TRUE
# flg_PCR <- TRUE
########################################
Q_update <- function(I, Q, C, O, Z, RDT_n, PCR_n, para, flg_ab, flg_multi_ab, flg_PCR, PCR_need_yesterday){
  Q <- Q + 1
  PCR_need_today <-0
  PCR_n_today <-0
  if (sum(!is.na(I))>para$infect_sz) { # take containment intervention only when the infected number passes the threshold
    # quarantine close contact
    for(i in 1:para$pop_sz){
      if(!is.na(I[i]) & (Z[i] == 2)){ # only when Z == 2 the case is reported at all
        u <- runif(4)

        # case detection strategy
        ok_onset_iso <- (I[i]-O[i]) == para$delay_symptom &
          u[1] < para$onsetiso
        ok_pcr_test <- flg_PCR &
          (PCR_n_today <= para$pcr_available) &
          ((I[i]-O[i]) == para$delay_pcr) &
          u[2] < para$pcr_test_rate
        ok_ab_test <-  flg_ab &
          ((I[i]-O[i]) == para$delay_ab) &
          u[3] < para$ab_test_rate
        ok_multiab_test <- flg_multi_ab &
          ((I[i]-O[i]) == 6 | (I[i]-O[i]) == 8 | (I[i]-O[i]) == 10) &
          u[4] < para$ab_test_rate

        # quarantine strategy
        # 1. test with RDT
        if(ok_ab_test & is.na(Q[i])){
          RDT_n <- RDT_n + 1
          u <- runif(2)
          seroconvert <- para$ab_rate(I[i],O[i]) # if the test shows positive
          if(u[1] < para$abiso & u[2] < para$sensitivity_ab*seroconvert){ # if the infected is isolated, quarantine their close contact
            Q[i] <- 1
            cc <- union(which(C[-i,i] >= 12 - para$tracing_cc_ab), which(C[i,-i] >= 12 - para$tracing_cc_ab))
            u <- runif(length(cc))
            # make sure the cc_success_rate for family is always 1
            tracing_rate <- ifelse(para$family_lable[i] == para$family_lable[cc], 1, para$cc_success_ab)
            Q[cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_ab]] <- 1
          }
        } else if( ok_multiab_test & is.na(Q[i])){
          RDT_n <- RDT_n + 1
          u <- runif(2)
          seroconvert <- para$ab_rate(I[i],O[i]) # if the test shows positive
          if(u[1] < para$abiso & u[2] < para$sensitivity_ab*seroconvert){ # if the infected is isolated, quarantine their close contact
            Q[i] <- 1
            cc <- union(which(C[-i,i] >= 12 - (I[i]-O[i])), which(C[i,-i] >= 12 - (I[i]-O[i])))
            u <- runif(length(cc))
            tracing_rate <- ifelse(para$family_lable[i] == para$family_lable[cc], 1, para$cc_success_ab)
            Q[cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_ab]] <- 1
          }
        }

        # 2. encourage non-positive to quarantine
        if(ok_onset_iso & is.na(Q[i])){
          Q[i] <- 1
          cc <- union(which(C[-i,i] >= 12 - para$tracing_cc_onset), which(C[i,-i] >= 12- para$tracing_cc_onset))
          u <- runif(length(cc))
          # make sure the cc_success_rate for family is always 1
          stopifnot(length(para$family_lable[cc]) == length(cc))
          tracing_rate <- ifelse(para$family_lable[i] == para$family_lable[cc], 1, para$cc_success_symptom)
          #####################################
          # for debug
          # print(length(cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_symptom]))
          #####################################
          Q[cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_symptom]] <- 1
        }

        # perform PCR on the rest
        if(ok_pcr_test & is.na(Q[i])){
          PCR_need_today <- PCR_need_today + 1
          pcravailable <- ifelse(para$setting==3,
                                 (runif(1)<(para$pcr_available/(PCR_need_yesterday+1)*2)),
                                 T) # in a slum, when PCR is not sufficient, smaller id people have a higher chance of geting the PCR test
          if (pcravailable) {
            PCR_n <- PCR_n + 1
            PCR_n_today <- PCR_n_today +1
            u <- runif(2)
            if(u[1] < para$pcriso & u[2] < para$sensitivity_pcr*(1-para$samplefailure_pcr)){ # if the infected is isolated, quarantine their close contact
              Q[i] <- 1
              cc <- union(which(C[-i,i] >= 12 - para$tracing_cc_pcr), which(C[i,-i] >= 12 - para$tracing_cc_pcr))
              u <- runif(length(cc))
              # make sure the cc_success_rate for family is always 1
              tracing_rate <- ifelse(para$family_lable[i] == para$family_lable[cc], 1, para$cc_success_pcr)
              Q[cc[is.na(Q[cc]) & u < tracing_rate * para$qrate_pcr]] <- 1

            }
          }
        }

      }
    }
    # release those after 14 days of quarantine
    for(i in 1:para$pop_sz){
      # we only release those who are do not have onset: if they do have onset,
      # they will be keep quarantined (this is our strategy to calculate medical needs)
      if(!is.na(Q[i]) & Q[i] >= 14){
        if(is.na(I[i])){
          Q[i] <- NA
        }else if(I[i] < O[i]){ # so this guys haven't develop any symptom, and will be free
          Q[i] <- NA
        }
      }
    }
  }
  return(list(Q,RDT_n,PCR_n,PCR_need_today))
}
