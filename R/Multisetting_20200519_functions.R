###############################################
# Preparation functions
###############################################

# incubation time: shape = 7.92, scale = 0.69
incub <- function(x) floor(pmax(pmin(rgamma(x, shape = 7.92, scale = 0.69), 14),1)) # maximum incubation time is 14 days

# here we use the shape of weibull from the science paper.
# The onset should be the mode of weibull, so we compute the scale from the onset .
# weibull has the property of keeping the shape, and mode change with scale
# physical property of weibull also make sense: single event happen with an rate of t^k

g_transmissibility <- function(t_onset) para$R0_adj/para$num_cc * dweibull(1:(len_sim+1), shape = 2.826, scale = t_onset/0.8568)
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


#' initiate default parameters for network simulation
#'
#' @param setting integer. rural simulation; 2 for urban, 3 for slum
#' @param country_id 1 stands for Uganda demographic
#' @param social_distancing_flg 1 stands for no social distancing
#'
#' @return
#' @export
#'
#' @examples
init_nw_para <- function(setting, country_id, social_distancing_flg){
  # library(tidyverse)
  suppressMessages(library(ergm))
  suppressMessages(library(ergm.ego))

  para <- list()
  para$setting <- setting # 1=rural 2=affluent 3=slum
  # Environment setup
  para$pop_sz <- 1000 # 5000
  para$ego.pop_sz <- 200 # using 200 to fit the egocentric network
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


  if(country_id == 1){
    age_mix <- read.csv("Uganda_age_mixture_population_level.csv")
    para$age_mix <- age_mix$total_age
    para$Home_age_mix <- age_mix$home_age
  }else{
    if(country_id == 2){
      AGE_matrix <- read.csv("south_africa_age.csv", header = F)
      HOME_age_matrix <- read.csv("home_south_africa_age.csv", header = F)
    }else if(country_id == 3){
      AGE_matrix <- read.csv("kenya_age.csv", header = F) # kenya
      HOME_age_matrix <- read.csv("home_kenya_age.csv", header = F) # kenya
    }
    else if(country_id == 4){
      AGE_matrix <- read.csv("Nigeria_age.csv", header = F) # nigeria
      HOME_age_matrix <- read.csv("home_Nigeria_age.csv", header = F) # Nigeria
    }
    # the age mixing matrix need to collpase into 3*3 since we only have age distribution for 3 groups.
    ##  ##  ##  ##
    ## this step is important: the contact matrix is not symetrical: row and columns represent different things
    ##  ##  ##  ##
    # the columns represent contact, we take colSums
    # the rows are the participant, we take the average number of contact for the participant: rowMeans
    age <- cbind(rowMeans(AGE_matrix[,1:3]),rowMeans(AGE_matrix[,4:5]),rowMeans(AGE_matrix[,6:16]))
    AGE_matrix <- rbind(colSums(age[1:3,]),colSums(age[4:5,]),colSums(age[6:16,]))
    # weight each column by the age distribution: we have to change the matrix to reflect the age distribution in the population
    AGE_matrix <- (as.matrix(rep(1,length(para$age_dist))) %*% t(para$age_dist)) * AGE_matrix

    # matrix in the end should be the number of contact for each age group pair
    AGE_matrix <- (AGE_matrix + t(AGE_matrix))/2
    para$age_mix <- as.matrix(AGE_matrix)[which(upper.tri(AGE_matrix,diag = T))]

    age <- cbind(rowMeans(HOME_age_matrix[,1:3]),rowMeans(HOME_age_matrix[,4:5]),rowMeans(HOME_age_matrix[,6:16]))
    HOME_age_matrix <- rbind(colSums(age[1:3,]),colSums(age[4:5,]),colSums(age[6:16,]))
    HOME_age_matrix <- (as.matrix(rep(1,length(para$age_dist))) %*% t(para$age_dist)) * HOME_age_matrix
    # matrix in the end should be a triangle one
    HOME_age_matrix <- (HOME_age_matrix + t(HOME_age_matrix))/2
    para$Home_age_mix <- as.matrix(HOME_age_matrix)[which(upper.tri(HOME_age_matrix,diag = T))]
  }

  # adjust the age matrix to represent the specified household contact rate
  para$age_mix <- para$Home_age_mix + (para$age_mix - para$Home_age_mix) *
    sum(para$Home_age_mix)/sum(para$age_mix - para$Home_age_mix) *  (1-para$percent_HH_cc)/para$percent_HH_cc

  # adjust R0 for 1) young people susceptibility 2) subclinical cases
  # adjust for the social distancing
  para$num_cc_scdst <- para$num_cc * ((1-para$percent_HH_cc)*para$Non_HH_CC_rate + para$percent_HH_cc) # reduce the number of cc
  para$age_mix_scdst <- para$Home_age_mix + (para$age_mix - para$Home_age_mix) * para$Non_HH_CC_rate
  para$percent_HH_cc_scdst <- para$percent_HH_cc/((1-para$percent_HH_cc)*para$Non_HH_CC_rate + para$percent_HH_cc)

}
###############################################
# load an existing the contact network
###############################################
#' load a simulated network with specific parameter settings
#'
#' @param setting integer. setting == 1 refers to rural setting with low contact rate; setting == 2 refers to urban setting with normal contact rate; setting == 3 refers to slum setting with high contact rate;
#' @param v_id integer. choosing a parameter value to perturb.
#' @param s_id integer. choosing a value for the target parameter.
#' @param para list. parameters that control the simulation
#'
#' @return list. Fist element of the list is a network object; second is the loaded para list.
#' @export
#'
#' @examples load_sample_nw(1, 8, 3 para)
load_sample_nw <- function(setting,v_id,s_id, para){
  NW_SIM <- NA
  while (is.na(NW_SIM)){
    try({
      idx <- sample(1:20, 1)
      if(v_id == 5 | v_id == 8 | v_id ==  12 | v_id ==  14){
        load(paste0("Networks/network_",idx,"_set",setting,"v",v_id,"s",s_id,".Rdata"))
      }else{
        social_dst_idx <- 1+(1-para$Non_HH_CC_rate)/0.2
        load(paste0("Networks/network_",idx,"_set",setting,"v8s", social_dst_idx,".Rdata"))
      }
    })
  }
  nw <- NW_SIM[[1]]
  para_nw <- nw[[2]]
  para$family_lable <- para_nw$family_lable
  para$clustering_x <- para_nw$clustering_x
  para$clustering_y <- para_nw$clustering_y
  para$clustering_effect <- para_nw$clustering_effect

  # load all the network parameter, include AGE
  para$pop_sz <- para_nw$pop_sz
  para$AGE <- para_nw$AGE
  para$age_dist <- para_nw$age_dist
  para$family_sz <- para_nw$family_sz
  para$HH_dist <- para_nw$HH_dist
  para$num_cc_scdst <- para_nw$num_cc_scds
  para$percent_HH_cc_scdst <- para_nw$percent_HH_cc_scdst
  para$age_mix_scdst <- para_nw$age_mix_scdst
  # non-social-distanced para
  para$age_mix <- para_nw$age_mix
  para$Home_age_mix <- para_nw$Home_age_mix
  para$percent_HH_cc <- para_nw$percent_HH_cc
  para$num_cc <- para_nw$num_cc

  return(list(nw[[1]], para))
}
###############################################
# #  generate the contact network
###############################################
initiate_nw <- function(para){
  para$AGE <- unlist(sapply(1:length(para$age_dist), function(x) rep(x,round(para$age_dist[x] * para$pop_sz))))
  stopifnot(length(para$AGE) == para$pop_sz)
  nw <- network.initialize(n = para$pop_sz , directed = FALSE)
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
  nw <- set.vertex.attribute(nw, "family", family_lable)
  nw <- set.vertex.attribute(nw, "age", para$AGE)
  # nw <- set.vertex.attribute(nw, "clustering", c(clustering_x,clustering_y))
  nw <- set.vertex.attribute(nw, "clustering_x", clustering_x)
  nw <- set.vertex.attribute(nw, "clustering_y", clustering_y)
  return(list(nw, para))
}

network_generate <- function(para, searched_clustering_number){
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
  ego.net.fitting <- ergm.ego(est.ego ~ edges  + nodematch("family") + mm("age", levels2 = -6) + absdiff("clustering_x", pow=2) + absdiff("clustering_y", pow=2),
                              control = control.ergm.ego(ergm = control.ergm(MCMLE.maxit = 400, SAN.maxit = 200)))


  return(list(ego.net.fitting,para,searched_clustering_number, deviation_target_statistics))

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
  edge_lst <- as.edgelist(simulate(est_nw, popsize =population.ego.sim))
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
      transmissibility <- g_transmissibility(O[i]) * ifelse(Z[i] == 1, para$asym_rate, 1)
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
