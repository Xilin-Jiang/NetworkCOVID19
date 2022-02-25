

test_that("number of contact scaling to large number", {
  para <- list()
  para$num_cc <- 20
  para$pop_sz <- 200
  para <- init_para(2, 1, 1, para = para)
  expect_equal(para$num_cc * para$cc_cyl, para$Tnum_cc)
  nw <- network_generate(para)

  NW_SIM <- nw
  len_sim <- 10
  PCR <- F
  RDT <- F

  ergm.fitting <- NW_SIM[[1]]
  para <- NW_SIM[[2]]
  para$len_sim <- len_sim

  para$onsetiso <- 0

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

    # add another round of infection to get the correct contact scaling
    if(!is.null(para$cc_cyl) & para$cc_cyl > 1){
      C_WD <- C_update(C, Q, ergm.fitting, para, WD_flg = T)
      lst <- I_O_update(I, Q, C_WD, O, Z, trace_inf_n, para, WD_flg = T)
      I <- lst[[1]]
      O <- lst[[2]]
      Z <- lst[[3]]
      trace_inf_n <- lst[[4]]
    }
    C <- C_update(C, Q, ergm.fitting, para)
    lst <- I_O_update(I, Q, C, O, Z, trace_inf_n, para)
    I <- lst[[1]]
    O <- lst[[2]]
    Z <- lst[[3]]
    trace_inf_n <- lst[[4]]

    C_eff <- C
    # combine the contact network from two rounds of simulation
    if(!is.null(para$cc_cyl) & para$cc_cyl > 1){
      C <- C_unify(C, C_WD)
    }
    #####################################
    # test the contact are merged properly
    #####################################
    expect_equal(sum(C == 12), union(which(C_eff == 12), which(C_WD == 12)) %>% length())

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
  # the cycles are correct (within day transmission would not cause problem in day records)
  expect(max(I,na.rm = T), para$len_sim)

})

