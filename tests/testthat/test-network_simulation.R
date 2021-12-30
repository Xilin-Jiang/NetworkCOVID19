para <- init_para(1, 1, 1)

test_that("initialise default parameter setting for network generation", {
  para <- init_para(2, 1, 1)
  para$ego.pop_sz <- 200 # testing for a small sample size
  nw <- network_generate(para)
  expect_equal(nw[[4]] < 0.05, T)
})

test_that("loading default parameter setting for simulation", {
  load("Networks/example.network.Rdata")
  sim_rslt <- simulate_transmission(NW_SIM = NW_SIM)
  expect_equal(c(sum(sim_rslt$quarantine_daily), sum(sim_rslt$RDT_used), sum(sim_rslt$PCR_used)), c(0,0,0))
})


