
test_that("initialise default parameter setting for network generation,
          loading default parameter setting for simulation", {
  para <- init_para(2, 1, 1)
  para$ego.pop_sz <- 200 # testing for a small sample size
  para$death_rate <- c(1,1,1)
  para$death_delay <- 0
  nw <- network_generate(para)
  expect_equal(nw[[4]] < 0.05, T)
  sim_rslt <- simulate_transmission(NW_SIM = nw)
  # there should be no test kits used
  expect_equal(c(sum(sim_rslt$quarantine_daily), sum(sim_rslt$RDT_used), sum(sim_rslt$PCR_used)), c(0,0,0))
  # death should be equal to case, when death_rate = 1
  expect_equal(sum(sim_rslt$death_daily)+para$E_0, sum(sim_rslt$new_daily_case))
})

test_that("contact number setting", {
  para <- init_para(2, 1, social_distancing_flg = 4, contact_number = 30)
  expect_equal(c(para$num_cc, para$num_cc_scdst), c(30, 16.14))
})

test_that("using egocentric parameter -- can we make sure the resulting network has the save family label defined?", {
  para <- init_para(2, 1, 1)
  para$ego.pop_sz <- 200 # testing for a small sample size
  nw <- network_generate(para)
  a <- simulate(nw[[1]], popsize = para$pop_sz, nsim = 10)
  para <- nw[[2]]
  expect_equal(network::get.vertex.attribute(a[[1]], "family"), para$family_lable - 40*(0:4))
})

test_that("wrapper function", {
  results <- network_covid_simulate(rep_num = 1, network_num = 1, output = "example", para = NA)
  plt <- plot_epidemic_curves(results, title_fig = "")
})


