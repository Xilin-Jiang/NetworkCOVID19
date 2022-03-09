
test_that("contact number setting", {
  para <- list()
  para$num_cc <- 30
  para <- init_para(2, 1, social_distancing_flg = 4, para = para)
  expect_equal(c(para$Tnum_cc, para$num_cc_scdst), c(30, 8.07))
})

test_that("two methods of fitting network", {
  # fitting a small ERGM and construct large network
  para <- list()
  para$community.pop_sz <- 200
  para$death_rate <- c(1,1,1)
  para$death_delay <- 0
  para <- init_para(2, 1, 1, para = para)
  # testing for joint simulation scheme
  nw <- network_generate(para)
  expect_equal(nw[[4]] < 0.05, T)
  expect_equal(network::network.size(aggregate_simulation(nw[[1]], para)), para$pop_sz)

  sim_rslt <- simulate_transmission(NW_SIM = nw)
  # there should be no test kits used
  expect_equal(c(sum(sim_rslt$quarantine_daily), sum(sim_rslt$RDT_used), sum(sim_rslt$PCR_used)), c(0,0,0))
  # death should be equal to case, when death_rate = 1
  expect_equal(sum(sim_rslt$death_daily)+para$E_0, sum(sim_rslt$new_daily_case))
  expect_equal(sum(sim_rslt$death_daily)+para$E_0, sum(sim_rslt$new_daily_case))
  # fitting a ERGM same size as target population
  para <- list()
  para$community.pop_sz <- 200
  para$pop_sz <- 200
  para <- init_para(2, 1, 1, para = para)
  nw <- network_generate(para)
  expect_equal(class(nw[[1]]), "ergm")
  para <- nw[[2]]
  expect_equal(length(para$family_lable), para$pop_sz)
  expect_equal(length(para$AGE), para$pop_sz)
  expect_equal(length(para$clustering_x), para$pop_sz)
  expect_equal(length(para$clustering_y), para$pop_sz)
})

test_that("population size setting", {
  # make sure there is no parameter setting error for population sampling
  for(pop_sz in c(100,150,250,300,350,450,500)){
    para <- list()
    para$community.pop_sz <- 100
    para$pop_sz <- pop_sz
    para <- init_para(2, 1, 1, para = para)
    para <- initiate_nw(para)[[2]]
    expect_equal(length(para$family_lable), para$pop_sz)
    expect_equal(length(para$AGE), para$pop_sz)
    expect_equal(length(para$clustering_x), para$pop_sz)
    expect_equal(length(para$clustering_y), para$pop_sz)
  }
})

test_that("wrapper function test no quarantine clost contact", {
  para <- list()
  para$community.pop_sz <- 100
  para$pop_sz <- 300
  para <- init_para(2, 1, 1, para= para)
  # testing for joint simulation scheme
  results <- network_covid_simulate(rep_num = 1, network_num = 1, output = "example", para = para)
  plt <- plot_epidemic_curves(results, title_fig = "")
})


