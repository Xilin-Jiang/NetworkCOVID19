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
A full explanation of parameter setting is in **Customise your simulation settings**.


```
para$death_rate <- c(1,1,1)
para$death_delay <- 0
nw <- network_generate(para)
expect_equal(nw[[4]] < 0.05, T)
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
