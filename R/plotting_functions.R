#' @import ggplot2
NULL

# use color blind friendly scheme
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# set the commonly used color
blue <- cbPalette[6]
red <- cbPalette[7]
green <- cbPalette[4]
orange <- cbPalette[2]
grey <- cbPalette[1]
yellow <- cbPalette[5]
purple <- cbPalette[8]
skyblue <- cbPalette[3]

#' Title
#'
#' @param results Output of the network_covid_simulate function
#' @param title_fig Give a title of the figure
#'
#' @return A ggplot figure object which plot the epidemic curves under four conditions.
#' @export
#'
#' @examples  results <- network_covid_simulate(rep_num = 1, network_num = 1, output = "example", para = NA)
#' plt <- plot_epidemic_curves(results, title_fig = "")
#'
plot_epidemic_curves <- function(results, title_fig = ""){

  len_sim <- dim(results$new_daily_case[[1]])[2]

  err_baseline <- apply(results$new_daily_case[[1]], 2, function(x) sd(x)/sqrt(length(x)))
  err_RDT <- apply(results$new_daily_case[[2]], 2, function(x) sd(x)/sqrt(length(x)))
  err_PCR <- apply(results$new_daily_case[[3]], 2, function(x) sd(x)/sqrt(length(x)))
  err_PCR_RDT <- apply(results$new_daily_case[[4]], 2, function(x) sd(x)/sqrt(length(x)))

  df_plt <- data.frame(day = 1:len_sim, RDT_new_daily = colMeans(results$new_daily_case[[2]]), new_daily = colMeans(results$new_daily_case[[1]]),
                       PCR_new_daily = colMeans(results$new_daily_case[[3]]), RDT_PCR_new_daily = colMeans(results$new_daily_case[[4]]),
                       err_baseline, err_RDT, err_PCR,err_PCR_RDT)
  plt <- ggplot(df_plt,aes(x=day)) +
    geom_line(aes( y = RDT_new_daily, color = "Strategy 1: RDT")) +
    geom_line(aes(y = new_daily, color = "Baseline")) +
    geom_line(aes(y = PCR_new_daily, color = "Strategy 2: PCR")) +
    geom_line(aes(y = RDT_PCR_new_daily, color = "Strategy 3: RDT & PCR")) +
    geom_ribbon(aes( ymax = RDT_new_daily + 2*err_RDT, ymin = RDT_new_daily - 2*err_RDT), fill = blue,alpha = 0.3) +
    geom_ribbon(aes( ymax = new_daily + 2*err_baseline, ymin = new_daily - 2*err_baseline), fill = red,alpha = 0.3) +
    geom_ribbon(aes( ymax = RDT_PCR_new_daily + 2*err_PCR_RDT, ymin = RDT_PCR_new_daily - 2*err_PCR_RDT), fill = purple,alpha = 0.3) +
    geom_ribbon(aes( ymax = PCR_new_daily + 2*err_PCR, ymin = PCR_new_daily- 2*err_PCR), fill = green,alpha = 0.3) +
    labs(x= "Days from first importing case", y = "New case per day",title = title_fig)+
    scale_color_manual(values = c("Baseline" = red, "Strategy 1: RDT" = blue,
                                  "Strategy 2: PCR" = green, "Strategy 3: RDT & PCR" = purple)) +
    theme(legend.position = "right", panel.background=element_blank())

  return(plt)
}
