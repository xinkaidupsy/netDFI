library(dplyr)
library(netDFI)

magna <- read.csv("/Users/xinkaidu/Documents/git/dfi_ggm/simulation/MAGNA_PTSD_weights.csv") %>%
  select(-cluster) %>%
  as.matrix

results_m <- dfi_ggm(magna, iter = 100, n = 250, ncores = 8, power = 0.95, n_misspec = 3, min_extra = 0.15)

net <- bootnet::genGGM(20, propPositive = 0.8, p = 0.4, graph = "random")
results_n <- dfi_ggm(net, iter = 100, n = 250, ncores = 8, power = 0.95, n_misspec = 5, min_extra = 0.15)

plot(results_n)[[5]]
qgraph::centrality(net, R2 = T)$R2
net2 <- bootnet::genGGM(20, propPositive = 0.8)
results_n2 <- dfi_ggm(net, iter = 100, n = 250, ncores = 8, power = 0.95, n_misspec = 5, min_extra = 0.15)

plot(results_n2)[[3]]

DASS <- read.csv("~/Documents/git/dfi_ggm/net/DASS_ggmModSelect_stepwise_weights.csv") %>%
  select(-cluster) %>%
  as.matrix

results_d <- dfi_ggm(DASS, iter = 100, n = 250, ncores = 8, power = 0.95, n_misspec = 5, min_extra = 0.15)

p <- plot(results_d)
p[[2]]
