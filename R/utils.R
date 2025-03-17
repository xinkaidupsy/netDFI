library(psychonetrics)
library(dplyr)

# Load the data:
data("StarWars")

# Observed variables:
obsvars <- c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10")

# Form GGM model:
ggm1 <- ggm(StarWars, vars = obsvars)

# Run model:
ggm1 <- ggm1 %>% runmodel

ggm1 <- ggm1 %>%
  prune(adjust = "fdr", alpha = 0.01) %>% runmodel

# Obtain network:
net1 <- getmatrix(ggm1, "omega")
sigma <- getmatrix(ggm1, "sigma")

eigen(sigma)$values
eigen(identity(nrow(net1)) - net1)$values

eigen(diag(nrow(net1)) - net1)$values





ordinal <- FALSE
nLevels <-  5
skewFactor <- 1
type <- "uniform"
missing <- 0
mod <- net1
n <- 500

# simulate data with seeds
ggm_dt_sim <- function(net, iter = 200, n = 500, ordinal = FALSE, nLevels = 4, skewFactor = 1, type =
                         c("uniform", "random"), missing = 0){

  # set up generator
  generator <- bootnet::ggmGenerator(ordinal = ordinal,
                                     nLevels = n_levels,
                                     skewFactor = skew_factor,
                                     type = type, missing = missing)

  dt <- lapply(seq_len(iter), function(i) {

    # set seed
    # seed <- i + 960308
    # set.seed(seed)

    # generate data
    data <- generator(n = 500, input = mod)
  })

  return(dt)

}



# a fit function that does not print output
silent_fit <- function(model) {
  # Capture and discard printed output, return the data frame
  capture.output(
    fit_result <- fit(model),
    file = NULL
  )
  return(fit_result)
}

# qgraph::qgraph(mod)

# ----- create misspecification -----

# add the additional parameter to the node with smallest predictability (most variance left to explain)
ggm_add <- function(mod){

  # calculate predictability
  r2 <- qgraph::centrality(mod, R2 = TRUE)$R2
  r2_sum <- matrix(r2, nrow = nrow(mod), ncol = ncol(mod)) + matrix(r2, nrow = nrow(mod), ncol = ncol(mod), byrow = TRUE)

  # avoid selecting diagnal and locations where there is already an edge
  r2_sum[mod != 0] <- 2
  diag(r2_sum) <- 2

  # location of smallest predictability
  loc <- which(r2_sum == min(r2_sum), arr.ind = TRUE)

  # remove duplication
  loc <- loc[(apply(loc, 1, sort) %>% t %>% duplicated()),]

  # in case there are multiple smallest, select one randomly
  # set.seed(960308)
  ind <- loc[sample(1:nrow(loc), 1),] %>% as.numeric

  # the smallest eigenvalue of the true is also negative
  # eigen(mod, symmetric = TRUE, only.values = TRUE)$values %>% min

  # the weight of the new edge is the same as the edge with the smallest absolute value in the network
  # purrr::map(data, function(dt) {
  edge_vec <- mod[mod!=0]
  add_edge_abs <- edge_vec[edge_vec %>% abs %>% which.min]
  # sign_edge <- ifelse(cor(dt[,ind[1]], dt[,ind[2]]) > 0, 1,-1)
  # * sign_edge
  #   return(mod)
  # })

  mod_misspec <- mod
  mod_misspec[ind[1], ind[2]] <- mod_misspec[ind[2], ind[1]] <- add_edge_abs

  set.seed(NULL)

  return(mod_misspec)
}

ggm_fit_misspec <- function(mod, iter) {

  # create adjacency matrix for CNA
  adj_net <- mod
  adj_net[adj_net != 0L] <- 1

  mod_misspec <- ggm_add(mod)

  # generate data for the misspecified model
  data <- ggm_dt_sim(net = mod_misspec, iter = iter)

  # getting the misspec model
  misspec_cna <- purrr::map(data, function(dt) {
    ggm(dt, omega = adj_net) %>% runmodel
  })

  # getting the fit of misspec model
  misspec_fit <- purrr::map_dfr(misspec_cna, function(m) {
    silent_fit(m) %>%
      filter(Measure %in% c("cfi", "rmsea", "tli")) %>%
      mutate(Measure = NULL, Value = round(Value, 3)) %>%
      t %>% as.data.frame %>%
      `colnames<-`(c("TLI_M","CFI_M","RMSEA_M")) %>%
      mutate(Model = "Misspec")
  }) %>% `rownames<-`(paste0("iter", 1:nrow(.))) %>% list

}

eigen(diag(nrow(mod)) - mod)$values

# true model --------------------------------------------------------------

ggm_fit_true <- function(mod) {

  # generate data from true
  data <- ggm_dt_sim(mod, iter)

  # create adjacency matrix for CNA
  adj_net <- mod
  adj_net[adj_net != 0L] <- 1

  # getting the true model
  true_cna <- purrr::map(data, function(dt) {
    ggm(dt, omega = adj_net) %>% runmodel
  })

  # getting the true model fit
  true_fit <- purrr::map_dfr(true_cna, function(m) {
    silent_fit(m) %>%
      filter(Measure %in% c("cfi", "rmsea", "tli")) %>%
      mutate(Measure = NULL, Value = round(Value, 3)) %>%
      t %>% as.data.frame %>%
      `colnames<-`(c("TLI_T","CFI_T","RMSEA_T")) %>%
      mutate(Model = "TRUE")
  }) %>% `rownames<-`(paste0("iter", 1:nrow(.))) %>% list

}

dt <- ggm_dt_sim(mod, iter = 100, n = 200, ordinal = FALSE)




misspec_fit <- ggm_fit_misspec(mod, iter = 100)

true_fit <- ggm_fit_true(mod)

misspec_sum <- purrr::map(misspec_fit,~dplyr::reframe(.,TLI_M=stats::quantile(TLI_M, c(seq(0.95,0,-0.01))),
                                                  RMSEA_M=stats::quantile(RMSEA_M, c(seq(0.05,1,0.01))),
                                                  CFI_M=stats::quantile(CFI_M, c(seq(0.95,0,-0.01)))))

#For the true model, compute the cutoffs (these will all be the same - just need in list form)
true_sum <- purrr::map(true_fit,~dplyr::reframe(.,TLI_T=stats::quantile(TLI_T, c(.05)),
                                               RMSEA_T=stats::quantile(RMSEA_T, c(.95)),
                                               CFI_T=stats::quantile(CFI_T, c(.05))))




