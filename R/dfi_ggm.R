utils::globalVariables(c(
  "TLI_M", "RMSEA_M", "CFI_M", "TLI_T", "RMSEA_T", "CFI_T", "Power"
))

#' Dynamic fit index cutoffs for Gaussian Graphical Model
#'
#' This function determines the dynamic fit index cutoffs for Gaussian Graphical Models (GGM)
#'
#' @param net The empirical network to estimate dynamic fit index cutoffs for; Input should be a matrix
#' @param power The power that the cutoff value aims for
#' @param iter The number of iterations used in your simulation.
#' @param n Your sample size
#' @param prop_pos Decide the proportion of positive in the extra edges added to the empirical network to create the mis-specified network
#' @param ordinal Logical; should ordinal data be generated?
#' @param n_levels Number of levels used in ordinal data.
#' @param skew_factor How skewed should ordinal data be? 1 indicates uniform data and higher values increase skewedness.
#' @param type Should thresholds for ordinal data be sampled at random or determined uniformly?
#' @param missing Proportion of data that should be simulated to be missing.
#' @param ncores How many cores you want to use in the simulation. Recommend to leave one core free so that other tasks in the system are not impacted.
#' @param n_misspec Number of mis-specified model you want in the simulation. Default to 5, meaning
#' there are five mis-specified models with 1, 2, ..., 5 extra edges respectively.
#' Avoid setting too large numbers. Otherwise the simulation might fail or take too long.
#'
#' @return An object of class dfi_ggm. Can use summary to view a summary of results
#' @author Xinkai Du
#' Maintainer: Xinkai Du <xinkai.du.xd@gmail.com>
#'
#' @import dplyr future future.apply
#' @importFrom psychonetrics ggm
#' @export dfi_ggm

dfi_ggm <- function(net, power = 0.8, n_misspec = 5, iter = 200, n = 500, prop_pos = 0.8,
                    ordinal = FALSE, n_levels = 4, skew_factor = 1, min_extra = 0.15,
                    type = c("uniform", "random"), missing = 0, ncores = 1) {

  type <- match.arg(type, c("uniform", "random"))

  # Set up parallel execution plan
  if (ncores > 1) {
    future::plan(future::multisession, workers = ncores)
    par_fun <- future.apply::future_lapply %>% suppressWarnings
    on.exit(future::plan(future::sequential), add = TRUE)
  } else {
    par_fun <- lapply
  }

  # create adjacency matrix for CNA
  adj_net <- net
  adj_net[adj_net != 0L] <- 1

  # list to store results
  res <- list()

  # fit of misspec model
  misspec <- ggm_fit_misspec(net = net, iter = iter, n = n,
                             adj_net = adj_net, prop_pos = prop_pos,
                             ordinal = ordinal, n_levels = n_levels,
                             skew_factor = skew_factor, type = type,
                             missing = missing, par_fun = par_fun,
                             n_misspec = n_misspec, min_extra = min_extra)

  misspec_fit <- misspec$misspec_fit

  # add the added edges to the return list
  res$added_edges <- misspec$added_edges %>% setNames(paste0("L", 1:n_misspec))

  # true fit
  true_fit <- ggm_fit_true(net = net, iter = iter, n = n, adj_net = adj_net,
                           ordinal = ordinal, n_levels = n_levels, skew_factor = skew_factor,
                           type = type, missing = missing, par_fun = par_fun)

  # store true and misspec fit
  res$fit <- fit_data(true_fit = true_fit, misspec_fit, par_fun = par_fun)

  # store misspec models
  res$mod_misspec <- misspec$mod_misspec %>%
    setNames(paste0("L", 1:length(misspec$mod_misspec)))

  # find the 0-95th percentile in the fit dist of misspec model;
  # this is to find the value that shows better fit than 95% of the values in the misspec dist
  # 95% percentile of TLI & CFI; 5% of RMSEA
  # ideally these values also performs worse than 95% of the values in the true distribution (if power = 0.95)
  misspec_sum <- par_fun(misspec_fit, function(df) {
    dplyr::reframe(df,TLI_M=stats::quantile(TLI_M, c(seq(0.95,0,-0.01))),
                   RMSEA_M=stats::quantile(RMSEA_M, c(seq(0.05,1,0.01))),
                   CFI_M=stats::quantile(CFI_M, c(seq(0.95,0,-0.01))))
  })

  # find in the true distribution the value that marks power = power
  true_sum <- par_fun(true_fit, function(df) {
    dplyr::reframe(df,TLI_T=stats::quantile(TLI_T, c(1-power)),
                   RMSEA_T=stats::quantile(RMSEA_T, c(power)),
                   CFI_T=stats::quantile(CFI_T, c(1-power)))
  })

  # calculate sensitivity based on misspec dist & power
  final <- par_fun(1:length(misspec_fit), function(i) {
    # Create data frame directly
    fit_i <- data.frame(misspec_sum[[i]], true_sum, Power = seq(.95, 0.0, -.01))

    # Add comparison columns
    fit_i$TLI <- ifelse(fit_i$TLI_M <= fit_i$TLI_T, 1, 0)
    fit_i$RMSEA <- ifelse(fit_i$RMSEA_M >= fit_i$RMSEA_T, 1, 0)
    fit_i$CFI <- ifelse(fit_i$CFI_M <= fit_i$CFI_T, 1, 0)

    # Process each fit index and combine directly
    tli_subset <- subset(fit_i, !duplicated(TLI) | Power == 0, c("TLI_M", "Power", "TLI"))
    tli <- tli_subset[tli_subset$TLI == 1 | tli_subset$Power == 0, 1:2]
    colnames(tli) <- c("TLI", "sensitivity_tli")

    rmsea_subset <- subset(fit_i, !duplicated(RMSEA) | Power == 0, c("RMSEA_M", "Power", "RMSEA"))
    rmsea <- rmsea_subset[rmsea_subset$RMSEA == 1 | rmsea_subset$Power == 0, 1:2]
    colnames(rmsea) <- c("RMSEA", "sensitivity_rmsea")

    cfi_subset <- subset(fit_i, !duplicated(CFI) | Power == 0, c("CFI_M", "Power", "CFI"))
    cfi <- cfi_subset[cfi_subset$CFI == 1 | cfi_subset$Power == 0, 1:2]
    colnames(cfi) <- c("CFI", "sensitivity_cfi")

    # Create and return only the final result
    cbind(tli[1,], rmsea[1,], cfi[1,])[c("TLI", "sensitivity_tli", "RMSEA", "sensitivity_rmsea", "CFI", "sensitivity_cfi")]
  })

  # dynamic cutoff that correspond to the chosen power
  L0 <- data.frame(cbind(true_sum[[1]]$TLI_T,power,
                         true_sum[[1]]$RMSEA_T,power,
                         true_sum[[1]]$CFI_T),power) %>%
    `colnames<-`(c("TLI","power_tli","RMSEA","power_rmsea","CFI","power_cfi")) %>%
    round(3) %>% mutate(mis_level = "L0")

  # fit value cutoffs & their sensitivity to the misspec models
  Ls <- bind_rows(final) %>% round(3) %>% mutate(mis_level = paste0("L",1:nrow(.)))

  # ditch the cutoff if its sensitivity < 0.5
  for (j in 1:length(misspec_sum)) {
    if(Ls[j,2]<.50){Ls[j,1]<-NA}
    if(Ls[j,4]<.50){Ls[j,3]<-NA}
    if(Ls[j,6]<.50){Ls[j,5]<-NA}
  }

  res$cutoff_true <- L0 # this is the cutoff that marks power = power
  res$cutoff_misspec <- Ls # these are the cutoff that marks sensitivity & power = power;
  # this is the best fitting value in misspec that is still lower than the worst fitting value in true

  # assign class to summarize
  class(res) <- "dfi_ggm"

  return(res)

}


#' #' The invariance partial pruning (IVPP) algorithm for panel GVAR models
#' #'
#' #' Plot the fit distribution
#' #'
#' #' @param net The empirical network to estimate dynamic fit index cutoffs for
#' #' @param power The power that the cutoff value aims for
#' #' @param iter The number of iterations used in your simulation.
#' #' @param n Your sample size
#' #' @param prop_pos Decide the proportion of positive in the extra edges added to the empirical network to create the mis-specified network
#' #' @param ordinal Logical; should ordinal data be generated?
#' #' @param n_levels Number of levels used in ordinal data.
#' #' @param skew_factor How skewed should ordinal data be? 1 indicates uniform data and higher values increase skewedness.
#' #' @param type Should thresholds for ordinal data be sampled at random or determined uniformly?
#' #' @param missing Proportion of data that should be simulated to be missing.
#' #' @param ncores How many cores you want to use in the simulation. Recommend to leave one core free so that other tasks in the system are not impacted.
#' #' @param n_misspec Number of mis-specified model you want in the simulation. Default to 5, meaning
#' #' there are five mis-specified models with 1, 2, ..., 5 extra edges respectively.
#' #' Avoid setting too large numbers. Otherwise the simulation might fail or take too long.
#' #'
#' #' @return An object of class dfi_ggm. Can use summary to view a summary of results
#' #'
#' #'
#' #' @import dplyr future future.apply
#' #' @importFrom bootnet ggmGenerator
#' #' @importFrom psychonetrics ggm
#' #' @export dfi_ggm
#' #'
#'
#' plot.dfi_ggm <- function() {
#'
#' }
