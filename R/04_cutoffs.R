

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
#' @param min_extra Minimum absolute edge weight used when adding extra edges to create misspecified models.
#' @param type Should thresholds for ordinal data be sampled at random or determined uniformly?
#' @param missing Proportion of data that should be simulated to be missing.
#' @param ncores How many cores you want to use in the simulation. Recommend to leave one core free so that other tasks in the system are not impacted.
#' @param n_misspec Number of mis-specified model you want in the simulation. Default to 3, meaning
#' there are five mis-specified models with 1, 2, ..., 5 extra edges respectively.
#' Avoid setting too large numbers. Otherwise the simulation might fail or take too long.
#'
#' @return An object of class dfi_ggm. Can use summary to view a summary of results
#' @author Xinkai Du
#' Maintainer: Xinkai Du <xinkai.du.xd@gmail.com>
#'
#' @export dfi_ggm
#' @examples
#' \dontrun{
#' library(psych)
#' library(psychonetrics)
#' library(dplyr)
#' library(netDFI)
#'
#' # get the big five inventory data from psych
#' data(bfi)
#'
#' # estimate ggm
#' bfi_mod <- ggm(bfi) %>% prune %>% runmodel
#'
#' # obtain the partial correlation matrix
#' bfi_net <- getmatrix(bfi_mod, "omega")
#'
#' # allow future_apply to use more memory
#' options(future.globals.maxSize = 2 * 1024^3)
#'
#' # run dfi
#' dfi_bfi <- dfi_ggm(
#'   bfi_net,
#'   ncores = parallel::detectCores(),
#'   power = 0.80,
#'   iter = 200,
#'   n_misspec = 2
#' )
#' dfi_bfi
#'
#' # plot results
#' p <- plot(dfi_bfi)
#' p[[1]]
#' }

dfi_ggm <- function(net, power = 0.95, n_misspec = 3, iter = 500, n = 500, prop_pos = 0.8,
                    ordinal = FALSE, n_levels = 4, skew_factor = 1, min_extra = 0.2,
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

  # assign row&colnames
  rownames(net) <- colnames(net) <- paste0("V", 1:nrow(net))

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
    reframe(df,TLI_M=stats::quantile(TLI_M, c(seq(0.95,0,-0.01)), na.rm = TRUE),
            RMSEA_M=stats::quantile(RMSEA_M, c(seq(0.05,1,0.01)), na.rm = TRUE),
            CFI_M=stats::quantile(CFI_M, c(seq(0.95,0,-0.01)), na.rm = TRUE))
  })

  # find in the true distribution the value that marks power = power
  true_cut <- true_fit %>%
    reframe(TLI_T=stats::quantile(TLI_T, c(1-power)),
            RMSEA_T=stats::quantile(RMSEA_T, c(power)),
            CFI_T=stats::quantile(CFI_T, c(1-power)))

  # calculate sensitivity based on misspec dist & power
  final <- par_fun(1:length(misspec_fit), function(i) {
    # Create data frame directly
    fit_i <- data.frame(misspec_sum[[i]], true_cut, Power = seq(.95, 0.0, -.01))

    # Add comparison columns
    fit_i$TLI_sens <- ifelse(fit_i$TLI_M <= fit_i$TLI_T, 1, 0)
    fit_i$RMSEA_sens <- ifelse(fit_i$RMSEA_M >= fit_i$RMSEA_T, 1, 0)
    fit_i$CFI_sens <- ifelse(fit_i$CFI_M <= fit_i$CFI_T, 1, 0)

    # HB cutoffs
    fit_i$TLI_HB <- 0.95
    fit_i$RMSEA_HB <- 0.05
    fit_i$CFI_HB <- 0.95

    # HB sensitivity
    fit_i$TLI_HB_sens <- ifelse(fit_i$TLI_M <= 0.95, 1, 0)
    fit_i$RMSEA_HB_sens <- ifelse(fit_i$RMSEA_M >= 0.05, 1, 0)
    fit_i$CFI_HB_sens <- ifelse(fit_i$CFI_M <= 0.95, 1, 0)

    # Calculate the sensitivity of DFI
    tli <- fit_i %>% filter(!duplicated(TLI_sens) | Power == 0) %>%
      select(TLI_M, Power, TLI_sens) %>%
      filter(TLI_sens == 1 | Power == 0) %>%
      select(TLI_M, Power) %>%
    `colnames<-`(c("TLI", "sensitivity_tli"))

    rmsea <- fit_i %>% filter(!duplicated(RMSEA_sens) | Power == 0) %>%
      select(RMSEA_M, Power, RMSEA_sens) %>%
      filter(RMSEA_sens == 1 | Power == 0) %>%
      select(RMSEA_M, Power) %>%
      `colnames<-`(c("RMSEA", "sensitivity_rmsea"))

    cfi <- fit_i %>% filter(!duplicated(CFI_sens) | Power == 0) %>%
      select(CFI_M, Power, CFI_sens) %>%
      filter(CFI_sens == 1 | Power == 0) %>%
      select(CFI_M, Power) %>%
      `colnames<-`(c("CFI", "sensitivity_cfi"))

    # calculate the sensitivity of HB cutoffs
    tli_HB <- fit_i %>% filter(!duplicated(TLI_HB_sens) | Power == 0) %>%
      select(TLI_HB, Power, TLI_HB_sens) %>%
      filter(TLI_HB_sens == 1 | Power == 0) %>%
      select(TLI_HB, Power) %>%
      `colnames<-`(c("TLI_HB", "sensitivity_tli_HB"))

    rmsea_HB <- fit_i %>% filter(!duplicated(RMSEA_HB_sens) | Power == 0) %>%
      select(RMSEA_HB, Power, RMSEA_HB_sens) %>%
      filter(RMSEA_HB_sens == 1 | Power == 0) %>%
      select(RMSEA_HB, Power) %>%
      `colnames<-`(c("RMSEA_HB", "sensitivity_rmsea_HB"))

    cfi_HB <- fit_i %>% filter(!duplicated(CFI_HB_sens) | Power == 0) %>%
      select(CFI_HB, Power, CFI_HB_sens) %>%
      filter(CFI_HB_sens == 1 | Power == 0) %>%
      select(CFI_HB, Power) %>%
      `colnames<-`(c("CFI_HB", "sensitivity_cfi_HB"))

    # Create and return only the final result
    cbind(
      tli[1,], rmsea[1,], cfi[1,],
      tli_HB[1,], rmsea_HB[1,], cfi_HB[1,]
      )[c("TLI", "sensitivity_tli",
          "RMSEA", "sensitivity_rmsea",
          "CFI", "sensitivity_cfi",
          "TLI_HB", "sensitivity_tli_HB",
          "RMSEA_HB", "sensitivity_rmsea_HB",
          "CFI_HB", "sensitivity_cfi_HB")]
  })

  # find the power percentile for true
  true_sum <- true_fit %>%
    reframe(TLI_T=stats::quantile(TLI_T, c(seq(0.05,1,0.01))),
            RMSEA_T=stats::quantile(RMSEA_T, c(seq(0.95,0,-0.01))),
            CFI_T=stats::quantile(CFI_T, c(seq(0.05,1,0.01)))) %>%
    mutate(Power = seq(.95, 0.0, -.01))

  # HB cutoffs
  true_sum$TLI_HB <- 0.95
  true_sum$RMSEA_HB <- 0.05
  true_sum$CFI_HB <- 0.95

  # HB power
  true_sum$TLI_HB_spec <- ifelse(true_sum$TLI_T >= 0.95, 1, 0)
  true_sum$RMSEA_HB_spec <- ifelse(true_sum$RMSEA_T <= 0.05, 1, 0)
  true_sum$CFI_HB_spec <- ifelse(true_sum$CFI_T >= 0.95, 1, 0)

  # calculate the sensitivity of HB cutoffs
  tli_HB <- true_sum %>% filter(!duplicated(TLI_HB_spec) | Power == 0) %>%
    select(TLI_HB, Power, TLI_HB_spec) %>%
    filter(TLI_HB_spec == 1 | Power == 0) %>%
    select(TLI_HB, Power) %>%
    `colnames<-`(c("TLI_HB", "power_tli_HB"))

  rmsea_HB <- true_sum %>% filter(!duplicated(RMSEA_HB_spec) | Power == 0) %>%
    select(RMSEA_HB, Power, RMSEA_HB_spec) %>%
    filter(RMSEA_HB_spec == 1 | Power == 0) %>%
    select(RMSEA_HB, Power) %>%
    `colnames<-`(c("RMSEA_HB", "power_rmsea_HB"))

  cfi_HB <- true_sum %>% filter(!duplicated(CFI_HB_spec) | Power == 0) %>%
    select(CFI_HB, Power, CFI_HB_spec) %>%
    filter(CFI_HB_spec == 1 | Power == 0) %>%
    select(CFI_HB, Power) %>%
    `colnames<-`(c("CFI_HB", "power_cfi_HB"))

  # L0 specificity to misspecification
  L0 <- cbind(
    true_cut$TLI_T,power,
    true_cut$RMSEA_T,power,
    true_cut$CFI_T,power,
    tli_HB[1,], rmsea_HB[1,], cfi_HB[1,]
  ) %>%
    `names<-`(c("TLI", "specificity_tli",
                "RMSEA", "specificity_rmsea",
                "CFI", "specificity_cfi",
                "TLI_HB", "specificity_tli_HB",
                "RMSEA_HB", "specificity_rmsea_HB",
                "CFI_HB", "specificity_cfi_HB")) %>%
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

  invisible(res)

}


