# obtain the fit distribution for the true and the misspecified model

# ----- ggm -----
# --- misspecified model ---
ggm_fit_misspec <- function(net, adj_net, iter, n = n, prop_pos,
                            ordinal, n_levels, skew_factor, type,
                            missing, par_fun, n_misspec, size_extra,
                            manual_size) {
  # create full list of candidate locations for the extra edge
  edge_cand_full <- edge_df(net = net)

  # select based on priority & non-replication
  edge_cand_ls <- select_edges(edge_cand_full,
    max_round = 10,
    n_misspec = n_misspec
  )

  # add
  mod_misspec <- ggm_add(
    net = net, edge_cand_ls = edge_cand_ls,
    n_misspec = n_misspec, prop_pos = prop_pos, n = n,
    size_extra = size_extra, manual_size = manual_size
  )

  # misspec model list
  m_net_ls <- mod_misspec$net_ls

  # progress bar across iter * n_misspec fits
  p <- progressr::progressor(
    steps = iter * length(m_net_ls),
    message = "Obtaining misspecified model fit distribution"
  )

  misspec_fit <-
    # generate data for the misspecified model
    par_fun(m_net_ls, function(net) {
      ggm_dt_sim(
        net = net, iter = iter, n = n, ordinal = ordinal,
        n_levels = n_levels, skew_factor = skew_factor,
        type = type, missing = missing, par_fun = par_fun
      ) %>% suppressWarnings()
    }) %>%
    # get cna model & fit based on the data
    par_fun(function(mod) {
      par_fun(mod, function(dt) {
        on.exit(p())
        fit_result <- psychonetrics::ggm(dt, omega = adj_net) %>%
          psychonetrics::runmodel(
            addMIs = FALSE,
            addSEs = FALSE,
            addInformation = FALSE
          ) %>%
          silent_fit()

        # in case model estimation fails, set fit to NA so other fits are returned
        tryCatch(
          {
            fit_result %>%
              filter(Measure %in% c("cfi", "rmsea", "tli", "srmr")) %>%
              mutate(Measure = NULL, Value = round(Value, 3)) %>%
              t() %>%
              as.data.frame() %>%
              `colnames<-`(c("TLI_M", "CFI_M", "RMSEA_M", "SRMR_M")) %>%
              mutate(
                Model = "misspec",
                TLI_M = case_when(TLI_M >= 1.2 | TLI_M < 0 ~ NA, TRUE ~ TLI_M),
                RMSEA_M = case_when(RMSEA_M >= 1 | RMSEA_M < 0 ~ NA, TRUE ~ RMSEA_M),
                CFI_M = case_when(CFI_M >= 1.2 | CFI_M < 0 ~ NA, TRUE ~ CFI_M),
                SRMR_M = case_when(SRMR_M >= 1 | SRMR_M < 0 ~ NA, TRUE ~ SRMR_M)
              )
          },
          error = function(e) {
            data.frame(
              TLI_M = NA,
              CFI_M = NA,
              RMSEA_M = NA,
              SRMR_M = NA,
              Model = "misspec"
            )
          }
        )
      }) %>%
        suppressWarnings() %>%
        bind_rows() %>%
        `rownames<-`(paste0("iter", 1:nrow(.)))
    }) %>% suppressWarnings()

  return(list(
    misspec_fit = misspec_fit,
    mod_misspec = m_net_ls,
    modified_edges = mod_misspec$modified_edges
  ))
}


# --- true model ---
ggm_fit_true <- function(net, adj_net, iter,
                         n, ordinal, n_levels,
                         skew_factor, type, missing,
                         par_fun) {
  # progress bar across iter fits
  p <- progressr::progressor(steps = iter,
                             message = "Obtaining true model fit distribution")

  true_fit <-
    # generate data from true
    ggm_dt_sim(
      net = net, iter = iter, n = n, ordinal = ordinal,
      n_levels = n_levels, skew_factor = skew_factor,
      type = type, missing = missing, par_fun = par_fun
    ) %>%
    # fit the model & obtain fit
    par_fun(function(dt) {
      on.exit(p())
      psychonetrics::ggm(dt, omega = adj_net) %>%
        psychonetrics::runmodel(
          addMIs = FALSE,
          addSEs = FALSE,
          addInformation = FALSE
        ) %>%
        silent_fit() %>%
        filter(Measure %in% c("cfi", "rmsea", "tli", "srmr")) %>%
        mutate(Measure = NULL, Value = round(Value, 3)) %>%
        t() %>%
        as.data.frame() %>%
        `colnames<-`(c("TLI_T", "CFI_T", "RMSEA_T", "SRMR_T")) %>%
        mutate(Model = "true")
    }) %>%
    suppressWarnings() %>%
    bind_rows() %>%
    `rownames<-`(paste0("iter", 1:nrow(.)))

  return(true_fit)
}
