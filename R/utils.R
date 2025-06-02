
utils::globalVariables(c(
  ".", "weights", "loc", "Measure", "Value"
))

#' @import dplyr future future.apply
#' @importFrom psychonetrics ggm
#' @importFrom stats setNames quantile
#' @importFrom utils capture.output
#' @importFrom bootnet ggmGenerator
#' @importFrom qgraph centrality

# ----- a fit function that does not print output -----
silent_fit <- function(model) {
  # Capture and discard printed output, return the data frame
  utils::capture.output(
    fit_result <- psychonetrics::fit(model),
    file = NULL
  )
  return(fit_result)
}

# ----- simulate data with seeds -----
ggm_dt_sim <- function(net, iter, n, ordinal, n_levels, skew_factor, type, missing, par_fun){

  # set up generator
  generator <- bootnet::ggmGenerator(ordinal = ordinal,
                                     nLevels = n_levels,
                                     skewFactor = skew_factor,
                                     type = type, missing = missing)

  # generate data
  par_fun(seq_len(iter), function(i) {
    generator(n = 500, input = net)
  })

}

# ----- net2vec -----
# obtain the lower.tri of a network, label each edge with the nodes that it connects, &
# output the vector of edges annotated with corresponding nodes
net2vec <- function(net){

  # assign row&colnames if there's none
  if (is.null(rownames(net))) {
    rownames(net) <- paste0("V", 1:nrow(net))
  }

  if (is.null(colnames(net))) {
    colnames(net) <- paste0("V", 1:ncol(net))
  }

  # obtain lower.tri
  vec <- data.frame(net[lower.tri(net)])
  names(vec) <- "weights"

  # annotate values with row&colnames
  names.Mat <- outer(rownames(net), colnames(net), paste, sep = "--")
  vec$loc <- names.Mat[lower.tri(names.Mat)]

  return(vec)

}

# ----- edge_df -----
# create a full list of locations that don't have an edge yet
# summed predictability of neighboring nodes that these locations connect are included
edge_df <- function(net){

  # calculate predictability of each edge
  r2 <- qgraph::centrality(net, R2 = TRUE)$R2

  # calculate the summed predictability of nodes connected by each edge
  r2_sum <- matrix(r2, nrow = nrow(net), ncol = ncol(net)) +
    matrix(r2, nrow = nrow(net), ncol = ncol(net), byrow = TRUE)

  r2_sum[net!=0] <- 2
  diag(r2_sum) <- 2

  # browser()

  # edge_df <- net2vec(r2_sum) %>%
  #   filter(weights != 2L) %>%
  #   arrange(desc(weights)) %>%
  #   mutate(from = as.numeric(gsub("V(\\d+)--V\\d+", "\\1", loc)),
  #          to = as.numeric(gsub("V\\d+--V(\\d+)", "\\1", loc)))

  # # avoid selecting diagnal and locations where there is already an edge
  # r2_sum[net != 0] <- -2
  # diag(r2_sum) <- -2
  #
  # # rank priorities based on summed predictability; here the weights are summed predictabilities
  # # edges that connect nodes of high predictability are prioritized because they impact the implied covariance to the largest extent
  edge_df <- net2vec(r2_sum) %>%
    filter(weights != 2L) %>%
    arrange(desc(weights)) %>%
    mutate(priority = dense_rank(desc(weights)),
           from = as.numeric(gsub("V(\\d+)--V\\d+", "\\1", loc)),
           to = as.numeric(gsub("V\\d+--V(\\d+)", "\\1", loc)))
  #
  # return(edge_df)

}

# ----- select_edges -----
# select candidate edges from the full edge list
select_edges <- function(edge_df, max_round, n_misspec) {

  edge_ls <- vector(mode = "list", length = max_round)

  for(i in 1:max_round){
    for(j in 1:n_misspec){
      # start from the ith row, select j rows
      edge_ls[[i]][[j]] <- edge_df[i:nrow(edge_df),][seq_len(j),]
    }
  }

  return(edge_ls)
}

# ----- create misspecification -----
# add the additional parameter to the node with smallest predictability (most variance left to explain)
ggm_add <- function(net, edge_cand_ls, n_misspec, prop_pos, min_extra) {
    # Extract non-zero elements from the original network
    edge_vec <- net[net != 0]

    # Get minimum absolute value of edge weights to use for new edges
    # set to be no smaller than min_extra
    # min_abs_edge <- max(min(abs(edge_vec)), min_extra)

    # Create a list to store resulting networks
    result_networks <- list()

    # Try each starting point
    for (start_idx in 1:length(edge_cand_ls)) {
      # Flag to check if all networks for this starting point are positive definite
      all_positive_definite <- TRUE

      # Temporary list to store networks for this starting point
      net_ls <- list()

      # a list to store added edges and their locations
      edges_to_add <-  list()

      # For each model size (1 to n_misspec)
      for (k in 1:n_misspec) {
        # Get the edges for this model
        edges_to_add[[k]] <- edge_cand_ls[[start_idx]][[k]]

        # Create a copy of the original network to modify
        mod_misspec <- net

        # Generate edge weights with random signs
        edges_to_add[[k]]$added_edges <- min_extra * sample(c(1, -1),
                                                            nrow(edges_to_add[[k]]),
                                                            prob = c(prop_pos, 1-prop_pos),
                                                            replace = TRUE)

        # Add each edge with the determined weight
        for (i in 1:nrow(edges_to_add[[k]])) {
          from_node <- edges_to_add[[k]]$from[i]
          to_node <- edges_to_add[[k]]$to[i]

          # Add the edge in both directions (symmetric matrix)
          mod_misspec[from_node, to_node] <- mod_misspec[to_node, from_node] <- edges_to_add[[k]]$added_edges[i]
        }

        # Check the positive definiteness of I - omega
        eigen_vals <- eigen((diag(nrow(mod_misspec)) - mod_misspec), symmetric = TRUE)$values

        if (all(eigen_vals > 0)) {
          # Store the positive definite network
          net_ls[[k]] <- mod_misspec
        } else {
          # Not all networks are positive definite, break and try next starting point
          all_positive_definite <- FALSE
          break
        }
      }

      # If all networks for this starting point are positive definite, return them
      if (all_positive_definite) {
        return(list(
          net_ls = net_ls,
          added_edges = lapply(edges_to_add, function(x) x %>% select(loc, added_edges))
        ))
      }
    }

    # If we've tried all starting points and none produced all positive definite networks
    stop("can't find networks that are all positive definite after 10 iter")
}

ggm_fit_misspec <- function(net, adj_net, iter, n = n, prop_pos, ordinal, n_levels, skew_factor, type, missing, par_fun, n_misspec, min_extra) {

  # create full list of candidate locations for the extra edge
  edge_cand_full <- edge_df(net = net)

  # select based on priority & non-replication (don't select edges that are connected to the same node)
  edge_cand_ls <- select_edges(edge_cand_full,
                               max_round = 10,
                               n_misspec = n_misspec)

  # add
  mod_misspec <- ggm_add(net = net, edge_cand_ls = edge_cand_ls,
                         n_misspec = n_misspec, prop_pos = prop_pos, min_extra = min_extra)

  # misspec model list
  m_net_ls <- mod_misspec$net_ls

  misspec_fit <-
    # generate data for the misspecified model
    par_fun(m_net_ls, function(net) {
      ggm_dt_sim(net = net, iter = iter, n = n, ordinal = ordinal,
                 n_levels = n_levels, skew_factor = skew_factor,
                 type = type, missing = missing, par_fun = par_fun) %>% suppressWarnings
    }) %>%
    # get cna model & fit based on the data
    par_fun(function(mod) {
      par_fun(mod, function(dt) {
        fit_result <- psychonetrics::ggm(dt, omega = adj_net) %>%
          psychonetrics::runmodel(addMIs = FALSE,
                                  addSEs = FALSE,
                                  addInformation = FALSE) %>% silent_fit

        # in case model estimation fails, set fit to NA so that other fit are returned
        tryCatch({
          fit_result %>%
            filter(Measure %in% c("cfi", "rmsea", "tli")) %>%
            mutate(Measure = NULL, Value = round(Value, 3)) %>%
            t %>% as.data.frame %>%
            `colnames<-`(c("TLI_M","CFI_M","RMSEA_M")) %>%
            mutate(Model = "misspec",
                   TLI_M = case_when(TLI_M >= 1.2 | TLI_M < 0 ~ NA, TRUE ~ TLI_M),
                   RMSEA_M = case_when(RMSEA_M >= 1 | RMSEA_M < 0 ~ NA, TRUE ~ RMSEA_M),
                   CFI_M = case_when(CFI_M >= 1.2 | CFI_M < 0 ~ NA, TRUE ~ CFI_M))
        }, error = function(e) {
          data.frame(TLI_M = NA, CFI_M = NA, RMSEA_M = NA, Model = "misspec")
        })
      }) %>% suppressWarnings %>% bind_rows %>%
        `rownames<-`(paste0("iter", 1:nrow(.)))
    }) %>% suppressWarnings

  return(list(
    misspec_fit = misspec_fit,
    mod_misspec = m_net_ls,
    added_edges = mod_misspec$added_edges
  ))

}

# true model --------------------------------------------------------------

ggm_fit_true <- function(net, adj_net, iter, n, ordinal, n_levels, skew_factor, type, missing, par_fun) {

  true_fit <-
    # generate data from true
    ggm_dt_sim(net = net, iter = iter, n = n, ordinal = ordinal,
               n_levels = n_levels, skew_factor = skew_factor,
               type = type, missing = missing, par_fun = par_fun) %>%
    # fit the model & obtain fit
    par_fun(function(dt) {
      psychonetrics::ggm(dt, omega = adj_net) %>%
        psychonetrics::runmodel(addMIs = FALSE,
                                addSEs = FALSE,
                                addInformation = FALSE) %>%
        silent_fit %>%
        filter(Measure %in% c("cfi", "rmsea", "tli")) %>%
        mutate(Measure = NULL, Value = round(Value, 3)) %>%
        t %>% as.data.frame %>%
        `colnames<-`(c("TLI_T","CFI_T","RMSEA_T")) %>%
        mutate(Model = "true")
  }) %>% suppressWarnings %>%
    bind_rows %>%
    `rownames<-`(paste0("iter", 1:nrow(.)))

  return(true_fit)

}

fit_data <- function(true_fit, misspec_fit, par_fun){

  df_results <- par_fun(misspec_fit,function(df) {
    cbind(df,true_fit)
  })

  #Create beginning of variable name for each
  dat_name <- rep(c("TLI_L","CFI_L","RMSEA_L","Type_L"),2)

  #Create vector of 0's for the True model
  dat_0 <- rep(0,4)

  #Get number of levels of misspecification
  dat_lev <- length(df_results)

  #Create combo of Level #'s and 0's to merge with variable names (in list form)
  dat_num <- list()
  for (i in 1:dat_lev){
    output <- c(rep(i,4),dat_0)
    dat_num[[i]] <- output
  }

  #Combine variable name with level # (in list form)
  var_names <- par_fun(dat_num, function(x){
    paste0(dat_name,x)
  })

  #Rename variables in dataset list
  dat_revised <- par_fun(seq_along(df_results), function(x){
    colnames(df_results[[x]]) <- var_names[[x]]
    #not sure why I need to mention it again but I do
    df_results[[x]]
  })

  #Combine into one dataset
  df_renamed <- do.call(cbind.data.frame,dat_revised)

  #Remove the duplicate L0 info
  df_renamed2 <- df_renamed[!duplicated(colnames(df_renamed))]

  #Remove the "type" variable (unnecessary)
  df_renamed3 <- df_renamed2[, -grep("Type",colnames(df_renamed2))]

  #Reorder variables
  df_renamed4 <- df_renamed3 %>%
    dplyr::relocate(order(colnames(df_renamed3))) %>%            #gets numbers in order
    dplyr::relocate(dplyr::starts_with("TLI"),dplyr::starts_with("RMSEA"))    #gets fit indices in order

  return(df_renamed4)
}
