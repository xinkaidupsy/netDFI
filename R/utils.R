
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
  dt <- par_fun(seq_len(iter), function(i) {
    data <- generator(n = 500, input = net)
  })

  return(dt)

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

  # avoid selecting diagnal and locations where there is already an edge
  r2_sum[net != 0] <- -2
  diag(r2_sum) <- -2

  # rank priorities based on summed predictability; here the weights are summed predictabilities
  # edges that connect nodes of low predictability are prioritized
  edge_df <- net2vec(r2_sum) %>%
    arrange(desc(weights)) %>%
    mutate(priority = dense_rank(desc(weights)),
           from = as.numeric(gsub("V(\\d+)--V\\d+", "\\1", loc)),
           to = as.numeric(gsub("V\\d+--V(\\d+)", "\\1", loc))) %>%
    filter(weights != -2L)

  return(edge_df)

}

# ----- select_edges -----
# select candidate edges from the full edge list
select_edges <- function(edge_df, max_round, n_misspec) {

  # Create a list to store results for different starting points
  all_results <- list()

  # For each starting row
  for (start_idx in 1:min(max_round, nrow(edge_df))) {
    # Create a sub-list for different edge counts
    model_results <- list()

    # Reorder edges to start from start_idx
    reordered_edges <- rbind(
      edge_df[start_idx:nrow(edge_df), ],
      if(start_idx > 1) edge_df[1:(start_idx-1), ] else NULL
    )

    # For each requested size (1 to n_misspec)
    for (k in 1:n_misspec) {
      # Initialize variables
      selected_edges <- edge_df[0, ]  # Empty data frame with same structure
      used_nodes <- c()

      # Go through each edge in the reordered dataframe
      for (i in 1:nrow(reordered_edges)) {
        # Get the current row
        current_edge <- reordered_edges[i, ]

        # Get the nodes of this edge
        from_node <- current_edge$from
        to_node <- current_edge$to

        # Check if either node is already used
        if (!(from_node %in% used_nodes || to_node %in% used_nodes)) {
          # Add this edge
          selected_edges <- rbind(selected_edges, current_edge)
          used_nodes <- c(used_nodes, from_node, to_node)
        }

        # If we've selected k edges, stop
        if (nrow(selected_edges) == k) {
          break
        }
      }

      # Add the selected edges to results
      model_results[[k]] <- selected_edges
    }

    # Store all models for this starting point
    all_results[[start_idx]] <- model_results
  }

  return(all_results)
}

# ----- create misspecification -----
# add the additional parameter to the node with smallest predictability (most variance left to explain)
ggm_add <- function(net, edge_cand_ls, n_misspec, prop_pos) {
    # Extract non-zero elements from the original network
    edge_vec <- net[net != 0]

    # Get minimum absolute value of edge weights to use for new edges
    min_abs_edge <- min(abs(edge_vec))

    # Create a list to store resulting networks
    result_networks <- list()

    # Try each starting point
    for (start_idx in 1:length(edge_cand_ls)) {
      # Flag to check if all networks for this starting point are positive definite
      all_positive_definite <- TRUE

      # Temporary list to store networks for this starting point
      net_ls <- list()

      # For each model size (1 to n_misspec)
      for (k in 1:n_misspec) {
        # Get the edges for this model
        edges_to_add <- edge_cand_ls[[start_idx]][[k]]

        # Create a copy of the original network to modify
        mod_misspec <- net

        # Generate edge weights with random signs
        edge_weights <- min_abs_edge * sample(c(1, -1), nrow(edges_to_add),
                                              prob = c(prop_pos, 1-prop_pos), replace = TRUE)

        # Add each edge with the determined weight
        for (i in 1:nrow(edges_to_add)) {
          from_node <- edges_to_add$from[i]
          to_node <- edges_to_add$to[i]

          # Add the edge in both directions (symmetric matrix)
          mod_misspec[from_node, to_node] <- mod_misspec[to_node, from_node] <- edge_weights[i]
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
        return(net_ls)
      }
    }

    # If we've tried all starting points and none produced all positive definite networks
    stop("can't find networks that are all positive definite after 10 iter")
}

ggm_fit_misspec <- function(net, adj_net, iter, n = n, prop_pos, ordinal, n_levels, skew_factor, type, missing, par_fun, n_misspec) {

  # create full list of candidate locations for the extra edge
  edge_cand_full <- edge_df(net = net)

  # select based on priority & non-replication (don't select edges that are connected to the same node)
  edge_cand_ls <- select_edges(edge_cand_full,
                               max_round = 10,
                               n_misspec = n_misspec)

  # add
  mod_misspec <- ggm_add(net = net, edge_cand_ls = edge_cand_ls,
                         n_misspec = n_misspec, prop_pos = prop_pos)

  # generate data for the misspecified model
  data <- par_fun(mod_misspec, function(net) {
    ggm_dt_sim(net = net, iter = iter, n = n, ordinal = ordinal,
               n_levels = n_levels, skew_factor = skew_factor,
               type = type, missing = missing, par_fun = par_fun) %>% suppressWarnings
  }) %>% suppressWarnings

  # getting the cna that fits the empirical model to the data simulated from the misspec model
  misspec_cna <- par_fun(data, function(dt) {
    m_c <- par_fun(dt, function(d) {
      m_c <- psychonetrics::ggm(d, omega = adj_net) %>% psychonetrics::runmodel()
    }) %>% suppressWarnings
  }) %>% suppressWarnings

  # getting the fit of misspec model
  misspec_fit <- par_fun(misspec_cna, function(mc) {
    mf <- par_fun(mc, function(m){
      silent_fit(m) %>%
        filter(Measure %in% c("cfi", "rmsea", "tli")) %>%
        mutate(Measure = NULL, Value = round(Value, 3)) %>%
        t %>% as.data.frame %>%
        `colnames<-`(c("TLI_M","CFI_M","RMSEA_M")) %>%
        mutate(Model = "misspec")
    }) %>% bind_rows %>%
      `rownames<-`(paste0("iter", 1:nrow(.)))
  })

  return(list(
    misspec_fit = misspec_fit,
    mod_misspec = mod_misspec
  ))

}

# true model --------------------------------------------------------------

ggm_fit_true <- function(net, adj_net, iter, n, ordinal, n_levels, skew_factor, type, missing, par_fun) {

  # generate data from true
  data <- ggm_dt_sim(net = net, iter = iter, n = n, ordinal = ordinal,
                     n_levels = n_levels, skew_factor = skew_factor,
                     type = type, missing = missing, par_fun = par_fun) %>% suppressWarnings

  # getting the true model
  true_cna <- par_fun(data, function(dt) {
    psychonetrics::ggm(dt, omega = adj_net) %>% psychonetrics::runmodel()
  }) %>% suppressWarnings

  # getting the true model fit
  true_fit <- par_fun(true_cna, function(m) {
    silent_fit(m) %>%
      filter(Measure %in% c("cfi", "rmsea", "tli")) %>%
      mutate(Measure = NULL, Value = round(Value, 3)) %>%
      t %>% as.data.frame %>%
      `colnames<-`(c("TLI_T","CFI_T","RMSEA_T")) %>%
      mutate(Model = "true")
  }) %>% bind_rows %>%
    `rownames<-`(paste0("iter", 1:nrow(.))) %>% list

}

fit_data <- function(true_fit, misspec_fit, par_fun){

  df_results <- par_fun(misspec_fit,function(df) {
    cbind(df,true_fit)
  })

  #Create beginning of variable name for each
  dat_name <- rep(c("TLI_L","RMSEA_L","CFI_L","Type_L"),2)

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


#
#
# fit1<-unlist(Fit)%>% matrix(nrow=(length(misspec_sum)+1), ncol=6) %>%
#   `colnames<-`(c("TLI","sensitivity_tli","RMSEA","sensitivity_rmsea","CFI","sensitivity_cfi"))
#
# # ditch the cutoff if its sensitivity < 0.5
# for (j in 2:(length(misspec_sum)+1)) {
#   if(fit1[j,2]<.50){fit1[j,1]<-NA}
#   if(fit1[j,4]<.50){fit1[j,3]<-NA}
#   if(fit1[j,6]<.50){fit1[j,5]<-NA}
# }
#
# #create blanks to make table easier to read
# pp<-c(rep("--",(length(misspec_sum)+1)))
# pp0<-c(rep("",(length(misspec_sum)+1)))
#
# # #matrix of cross-loadings added at each level
# # mag <- ggm_add(net1, prop_pos = 0.6) %>%
# #   tidyr::separate(V1,into=c("a","b","Magnitude","d","e"),sep=" ") %>%
# #   select(Magnitude) %>%
# #   mutate(Magnitude=as.numeric(Magnitude),
# #          Magnitude=round(Magnitude,digits=3))
#
# # #Create column of cross-loadings added at each level
# # mname<-c("NONE","","")
# # for (i in 1:length(misspec_sum)){
# #   mi<-c(mag[i,],"","")
# #   mname<-cbind(mname,mi)
# # }
#
# #format column for each index and the misspecification magnitude
# TT<-noquote(matrix(rbind(fit1[,1],fit1[,2],pp0),ncol=1))
# RR<-noquote(matrix(rbind(fit1[,3],fit1[,4],pp0),ncol=1))
# CC<-noquote(matrix(rbind(fit1[,5],fit1[,6],pp0),ncol=1))
# # MM<-noquote(matrix(mname, ncol=1))
#
# # #bind columns together into one table
# # Table<-noquote(cbind(TT,RR,CC,MM) %>%
# #                  `colnames<-`(c("TLI","RMSEA","CFI","Magnitude")))
#
# #bind columns together into one table
# Table<-noquote(cbind(TT,RR,CC) %>%
#                  `colnames<-`(c("TLI","RMSEA","CFI")))
#
# #update row names
# rname<-c("Level-0","Specificity", "")
# for (i in 1:length(misspec_sum)){
#   ri<-c(paste("Level-",i,sep=""),"Sensitivity","")
#   rname<-cbind(rname,ri)
# }
# rownames(Table)<-rname
#
# #delete last blank row
# Table<-Table[1:(nrow(Table)-1),]
#
# #Put into list
# res$cutoffs <- Table

