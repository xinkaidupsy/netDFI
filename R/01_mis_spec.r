# create misspecification by
# 1. identifying the candidate location
# 2. selecting edges based on priority & non-replication
# 3. adding the edges & checking positive definiteness

# ----- edge_df -----
# create a full list of locations that don't have an edge yet
# summed predictability of neighboring nodes that these locations connect are included
edge_df <- function(net){

  # calculate predictability of each edge
  r2 <- qgraph::centrality(net, R2 = TRUE)$R2

  # find the nodes with highest r2
  r2_df_from <- data.frame(
    r2_from = r2,
    from = as.numeric(gsub("V(\\d+)", "\\1", colnames(net)))
  ) %>%
    arrange(desc(r2_from)) %>%
    mutate(priority = 1:nrow(.))

  # calculate the summed predictability of nodes connected by each edge
  # this is to find zero-edges
  r2_sum <- matrix(r2, nrow = nrow(net), ncol = ncol(net)) +
    matrix(r2, nrow = nrow(net), ncol = ncol(net), byrow = TRUE)

  # fix non-zero & diagnal positions to be 2 so that only zero-edges remain
  r2_sum[net!=0] <- 2
  diag(r2_sum) <- 2

  rownames(r2_sum) <- colnames(r2_sum) <- rownames(net)

  # find positions of zero edges
  edge_df_1 <- net2vec(r2_sum) %>%
    mutate(from = as.numeric(gsub("V\\d+--V(\\d+)", "\\1", loc)),
           to = as.numeric(gsub("V(\\d+)--V\\d+", "\\1", loc)))

  edge_df_2 <- net2vec(r2_sum) %>%
    mutate(to = as.numeric(gsub("V\\d+--V(\\d+)", "\\1", loc)),
           from = as.numeric(gsub("V(\\d+)--V\\d+", "\\1", loc)))

  edge_df <- rbind(edge_df_1, edge_df_2) %>%
    filter(weights != 2)

  # join r2 df and edge df, order by priority & weights
  # such that the search would start from nodes of highest r2 & try to connect
  # to nodes of lowest r2 that the node of highest r2 has yet connected to
  edge_df_final <- inner_join(r2_df_from, edge_df, by = "from") %>%
    arrange(priority, weights) %>%
    distinct(loc, .keep_all = TRUE)

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
# add the additional parameter to selected edge locations &
# return if positive definite
ggm_add <- function(net, edge_cand_ls, n_misspec, prop_pos, n, size_extra, manual_size) {

  # number of nodes
  p <- nrow(net)

  # for the beta_min criterion, compute the node-specific minimum detectable
  # partial correlation (Buhlmann & Van De Geer, 2011): sqrt(log(p) / n) / Theta_ii,
  # where Theta_ii is the diagonal of the precision matrix. The package only receives
  # the network (omega), so we reconstruct the model-implied correlation matrix
  # (K proportional to I - omega) and take diag(solve(R)) as Theta_ii.
  if (size_extra == "beta_min") {
    R_implied <- stats::cov2cor(solve(diag(p) - net))
    conditional_precision <- diag(solve(R_implied))          # Theta_ii, node-specific
    beta_min_vec <- sqrt(log(p) / n) / conditional_precision # length p, indexed V1..Vp
  }

  # Create a list to store resulting networks
  result_networks <- list()

  # Try each starting point
  for (start_idx in 1:length(edge_cand_ls)) {

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

      # Determine the magnitude of each added edge, then apply a random sign.
      # "manual": a fixed magnitude for every edge.
      # "beta_min": node-specific magnitude sized by the edge's "from" node.
      edge_mag <- if (size_extra == "manual") {
        rep(manual_size, nrow(edges_to_add[[k]]))
      } else {
        beta_min_vec[edges_to_add[[k]]$from]
      }
      edge_sign <- sample(c(1, -1),
                          nrow(edges_to_add[[k]]),
                          prob = c(prop_pos, 1 - prop_pos),
                          replace = TRUE)
      edges_to_add[[k]]$added_edges <- edge_mag * edge_sign

      # Add each edge with the determined weight
      for (i in 1:nrow(edges_to_add[[k]])) {

        from_node <- edges_to_add[[k]]$from[i]
        to_node <- edges_to_add[[k]]$to[i]

        # Add the edge in both directions (symmetric matrix)
        mod_misspec[from_node, to_node] <- mod_misspec[to_node, from_node] <- edges_to_add[[k]]$added_edges[i]

      }

      # Check the positive definiteness of I - omega
      eigen_vals <- eigen((diag(nrow(mod_misspec)) - mod_misspec), symmetric = TRUE)$values

      # Flag to check if all networks for this starting point are positive definite
      all_positive_definite <- TRUE

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
