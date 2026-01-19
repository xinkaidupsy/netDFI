

# ----- a fit function that does not print output -----
silent_fit <- function(model) {
  # Capture and discard printed output, return the data frame
  utils::capture.output(
    fit_result <- psychonetrics::fit(model),
    file = NULL
  )
  return(fit_result)
}

# ----- net2vec -----
# obtain the lower.tri of a network, label each edge with the nodes that it connects, &
# output the vector of edges annotated with corresponding nodes
net2vec <- function(net){

  # obtain lower.tri
  vec <- data.frame(net[lower.tri(net)])
  names(vec) <- "weights"

  # annotate values with row&colnames
  names.Mat <- outer(rownames(net), colnames(net), paste, sep = "--")
  vec$loc <- names.Mat[lower.tri(names.Mat)]

  return(vec)

}

# merge true and misspecified fit dataframes into one dataframe with appropriate names
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
