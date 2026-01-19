
#' Plot functions
#'
#' Plot the fit distribution
#'
#' @param x an object of class dfi_ggm
#' @param ... further arguments passed from plot(); currently ignored.
#' @return Plots
#'
#' @export
#'
#'
plot.dfi_ggm <- function(x, ...) {

  # Check if object has the right class
  if (!inherits(x, "dfi_ggm")) {
    stop("Object must be of class 'dfi_ggm'")
  }

  # obtain L0 & Ls data
  L0 <- x$fit %>%
    select(matches("L0")) %>%
    `colnames<-`(c("TLI", "RMSEA", "CFI")) %>%
    `rownames<-`(NULL) %>%
    mutate(Model = "true")

  Ls <- lapply(seq_len(nrow(x$cutoff_misspec)), function(i){
    x$fit %>%
      select(matches(paste0("L",i))) %>%
      `colnames<-`(c("TLI", "RMSEA", "CFI")) %>%
      `rownames<-`(NULL) %>%
      mutate(Model = "misspecified")
  })

  # merge to create plot data
  plot_dt <- lapply(Ls, function(x){
    rbind(x, L0)
  })

  # plot
  p_patched <- lapply(seq_along(plot_dt), function(i) {

    p <- lapply(c("TLI", "RMSEA", "CFI"), function(index) {
      ggplot2::ggplot(plot_dt[[i]], aes(x = get(index), fill = Model)) +
        geom_histogram(position = "identity", bins=30, alpha = 0.5) +
        scale_fill_brewer(palette = "Set1") +
        geom_vline(aes(xintercept = x[["cutoff_misspec"]][[index]][i],
                       linetype = "Dynamic cutoff")) +
        geom_vline(aes(xintercept = x[["cutoff_misspec"]][[paste0(index,"_HB")]][i],
                       linetype = "Hu & Bentler cutoff")) +
        scale_linetype_manual(values = c("Dynamic cutoff" = "longdash",
                                         "Hu & Bentler cutoff" = "dotted")) +
        labs(y = "", x = index, linetype = "Cutoffs") +
        theme_classic() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }) %>% setNames(c("TLI", "RMSEA", "CFI"))

    patchwork::wrap_plots(p) +
      plot_layout(guides = "collect") +
      plot_annotation(paste("Level", i)) &
      theme(legend.position = 'bottom')

  })

  return(p_patched)

}
