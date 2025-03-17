#' The invariance partial pruning (IVPP) algorithm for panel GVAR models
#'
#' This function implements the IVPP algorithm to compare networks in the multi-group panelGVAR models.
#' The IVPP algorithm is a two-step procedure that first conducts an global invariance test of network difference and then performs partial pruning for the specific edge-level differences.
#' Currently supports the comparison of temporal and contemporaneous networks.
#'
#' @param data A data frame containing the long-formatted panel data
#' @param vars A character vector of variable names
#' @param idvar A character string specifying subject IDs
#' @param beepvar A character string specifying the name of wave (time) variable
#' @param groups A character string specifying the name of group variable
#' @param test A character vector specifying the network you want to test group-equality on in the global invariance test.
#' Specify "both" if you want to test on both temporal or contemporaneous networks.
#' Specify "temporal" if you want to test only on the temporal network.
#' Specify "contemporaneous" if you want to test only on the contemporaneous network.
#' See the Details section for more information.
#' @param net_type A character vector specifying the type of networks to be compared in the global test.
#' Specify "saturated" if you want to estimate and compare the saturated networks.
#' Specify "sparse" if you want to estimate and compare the pruned networks.
#' @param prune_alpha A numeric value specifying the alpha level for the pruning (if net_type = "sparse").
#' @param partial_prune A logical value specifying whether to conduct partial pruning test or not.
#' @param prune_net A character vector specifying the network you want to partial prune on. Only works when partial_prune = TRUE.
#' @param p_prune_alpha A numeric value specifying the alpha level for the partial pruning (if partial_prune = TRUE).
#' @param estimator A character string specifying the estimator to be used. Must be "FIML"
#' @param standardize A character string specifying the type of standardization to be used.
#' "none" (default) for no standardization, "z" for z-scores,
#' and "quantile" for a non-parametric transformation to the quantiles of the marginal standard normal distribution.
#' @param ... Additional arguments to be passed to the \code{\link[psychonetrics]{dlvm1}} function
#'
#' @details
#' The comparison between the fully unconstrained (free) model and tempEq model is a test for group equality in temporal networks.
#' The comparison between fully constrained model (bothEq) and tempEq is a test for group equality in contemporaneous networks.
#' Similarly, the comparison between the free model and contEq model is a test for group equality in contemporaneous networks,
#' and the comparison between bothEq and contEq is a test for group equality in temporal networks.
#'
#' @return A list containing the results of IVPP and networks of all groups.
#'
#' @import dplyr
#' @import psychonetrics
#' @export ggm_dfi
#' @examples
