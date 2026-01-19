
# simulate data from either the true or the misspecified model

# ----- ggm -----
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
