# Fast end-to-end smoke test for dfi_ggm() on a minimal 3-node GGM.
# Designed to run in well under 5 seconds (small n_misspec / iter / n).

test_that("dfi_ggm runs end-to-end on a small 3-node network", {
  skip_on_cran()

  # minimal 3-node partial-correlation network: a single V1--V2 edge,
  # leaving two zero-edges as mis-specification candidates.
  net <- matrix(0, 3, 3)
  net[1, 2] <- net[2, 1] <- 0.3

  set.seed(1)
  res <- dfi_ggm(
    net,
    n_misspec   = 1,
    iter        = 75,
    n           = 200,
    ncores      = 1,
    progressbar = FALSE
  )

  # class + core structure
  expect_s3_class(res, "dfi_ggm")
  expect_named(
    res,
    c("modified_edges", "fit", "mod_misspec", "cutoff_true", "cutoff_misspec")
  )

  # one mis-spec level was produced
  expect_length(res$mod_misspec, 1L)
  expect_equal(nrow(res$cutoff_misspec), 1L)

  # L0 (true) cutoff row present and labelled
  expect_equal(res$cutoff_true$mis_level, "L0")

  # cutoff table exposes the four fit indices
  expect_true(all(c("TLI", "RMSEA", "CFI", "SRMR") %in% names(res$cutoff_misspec)))
})
