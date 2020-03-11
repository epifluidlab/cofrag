library(testthat)

source("../calc_distance.R", chdir = TRUE)

test_that("Have a test", {
  expect_equal(1, 1)
  expect_equal(2, 2)
})

test_that("Test calc_bfp", {
  # Prepare the simulation data
  n <- 10000
  bin_size <- 100
  gr_start <- 16287
  gr_end <- gr_start + n - 1
  gr <- GRanges(sprintf("22:%d-%d", gr_start, gr_end))

  min_isize <- -112
  max_isize <- 187
  aligned_reads <-
    list(pos = sort(round(runif(
      n, gr_start, gr_end
    ))), isize = round(runif(n, min_isize, max_isize)))

  # Calculated value
  val <- calc_bfp(aligned_reads, gr, bin_size)

  # Expected value
  bin_num <- n / bin_size
  ex_val <- lapply(1:bin_num, function(bin_idx) {
    bin_start <- gr_start + (bin_idx - 1) * bin_size
    bin_end <- gr_start + bin_idx * bin_size - 1
    abs(aligned_reads$isize[aligned_reads$pos >= bin_start &
                              aligned_reads$pos <= bin_end])
  })

  expect_true(identical(ex_val, val))
  # expect_true(TRUE)
})