# Test hyperparameter configurations
library(SDMtune)

# Simulate different sample sizes
test_sizes <- c(92, 252, 300, 400, 600)

for (n in test_sizes) {
  cat("\n=== Sample size:", n, "===\n")

  if (n < 300) {
    h <- list(
      reg = c(1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10),
      fc = c("lq", "lh", "qh", "lqh", "lp", "qp")
    )
    cat("Configuration: CONSERVATIVE\n")
  } else if (n < 500) {
    h <- list(
      reg = seq(0.5, 10, 0.5),
      fc = c("lq", "lp", "lh", "qp", "qh", "lqp", "lqh", "qph")
    )
    cat("Configuration: MODERATE\n")
  } else {
    h <- list(
      reg = seq(0.5, 10, 0.5),
      fc = c("lq", "lp", "lh", "qp", "qh", "ph", "lqp", "lqh", "lph", "qph", "lqph")
    )
    cat("Configuration: FULL\n")
  }

  cat("Regularization values:", length(h$reg), "\n")
  cat("Feature classes:", length(h$fc), "\n")
  cat("Total combinations:", length(h$reg) * length(h$fc), "\n")
  cat("Feature classes:", paste(h$fc, collapse = ", "), "\n")
}
