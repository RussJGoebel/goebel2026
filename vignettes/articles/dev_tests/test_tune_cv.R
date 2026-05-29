# ------------------------------------------------------------
# tiny helper
# ------------------------------------------------------------
time_it <- function(expr, label = deparse(substitute(expr))) {
  gc()
  t <- system.time(val <- eval.parent(substitute(expr)))
  cat(sprintf("%-35s  %.2f sec\n", label, t["elapsed"]))
  invisible(list(time = t, value = val))
}

score_fn <- fastblm:::.make_score_fn("mse")
prior0   <- Q_fun_soft(c(rho = 0.9, log_alpha = 0))
phi0     <- 82
fold0    <- 1L

train_idx0 <- which(fold_assignments != fold0)
test_idx0  <- which(fold_assignments == fold0)

A_train0 <- A_aug[train_idx0, , drop = FALSE]
y_train0 <- y_alb[train_idx0]
A_test0  <- A_aug[test_idx0, , drop = FALSE]
y_test0  <- y_alb[test_idx0]
