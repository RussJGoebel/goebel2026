par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

for (surf_name in c("albedo", "sif")) {
  surf <- if (surf_name == "albedo") likelihood_surface_albedo else likelihood_surface_sif
  lp   <- surf$log_phi_seq

  # --- ML curve (normalized) ---
  ml  <- surf$ml_ll
  ml  <- (ml - max(ml, na.rm=TRUE)) / abs(diff(range(ml, na.rm=TRUE)))

  # --- CV curves (normalized, flipped so higher=better) ---
  cvr <- surf$cv_mse_random
  cvr <- -(cvr - min(cvr, na.rm=TRUE)) / abs(diff(range(cvr, na.rm=TRUE)))

  cvb <- surf$cv_mse_blocked
  cvb <- -(cvb - min(cvb, na.rm=TRUE)) / abs(diff(range(cvb, na.rm=TRUE)))

  # Panel 1: ML
  plot(lp, ml, type="l", col="tomato", lwd=2,
       ylim=c(-1,0), xlab="log phi", ylab="Scaled objective",
       main=sprintf("%s: ML", surf_name))
  abline(v=log(surf$ml_phi_opt), col="tomato", lty=2)
  text(log(surf$ml_phi_opt), -0.9,
       sprintf("phi=%.3g", surf$ml_phi_opt),
       col="tomato", adj=c(-0.1,0), cex=0.8)

  # Panel 2: CV random
  plot(lp, cvr, type="l", col="steelblue", lwd=2,
       ylim=c(-1,0), xlab="log phi", ylab="Scaled objective",
       main=sprintf("%s: CV random", surf_name))
  abline(v=log(surf$cv_phi_random), col="steelblue", lty=2)
  text(log(surf$cv_phi_random), -0.9,
       sprintf("phi=%.3g", surf$cv_phi_random),
       col="steelblue", adj=c(-0.1,0), cex=0.8)

  # Panel 3: CV blocked
  plot(lp, cvb, type="l", col="darkgreen", lwd=2,
       ylim=c(-1,0), xlab="log phi", ylab="Scaled objective",
       main=sprintf("%s: CV blocked", surf_name))
  abline(v=log(surf$cv_phi_blocked), col="darkgreen", lty=2)
  text(log(surf$cv_phi_blocked), -0.9,
       sprintf("phi=%.3g", surf$cv_phi_blocked),
       col="darkgreen", adj=c(-0.1,0), cex=0.8)
}
