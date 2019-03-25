# R code for simulating the conditional distribution of the kernel estimator of pdf_X and E(Y|X).
# Code: Andreï V. Kostyrka
# March 18, 2019
# v1.0: 2019-03-18

# Make sure you have ImageMagick installed on your computer, otherwise the animation magic will not happen!
# Go to https://www.imagemagick.org/script/download.php and install the latest version.

rm(list = ls()) # Clear workspace.
source("smoothing-functions.R") # Load the functions.
n <- 200 # Number of observations.
MC <- 10000 # Number of Monte-Carlo simulations
Xgrid <- seq(0, qchisq(0.99, 3), length.out = 101)
fhat <- muhat <- matrix(NA, nrow = length(Xgrid), ncol = MC)
mu <- function(x) 0.2 * x^2 + 2 * sin(x)
mybw <- 0.7
write.img <- TRUE # We shall be using PNG and PDF image formats because PDF cannot handle 10000 lines well

# We are using the following DGP:
# X is chi-squared with 3 df
# U is centred Poisson with lambda 4 (discrete!)
# Y = 0.2*X^2 + 2*sin(X)

# This simulation takes about 20 seconds
system.time({
  for (i in 1:MC) {
    set.seed(i)
    X <- rchisq(n, 3)
    U <- rpois(n, lambda = 4) - 4
    Y <- mu(X) + U
    # plot(X, Y)
    # curve(mu(x), 0, 10, add=TRUE)
    myfhat <- kernelDensity(X, Xgrid, bw = mybw)
    mymuhat <- kernelSmooth(X, Y, Xgrid, bw = mybw)
    fhat[, i] <- myfhat
    muhat[, i] <- mymuhat
    if (i %% 100 == 0) cat(i, "\n")
  }
})

# Substantial productivity gains can be achieved if parallslisation (unavailable on Windows) is used!
library(parallel)
ncores <- if (.Platform$OS.type == "windows") 1 else detectCores()
if (.Platform$OS.type != "windows") {
  system.time({
    a <- mclapply(1:MC, function(i) {
      set.seed(i)
      X <- rchisq(n, 3)
      U <- rpois(n, lambda = 4) - 4
      Y <- mu(X) + U
      myfhat <- kernelDensity(X, Xgrid, bw = mybw)
      mymuhat <- kernelSmooth(X, Y, Xgrid, bw = mybw)
      if (i %% 100 == 0) cat(i, "\n")
      return(list(myfhat, mymuhat))
    }, mc.cores = ncores)
    fhat.parallel <- matrix(unlist(lapply(a, "[[", 1)), ncol = MC)
    muhat.parallel <- matrix(unlist(lapply(a, "[[", 2)), ncol = MC)
  })
  print(all.equal(fhat, fhat.parallel))
  print(all.equal(muhat, muhat.parallel))
}

# Quickly look at the optimal bandwidths
opt.bw.dcv <- optimise(function(b) DCV(X, bw = b), c(0.1, 5))$minimum
opt.bw.lscv <- optimise(function(b) LSCV(X, Y, bw = b), c(0.1, 5))$minimum
cat("Optimal DCV BW is", round(opt.bw.dcv, 3), "and LSCV BW is", round(opt.bw.lscv, 3), ".\n")

# Quickly look at what one density and conditional expectation estimation looked like
if (write.img) pdf("30-estimation-example.pdf", 9, 5)
par(mfrow = c(1, 2))
par(oma = c(0, 0, 0, 0))
par(mar = c(4, 4, 2, 1))
plot(Xgrid, myfhat, ylim = c(0, 0.25), main = "Estimated density in simulation #10000", bty = "n", ylab = "Density", xlab = "x")
curve(dchisq(x, 3), 0, 11, add = TRUE, lwd = 2, lty = 2, col = "blue")
rug(X)
plot(Xgrid, mymuhat, ylim = c(-2, 27), main = "Estimated mu in simulation #10000", bty = "n", ylab = expression(mu(x)), xlab = "x")
curve(0.2 * x^2 + 2 * sin(x), 0, 11, add = TRUE, lwd = 2, lty = 2, col = "red")
rug(X)
points(X, Y, pch = 16, cex = 0.5, col = "red")
dev.off()

# Plot the distribution of the density and regression function estimator
if (write.img) pdf("31-estimator-bands.pdf", 9, 5)
par(mfrow = c(1, 2))
par(oma = c(0, 0, 0, 0))
par(mar = c(4, 4, 2, 1))
plot(Xgrid, dchisq(Xgrid, 3), ylim = c(0, 0.3), col = "blue", ylab = "Density", xlab = "x", type = "l", lty = 2, lwd = 2, bty = "n", main = "Distribution of the density estimator")
lines(Xgrid, apply(fhat, 1, mean), lwd = 2)
for (q in c(0.250, 0.750)) lines(Xgrid, apply(fhat, 1, function(x) quantile(x, q)), lwd = 1, lty = 2)
for (q in c(0.025, 0.975)) lines(Xgrid, apply(fhat, 1, function(x) quantile(x, q)), lwd = 1, lty = 3)
for (q in c(0.005, 0.995)) lines(Xgrid, apply(fhat, 1, function(x) quantile(x, q)), lwd = 1, lty = 4)
legend("topright", c("True density", "Average estimate", "50% (IQ) range", "95% range", "99% range"), lty = c(2, 1, 2, 3, 4), lwd = c(2, 2, 1, 1, 1), col = c("blue", "black", "black", "black", "black"), bty = "n")
plot(Xgrid, mu(Xgrid), ylim = c(-2, 27), col = "red", type = "l", lwd = 2, lty = 2, bty = "n", ylab = expression(mu(x)), xlab = "x", main = "Distribution of the estimator of mu")
lines(Xgrid, apply(muhat, 1, mean), lwd = 2)
for (q in c(0.250, 0.750)) lines(Xgrid, apply(muhat, 1, function(x) quantile(x, q)), lwd = 1, lty = 2)
for (q in c(0.025, 0.975)) lines(Xgrid, apply(muhat, 1, function(x) quantile(x, q)), lwd = 1, lty = 3)
for (q in c(0.005, 0.995)) lines(Xgrid, apply(muhat, 1, function(x) quantile(x, q)), lwd = 1, lty = 4)
legend("topleft", c("True mu(x)", "Average estimate", "50% (IQ) range", "95% range", "99% range"), lty = c(2, 1, 2, 3, 4), lwd = c(2, 2, 1, 1, 1), col = c("red", "black", "black", "black", "black"), bty = "n")
dev.off()

# We can also plot all lines with semi-transparency at once and look at the shades
if (write.img) png("32-estimation-lines.png", 1200, 600, pointsize = 14, type = "cairo")
par(mfrow = c(1, 2))
par(oma = c(0, 0, 0, 0))
par(mar = c(4, 4, 2, 1))
plot(NULL, NULL, xlim = range(Xgrid), ylim = c(0, 0.3), ylab = "Density", xlab = "x", bty = "n", main = "Distribution of the density estimator")
for (l in 1:MC) lines(Xgrid, fhat[, l], lwd = 1, col = "#00000002")
lines(Xgrid, dchisq(Xgrid, 3), col = "blue", lty = 2, lwd = 2)
lines(Xgrid, apply(fhat, 1, mean), lwd = 2, col = "black")
plot(NULL, NULL, xlim = range(Xgrid), ylim = c(-2, 27), bty = "n", ylab = expression(mu(x)), xlab = "x", main = "Distribution of the estimator of mu")
for (l in 1:MC) lines(Xgrid, muhat[, l], lwd = 1, col = "#00000002")
lines(Xgrid, mu(Xgrid), col = "red", lwd = 2, lty = 2)
lines(Xgrid, apply(muhat, 1, mean), lwd = 2, col = "black")
curve(dchisq(x, 3) * 50, 0, 11, add = TRUE, lty = 3, col = "blue")
dev.off()

# And now---film
den.grid <- seq(0, 0.3, length.out = 401)
if (write.img) png("33-density-estimator-distribution.png", 1200, 600, pointsize = 14, type = "cairo")
par(mfrow = c(1, 2))
par(oma = c(0, 0, 0, 0))
par(mar = c(4, 4, 2, 1))
i <- 1
plot(den.grid, kernelDensity(fhat[i, ], den.grid, 0.0004), type = "l", main = paste0("Distribution of KDE at X = ", round(Xgrid[i], 2)), ylab = "Frequency across simulations", xlab = "Density estimator", bty = "n", lwd = 2)
abline(v = dchisq(Xgrid[i], 3), lwd = 2, col = "blue", lty = 2)
rug(fhat[i, ], col = "#00000011")
legend("topright", "True density value", col = "blue", lwd = 2, lty = 2, bty = "n")
i <- 30
plot(den.grid, kernelDensity(fhat[i, ], den.grid, 0.0004), type = "l", main = paste0("Distribution of KDE at X = ", round(Xgrid[i], 2)), ylab = "Frequency across simulations", xlab = "Density estimator", bty = "n", lwd = 2)
abline(v = dchisq(Xgrid[i], 3), lwd = 2, lty = 2, col = "blue")
rug(fhat[i, ], col = "#00000011")
dev.off()

if (!("animation" %in% rownames(installed.packages()))) install.packages("animation")
library(animation) # We can make a wireframe plot using this extra library
options(bitmapType = "cairo")

# Density estimator
saveGIF({ # Might takes from 30 seconds to 5 minutes depenging on the performance of ImageMagick
  for (i in 1:101) {
    plot(den.grid, kernelDensity(fhat[i, ], den.grid, 0.0004), type = "l", main = paste0("Distribution of KDE at X=", round(Xgrid[i], 2)), bty = "n", lwd = 2, ylab = "Frequency across simulations", xlab = "Density")
    abline(v = dchisq(Xgrid[i], 3), lwd = 2, col = "blue", lty = 2)
    rug(fhat[i, ], col = "#00000011")
    legend("topright", legend = c(expression(f[X](x)), expression(hat(f)[X](x))), col = c("blue", "black"), lty = c(2, 1), lwd = 2, bty = "n")
    if (i %% 10 == 0) cat(i, "\n")
  }
}, movie.name = "34-density-estimators-the-film.gif", interval = 0.25, ani.width = 640, ani.height = 480, ani.dev = if (.Platform$OS.type=="windows") function(...) png(..., type="cairo") else "png")

# Regression estimator
mu.grid <- seq(0, 27, length.out = 401)
saveGIF({ # Might takes from 50 seconds to 5 minutes depenging on the performance of ImageMagick
  for (i in 1:101) {
    plot(mu.grid, kernelDensity(muhat[i, ], mu.grid, bw.nrd(muhat[i, ]) / 5), type = "l", main = paste0("Distribution of NW estimator at X =", round(Xgrid[i], 2)), lwd = 2, ylab = "Frequency across simulations", xlab = expression(mu))
    abline(v = mu(Xgrid[i]), lwd = 2, col = "red")
    rug(muhat[i, ], col = "#00000011")
    legend("topright", legend = expression(mu(x)), col = "red", lwd = 2, bty = "n")
    if (i %% 10 == 0) cat(i, "\n")
  }
}, movie.name = "35-regression-estimators-the-film.gif", interval = 0.25, ani.width = 640, ani.height = 480, ani.dev = if (.Platform$OS.type=="windows") function(...) png(..., type="cairo") else "png")

# In 3D, we can use a coarser grig, but generating 101 densities still takes a while, so we parallelise while we can
den.grid <- seq(0, 0.3, length.out = 201)
kdes <- mclapply(1:101, function(i) kernelDensity(fhat[i, ], den.grid, bw.nrd(fhat[i, ]) / 5), mc.cores = ncores)
kdes <- lapply(kdes, function(x) x / max(x)) # Normalising them to one for nice plotting
saveGIF({ # Might takes from 20 seconds to 5 minutes depenging on the performance of ImageMagick
  for (a in 1:360) {
    par(mar = c(2, 2, 2, 2))
    p <- persp(c(0, 11), c(0, 0.3), matrix(rep(0, 4), 2), col = "white", zlim = c(0, 1.01), theta = a, phi = 40, xlab = "x", ylab = "Density", zlab = "KDE concentration")
    lines(trans3d(Xgrid, dchisq(Xgrid, 3), rep(0, 101), pmat = p), lwd = 1, col = "blue")
    for (i in 1:101) lines(trans3d(rep(Xgrid[i], length(den.grid)), den.grid, kdes[[i]], pmat = p), lwd = 1, col = "#00000077")
    if (a %% 5 == 0) cat(a, "/ 360\n")
  }
}, movie.name = "36-kde-3d.gif", interval = 1 / 25, ani.width = 480, ani.height = 480, ani.dev = if (.Platform$OS.type=="windows") function(...) png(..., type="cairo") else "png")

mu.grid <- seq(0, 27, length.out = 201)
mus <- mclapply(1:101, function(i) kernelDensity(muhat[i, ], mu.grid, bw.nrd(muhat[i, ]) / 5), mc.cores = ncores)
mus <- lapply(mus, function(x) x / max(x))
saveGIF({ # Might takes from 20 seconds to 5 minutes depenging on the performance of ImageMagick
  for (a in 1:360) {
    par(mar = c(2, 2, 2, 2))
    p <- persp(c(0, 11), c(0, 25), matrix(rep(0, 4), 2), col = "white", zlim = c(0, 1.01), theta = a, phi = 40, xlab = "x", ylab = "mu", zlab = "Estimator density")
    lines(trans3d(Xgrid, mu(Xgrid), rep(0, 101), pmat = p), lwd = 2, col = "red")
    for (i in 1:101) lines(trans3d(rep(Xgrid[i], length(mu.grid)), mu.grid, mus[[i]], pmat = p), lwd = 1, col = "#00000077")
    if (a %% 5 == 0) cat(a, "/ 360\n")
  }
}, movie.name = "37-mu-3d.gif", interval = 1 / 25, ani.width = 480, ani.height = 480, ani.dev = if (.Platform$OS.type=="windows") function(...) png(..., type="cairo") else "png")
