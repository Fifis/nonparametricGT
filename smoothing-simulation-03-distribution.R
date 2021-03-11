# R code for simulating the conditional distribution of the kernel estimator of pdf_X and E(Y|X).
# Code: Andreï V. Kostyrka
# March 18, 2019
# v1.0: 2019-03-18
# v1.2: 2021-03-10 (improved the animation, added the HD film)

# Make sure you have ImageMagick installed on your computer, otherwise the animation magic will not happen!
# Go to https://www.imagemagick.org/script/download.php and install the latest version.

rm(list = ls()) # Clear the workspace.
source("nonparametricGT/R/smoothing-functions.R") # Load the functions.
n <- 200 # Number of observations.
MC <- 10000 # Number of Monte-Carlo simulations
Xgrid <- unique(sort(c(seq(0, 10, 0.1), seq(10, 12, 0.05))))
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

# Substantial productivity gains can be achieved if parallelisation (unavailable on Windows) is used!
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
opt.bw.dcv <- bw.CV(x = X, CV = "DCV")
opt.bw.lscv <- bw.CV(x = X, y = Y, CV = "LSCV")
cat("The optimal BW_DCV is ", sprintf("%1.3f", opt.bw.dcv), " and BW_LSCV is ", sprintf("%1.3f", opt.bw.lscv), ".\n", sep = "")

# Quickly look at what one density and conditional expectation estimation looked like
if (write.img) pdf("30-estimation-example.pdf", 9, 5)
par(mfrow = c(1, 2))
par(oma = c(0, 0, 0, 0))
par(mar = c(4, 4, 2, 1))
plot(Xgrid, myfhat, ylim = c(0, 0.25), main = "Estimated density in simulation #10000", bty = "n", ylab = "Density", xlab = "x")
curve(dchisq(x, 3), 0, max(Xgrid), add = TRUE, lwd = 2, lty = 2, col = "blue")
rug(X)
plot(Xgrid, mymuhat, ylim = c(-2, 27), main = "Estimated mu in simulation #10000", bty = "n", ylab = expression(mu(x)), xlab = "x")
curve(mu(x), 0, max(Xgrid), add = TRUE, lwd = 2, lty = 2, col = "red")
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
for (q in c(0.250, 0.750)) lines(Xgrid, apply(fhat, 1, function(x) quantile(x, q)), lwd = 1.5, lty = 3, col = "#000000FF")
for (q in c(0.025, 0.975)) lines(Xgrid, apply(fhat, 1, function(x) quantile(x, q)), lwd = 1.5, lty = 3, col = "#00000099")
for (q in c(0.005, 0.995)) lines(Xgrid, apply(fhat, 1, function(x) quantile(x, q)), lwd = 1.5, lty = 3, col = "#00000066")
legend("topright", c("True density", "Average estimate", "50% (IQ) range", "95% range", "99% range"), lty = c(2, 1, 3, 3, 3), lwd = c(2, 2, 1.5, 1.5, 1.5), col = c("blue", "black", "#000000FF", "#00000099", "#00000066"), bty = "n")
plot(Xgrid, mu(Xgrid), ylim = c(-2, 27), col = "red", type = "l", lwd = 2, lty = 2, bty = "n", ylab = expression(mu(x)), xlab = "x", main = "Distribution of the estimator of mu")
lines(Xgrid, apply(muhat, 1, mean), lwd = 2)
for (q in c(0.250, 0.750)) lines(Xgrid, apply(muhat, 1, function(x) quantile(x, q)), lwd = 1.5, lty = 3, col = "#000000FF")
for (q in c(0.025, 0.975)) lines(Xgrid, apply(muhat, 1, function(x) quantile(x, q)), lwd = 1.5, lty = 3, col = "#00000099")
for (q in c(0.005, 0.995)) lines(Xgrid, apply(muhat, 1, function(x) quantile(x, q)), lwd = 1.5, lty = 3, col = "#00000066")
legend("topleft", c(expression("True" ~ mu(x)), "Average estimate", "50% (IQ) range", "95% range", "99% range"), lty = c(2, 1, 3, 3, 3), lwd = c(2, 2, 1.5, 1.5, 1.5), col = c("red", "black", "#000000FF", "#00000099", "#00000066"), bty = "n")
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
den.grid <- seq(0, 0.3, length.out = 201)
if (write.img) png("33-density-estimator-distribution.png", 1200, 600, pointsize = 14, type = "cairo")
par(mfrow = c(1, 2))
par(oma = c(0, 0, 0, 0))
par(mar = c(4, 4, 2, 1))
i <- 2
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
  for (i in 1:length(Xgrid)) {
    kd <- kernelDensity(fhat[i, ], den.grid, bw = bw.rot(fhat[i, ]) / 3)
    plot(NULL, NULL, main = paste0("Distribution of the KDE at X = ", sprintf("%1.2f", Xgrid[i])), bty = "n", lwd = 2, ylab = "Density in simulations", xlab = "x", xlim = range(Xgrid), ylim = c(0, max(fhat)))
    lines(kd / max(X) * 3, den.grid, lwd = 2)
    rug(fhat[i, ], side = 2, col = "#00000011")
    abline(h = quantile(fhat[i, ], c(0.025, 0.975)), lty = 3)
    legend("topright", legend = c(expression(f[X](x)), expression(hat(f)[X](x) ~ "density"), "Point-wise 95% CI"), col = c("blue", "black", "black"), lty = c(2, 1, 3), lwd = c(2, 2, 1), bty = "n")
    lines(Xgrid, dchisq(Xgrid, 3), col = "blue", lty = 2, lwd = 2)
    if (i > 1) for (q in c(0.025, 0.975)) lines(Xgrid[1:i], apply(fhat[1:i, ], 1, function(x) quantile(x, q)), lty = 3, col = "#000000")
    abline(h = dchisq(Xgrid[i], 3), lwd = 2, col = "blue")
    abline(v = Xgrid[i], lwd = 1, col = "blue")
    points(Xgrid[i], dchisq(Xgrid[i], 3), pch = 16, col = "blue")
    if (i %% 10 == 0) cat(i, "\n")
  }
}, movie.name = "34-density-estimators-the-film.gif", interval = 0.25, ani.width = 640, ani.height = 480, ani.dev = if (.Platform$OS.type=="windows") function(...) png(..., type = "cairo") else "png")

# Regression estimator
max.mu <- quantile(muhat, 0.995)
mu.grid <- seq(0, max.mu, length.out = 201)
saveGIF({ # Might takes from 50 seconds to 5 minutes depenging on the performance of ImageMagick
  for (i in 1:length(Xgrid)) {
    md <- kernelDensity(muhat[i, ], mu.grid, bw = bw.rot(muhat[i, ]) / 3)
    plot(NULL, NULL, main = paste0("Distribution of the Nadaraya—Watson estimator at X = ", sprintf("%1.2f", Xgrid[i])), bty = "n", lwd = 2, ylab = "Density in simulations", xlab = "x", xlim = range(Xgrid), ylim = c(0, max.mu))
    lines(md / max(X) * 80, mu.grid, lwd = 2)
    lines(Xgrid, dchisq(Xgrid, df = 3) * max.mu * 4, lty = 2, col = "#0000FF88")
    rug(muhat[i, ], side = 2, col = "#00000011")
    abline(h = quantile(muhat[i, ], c(0.025, 0.975)), lty = 3)
    legend("topright", legend = c(expression(mu(x)), expression(hat(mu)(x) ~ "density"), "Point-wise 95% CI", expression(f[X](x))), col = c("red", "black", "black", "#0000FF88"), lty = c(2, 1, 3, 2), lwd = c(2, 2, 1, 1), bty = "n")
    lines(Xgrid, mu(Xgrid), col = "red", lty = 2, lwd = 2)
    if (i > 1) for (q in c(0.025, 0.975)) lines(Xgrid[1:i], apply(muhat[1:i, ], 1, function(x) quantile(x, q)), lty = 3, col = "#000000")
    abline(h = mu(Xgrid[i]), lwd = 2, col = "red")
    abline(v = Xgrid[i], lwd = 1, col = "red")
    points(Xgrid[i], mu(Xgrid[i]), pch = 16, col = "red")
    if (i %% 10 == 0) cat(i, "\n")
  }
}, movie.name = "35-regression-estimators-the-film.gif", interval = 0.25, ani.width = 640, ani.height = 480, ani.dev = if (.Platform$OS.type=="windows") function(...) png(..., type="cairo") else "png")

X.inds <- which(round(Xgrid, 2) %in% round(seq(0, max(Xgrid), 0.1), 2)) # For uniform plotting

# In 3D, we can use a coarser grig, but generating more than 100 densities still takes a while, so we parallelise while we can
kdes <- mclapply(1:length(Xgrid), function(i) kernelDensity(fhat[i, ], den.grid, bw.rot(fhat[i, ]) / 2), mc.cores = ncores)
kde.peaks <- unlist(lapply(kdes, max))
kde.peaks <- (kde.peaks - min(kde.peaks)) / diff(range(kde.peaks)) / 2 + 0.5
# We want to make it visually appealing by making the vertical range of peaks fixed between 0.5 and 1
kdes <- lapply(1:length(kdes), function(i) kdes[[i]] / max(kdes[[i]]) * kde.peaks[[i]])
saveGIF({ # Might takes from 20 seconds to 5 minutes depending on the performance of ImageMagick
  for (a in seq(0, 359, 2)) {
    par(mar = c(2, 2, 2, 2))
    p <- persp(range(Xgrid), range(den.grid), matrix(rep(0, 4), 2), col = "white", zlim = c(0, 1.01), theta = a, phi = 40, xlab = "x", ylab = "Density", zlab = "KD estimator density")
    lines(trans3d(Xgrid, dchisq(Xgrid, 3), rep(0, length(Xgrid)), pmat = p), lwd = 1, col = "blue")
    for (i in X.inds) lines(trans3d(rep(Xgrid[i], length(den.grid)), den.grid, kdes[[i]], pmat = p), lwd = 1, col = "#00000077")
    if (a %% 5 == 0) cat(a, "/", 360, "\n")
  }
}, movie.name = "36-kde-3d.gif", interval = 1 / 10, ani.width = 480, ani.height = 480, ani.dev = if (.Platform$OS.type=="windows") function(...) png(..., type="cairo") else "png")

mus <- mclapply(1:length(Xgrid), function(i) kernelDensity(muhat[i, ], mu.grid, bw.nrd(muhat[i, ]) / 2), mc.cores = ncores)
mu.peaks <- unlist(lapply(mus, max))
mu.peaks <- (mu.peaks - min(mu.peaks)) / diff(range(mu.peaks)) / 2 + 0.5
mus <- lapply(1:length(mus), function(i) mus[[i]] / max(mus[[i]]) * mu.peaks[[i]])
saveGIF({ # Might takes from 20 seconds to 5 minutes depenging on the performance of ImageMagick
  for (a in seq(0, 359, 2)) {
    par(mar = c(2, 2, 2, 2))
    p <- persp(range(Xgrid), range(mu.grid), matrix(rep(0, 4), 2), col = "white", zlim = c(0, 1.01), theta = a, phi = 40, xlab = "x", ylab = "mu", zlab = "NW estimator density")
    lines(trans3d(Xgrid, mu(Xgrid), rep(0, length(Xgrid)), pmat = p), lwd = 2, col = "red")
    for (i in X.inds) lines(trans3d(rep(Xgrid[i], length(mu.grid)), mu.grid, mus[[i]], pmat = p), lwd = 1, col = "#00000077")
    if (a %% 5 == 0) cat(a, "/ 360\n")
  }
}, movie.name = "37-mu-3d.gif", interval = 1 / 10, ani.width = 480, ani.height = 480, ani.dev = if (.Platform$OS.type=="windows") function(...) png(..., type="cairo") else "png")

# Exporting to a better and smoother format with FFmpeg
j <- 1000
for (a in seq(0, 359.5, 0.5)) {
  j <- j + 1
  png(file = paste0("/tmp/example", j, ".png"), width = 960, height = 720, type = "cairo", pointsize = 18)
  par(mar = c(2, 2, 2, 2))
  p <- persp(range(Xgrid), range(den.grid), matrix(rep(0, 4), 2), col = "white", zlim = c(0, 1.01), theta = a, phi = 40, xlab = "x", ylab = "Density", zlab = "KD estimator density")
  lines(trans3d(Xgrid, dchisq(Xgrid, 3), rep(0, length(Xgrid)), pmat = p), lwd = 1, col = "blue")
  for (i in X.inds) lines(trans3d(rep(Xgrid[i], length(den.grid)), den.grid, kdes[[i]], pmat = p), lwd = 1.5, col = "#00000077")
  dev.off()
  if (a %% 10 == 0) cat(a, "/ 360\n")
}
system("ffmpeg -y -framerate 30 -pattern_type glob -i '/tmp/example*.png' -c:v libx265 -crf 29 -preset slower -tune animation -pix_fmt yuv420p ./38-kde-3d.mp4")

j <- 1000
for (a in seq(0, 359.5, 0.5)) {
  j <- j + 1
  png(file = paste0("/tmp/mu", j, ".png"), width = 960, height = 720, type = "cairo", pointsize = 18)
  par(mar = c(2, 2, 2, 2))
  p <- persp(range(Xgrid), range(mu.grid), matrix(rep(0, 4), 2), col = "white", zlim = c(0, 1.01), theta = a, phi = 40, xlab = "x", ylab = "Density", zlab = "NW estimator density")
  lines(trans3d(Xgrid, mu(Xgrid), rep(0, length(Xgrid)), pmat = p), lwd = 1, col = "red")
  for (i in X.inds) lines(trans3d(rep(Xgrid[i], length(mu.grid)), mu.grid, mus[[i]], pmat = p), lwd = 1.5, col = "#00000077")
  dev.off()
  if (a %% 10 == 0) cat(a, "/ 360\n")
}
system("ffmpeg -y -framerate 30 -pattern_type glob -i '/tmp/mu*.png' -c:v libx265 -crf 29 -preset slower -tune animation -pix_fmt yuv420p ./39-mu-3d.mp4")
file.remove(list.files("/tmp/", pattern=".png", full.names = TRUE))

print("Done!")

