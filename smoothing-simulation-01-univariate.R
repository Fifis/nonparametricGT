# R code for kernel estimation of pdf_X and E(Y|X).
# Original code and idea: Gautam Tripathi
# Rewrite: Andreï V. Kostyrka
# v1.0: 2019-03-07
# v1.1: 2019-03-18 (made plotting optional)
# v1.2: 2021-03-10 (implemented a convenient interface for data-driven bandwidth selection)

# This file is expected to run in 2 seconds for n = 100, in 6 s for n = 200, in 40 s for n = 400, and 180 s for n = 1000
rm(list = ls()) # Clear the workspace.

devtools::install_github("Fifis/nonparametricGT", subdir = "package")
library(nonparametricGT)
# In case installation from GitHub fails, try downloading the package manually and running:
# source("nonparametricGT/package/R/smoothing-functions.R")
start.time <- Sys.time() # Start clock.
set.seed(12345678) # Set seed for replication.
write.pdf <- FALSE # Do we want plots to be shown on the screen, or to be written in PDF files?

n <- 100 # Number of observations.

all.kernels <- c("uniform", "triangular", "epanechnikov", "quartic", "gaussian")
my.colours <- c("#000000CC", "#0000CCCC", "#CC0000CC", "#00AA00CC", "#BB8800CC")

# Visualise the kernels. We want to save the plots in PDF.
if (write.pdf) pdf(file = "10-kernels.pdf", width = 7, height = 5) # Open PDF file for writing.
curve(kernelFun(x, kernel = "uniform", rescale = FALSE), -2, 2, n = 501, ylim = c(0, 1.1), col = my.colours[1], lwd = 2, ylab = "Kernel", main = "Various kernels used in smoothing", bty = "n")
for (i in 2:5) curve(kernelFun(x, kernel = all.kernels[i], rescale = FALSE), -2, 2, add = TRUE, col = my.colours[i], lwd = 2)
legend("topright", legend = all.kernels, lwd = 2, col = my.colours, bty = "n")
dev.off() # Close the graphical device (if a PDF is being written, then finalises it).

# Visualise the kernels so that the integral of x^2 k(x) dx over R be 1.
if (write.pdf) pdf(file = "11-kernels-rescaled.pdf", width = 7, height = 5)
curve(kernelFun(x, kernel = "uniform", rescale = TRUE), -3, 3, n = 501, ylim = c(0, 1.1), col = my.colours[1], lwd = 2, ylab = "Rescaled kernel", bty = "n", main = "Rescaled kernels")
for (i in 2:5) curve(kernelFun(x, kernel = all.kernels[i], rescale = TRUE), -3, 3, add = TRUE, col = my.colours[i], lwd = 2)
legend("topright", legend = all.kernels, lwd = 2, col = my.colours, bty = "n")
dev.off()

for (k in all.kernels) print(integrate(function(x) kernelFun(x, kernel = k, rescale = FALSE), -Inf, Inf, rel.tol = 1e-8)) # Integrates to one.
for (k in all.kernels) print(integrate(function(x) kernelFun(x, kernel = k, rescale = TRUE), -Inf, Inf, rel.tol = 1e-8)) # Integrates to one.
for (k in all.kernels) print(integrate(function(x) x^2 * kernelFun(x, kernel = k, rescale = FALSE), -Inf, Inf, rel.tol = 1e-8)) # Depends on the shape of the kernel.
for (k in all.kernels) print(integrate(function(x) x^2 * kernelFun(x, kernel = k, rescale = TRUE), -Inf, Inf, rel.tol = 1e-8)) # Integrates to 1.

if (write.pdf) pdf(file = "12-convolutions.pdf", width = 7, height = 5)
par(oma = c(0, 0, 3, 0)) # Set margins: bottom = left = right = 0, top = 3 for the main title.
par(mfrow = c(2, 3)) # Plot a 2x3 array of plots.
par(mar = c(2, 1, 4, 1)) # Set margins: bottom = 2, left = 1, top = 4, right = 1 for the sub-plots.
curve(kernelFun(x, kernel = "uniform", rescale = FALSE, convolution = FALSE), -3, 3, ylim = c(0, 1.01), col = my.colours[1], lwd = 1, lty = 2, xlab = "", ylab = "", main = all.kernels[1], font.main = 1, bty = "n", yaxt = "n") # Individual title with normal (roman, not bold) title.
curve(kernelFun(x, kernel = "uniform", rescale = FALSE, convolution = TRUE), -3, 3, col = my.colours[1], lwd = 2, add = TRUE)
for (i in 2:5) {
  curve(kernelFun(x, kernel = all.kernels[i], rescale = FALSE, convolution = FALSE), ylim = c(0, 1.01), -3, 3, add = FALSE, col = my.colours[i], lwd = 1, lty = 2, ylab = "", main = all.kernels[i], font.main = 1, bty = "n", yaxt = "n")
  curve(kernelFun(x, kernel = all.kernels[i], rescale = FALSE, convolution = TRUE), -3, 3, add = TRUE, col = my.colours[i], lwd = 2)
}
mtext("Convolutions of kernels", outer = TRUE, font = 2) # Main title on top.
dev.off()

# Generate the data.
X <- rnorm(n, mean = 0, sd = 1)
U <- rnorm(n, mean = 0, sd = 1)
Y <- X^2 + U
# Generate the values of x at which pdf_X(x) and E(Y|X=x) are to be estimated.
Xgrid <- seq(-5, 5, length.out = 301) # Array of length 301.

# Plot the data.
plot(X, Y, bty = "n")
points(X, X^2, col = "red", pch = 16, cex = 0.75)

# Estimate the basic kernel density with a good value (Silverman's RoT), with one too low and one too high
bws <- c(bw.nrd(X), bw.nrd(X) / 10, bw.nrd(X) * 10)
myfhat <- kernelDensity(X, Xgrid, bws[1], "gaussian")
myfhat.under <- kernelDensity(X, Xgrid, bws[2], "gaussian")
myfhat.over <- kernelDensity(X, Xgrid, bws[3], "gaussian")

if (write.pdf) pdf(file = "13-bandwidth-choice-consequences.pdf", width = 7, height = 5)
plot(Xgrid, myfhat, type = "l", ylim = c(0, max(myfhat.under)), ylab = "Density", xlab = "x", bty = "n", lwd = 2, main = "Effect of bandwidth on the density estimator")
lines(Xgrid, dnorm(Xgrid), lty = 2, lwd = 2)
rug(X)
lines(Xgrid, myfhat.under, type = "l", col = "red", lwd = 2)
lines(Xgrid, myfhat.over, type = "l", col = "blue", lwd = 2)
legend("topright", paste0(c("True density", "Optimal smth (b=", "Under-smth (b=", "Over-smth (b="), c("", sprintf("%1.3f", bws)), c("", rep(")", 3))), bty = "n", col = c("black", "black", "red", "blue"), lty = c(2, 1, 1, 1), lwd = 2)
dev.off()

# Checking equality to one via numeric integration
stepsize <- Xgrid[2] - Xgrid[1]
sum(myfhat) * stepsize
sum(myfhat.under) * stepsize

# Cross-validating the bandwidth for the density---using exponential spacing for bandwidths and constructing the grid aroud our initial RoT guess plus or minus four times
bw.grid <- exp(seq(log(bws[1] / 4), log(bws[1] * 4), length.out = 101))
DCV.values <- lapply(all.kernels, function(k) DCV(X, bw = bw.grid, kernel = k))
DCV.values <- matrix(unlist(DCV.values), ncol = 5)
min.dcv.indices <- apply(DCV.values, 2, function(x) which.min(x))
opt.bw.dcv <- bw.grid[min.dcv.indices]
names(opt.bw.dcv) <- all.kernels
min.cv <- apply(DCV.values, 2, function(x) min(x, na.rm = TRUE))
print(opt.bw.dcv)

if (write.pdf) pdf(file = "14-DCV.pdf", width = 7, height = 5)
par(oma = c(0, 0, 3, 0))
par(mfrow = c(2, 3))
par(mar = c(4, 4, 4, 1))
for (i in 1:5) {
  plot(bw.grid, DCV.values[, i], ylim = range(DCV.values), col = my.colours[i], type = "l", lwd = 2, xlab = "Bandwidth b", ylab = "DCV(b)", main = all.kernels[i], bty = "n", font.main = 1, log = "x")
  points(opt.bw.dcv[i], min.cv[i], cex = 2, col = "black", pch = 16)
}
mtext("Density cross-validation", outer = TRUE, font = 2)
dev.off()

# Or can find the exact minimum using the minimiser, but multiple minima might be an issue
opt.bw.dcv.num <- sapply(all.kernels, function(k) bw.CV(X, kernel = k))
opt.bw.dcv.num.check <- sapply(all.kernels, function(k) bw.CV(X, kernel = k, opt.fun = "optim", method = "BFGS"))
round(cbind(Grid = opt.bw.dcv, Numerical_NLM = opt.bw.dcv.num, Numerical_BFGS = opt.bw.dcv.num.check), 3)

density.optimal <- lapply(all.kernels, function(k) kernelDensity(X, Xgrid, bw = opt.bw.dcv[k], kernel = k))
density.optimal <- matrix(unlist(density.optimal), ncol = 5)

if (write.pdf) pdf(file = "15-density-optimal.pdf", width = 7, height = 5)
par(oma = c(0, 0, 3, 0))
par(mfrow = c(2, 3))
par(mar = c(4, 5, 4, 1))
for (i in 1:5) {
  plot(Xgrid, density.optimal[, i], ylim = range(density.optimal, 0.5), col = my.colours[i], type = "l", lwd = 2, xlab = "x", ylab = expression(hat(f)[X](x)), main = all.kernels[i], bty = "n", font.main = 1)
  lines(Xgrid, dnorm(Xgrid), lty = 2)
  rug(X)
  legend("topright", legend = bquote(b[DCV] == .(sprintf("%1.2f", opt.bw.dcv[i]))), bty = "n")
}
mtext("Density estimators with cross-validated bandwidths", outer = TRUE, font = 2)
dev.off()

# Now, we estimate the best predictor of Y given X by Nadaraya---Watson estimator for
mymuhat <- kernelSmooth(X, Y, Xgrid, opt.bw.dcv["gaussian"], "gaussian")

if (write.pdf) pdf(file = "16-NW-estimator-optimal.pdf", width = 7, height = 5)
par(mar = c(4, 5, 4, 1))
plot(X, Y, xlab = "x", ylab = expression(hat(E)(Y ~ "|" ~ X == x)), bty = "n", main = paste0("Nadaraya-Watson estimator (b = ", sprintf("%1.2f", opt.bw.dcv["gaussian"]), ")"))
lines(Xgrid, mymuhat, lwd = 2)
dev.off()

# Computing least-squares cross-validation function value for bandwidth 0.42
LSCV(X, Y, opt.bw.dcv["gaussian"], "gaussian")
LSCV(X, Y, opt.bw.dcv["gaussian"] * 0.8, "gaussian") # It is smaller!

# Cross-validating the bandwidth for the Nadaraya---Watson estimator using Least Squares
LSCV.values <- lapply(all.kernels, function(k) LSCV(X, Y, bw = bw.grid, kernel = k))
LSCV.values <- matrix(unlist(LSCV.values), ncol = 5)
LSCV.values[!is.finite(LSCV.values)] <- NA
min.cv.indices <- apply(LSCV.values, 2, function(x) which.min(x))
opt.bw.lscv <- bw.grid[min.cv.indices]
names(opt.bw.lscv) <- all.kernels
min.cv <- apply(LSCV.values, 2, function(x) min(x, na.rm = TRUE))

if (write.pdf) pdf("17-LSCV.pdf", width = 7, height = 5)
plot(bw.grid, LSCV.values[, 1], type = "l", ylim = range(LSCV.values, na.rm = TRUE), ylab = "CV", xlab = "Bandwidth", lwd = 2, main = "LS cross-validation function for the Nadaraya-Watson estimator", bty = "n", log = "x")
for (i in 2:5) lines(bw.grid, LSCV.values[, i], col = my.colours[i], lwd = 2)
points(opt.bw.lscv, min.cv, cex = 2, col = my.colours, pch = 16)
legend("topleft", legend = all.kernels, lwd = 2, col = my.colours, bty = "n")
dev.off()

# Getting the exact minimum is troublesome in this case since the solution is on the boundary, but we can find one for the normal kernel
opt.bw.lscv["gaussian"] <- bw.CV(x = X, y = Y, kernel = "gaussian", CV = "LSCV")
# It can always be done manually
optimise(function(b) LSCV(X, Y, bw = b, kernel = "gaussian"), c(0.1, 0.55))$minimum
opt.bw.lscv["gaussian"]

fhat.opt <- sapply(all.kernels, function(k) kernelDensity(X, Xgrid, opt.bw.dcv[k], kernel = k))
muhat.opt <- sapply(all.kernels, function(k) kernelSmooth(X, Y, Xgrid, opt.bw.lscv[k], kernel = k))

end.time <- Sys.time() # Start clock.
seconds.taken <- difftime(end.time, start.time, units = "s")

if (write.pdf) pdf(file = "18-density-and-regression.pdf", width = 7, height = 6)
par(oma = c(0, 0, 3, 0))
par(mfrow = c(2, 2))
par(mar = c(2, 5, 4, 2))

plot(Xgrid, fhat.opt[, "gaussian"], type = "l", lwd = 2, xlab = "x", ylab = expression(hat(f)[X](x)), ylim = range(fhat.opt, na.rm = TRUE), bty = "n")
lines(Xgrid, dnorm(Xgrid), lty = 2)
rug(X)
title(main = "Gaussian kernel", font.main = 2, line = 3)
title(main = paste0("bandwidth = ", sprintf("%1.3f", opt.bw.dcv["gaussian"])), font.main = 1, line = 0)

plot(Xgrid, fhat.opt[, "uniform"], type = "l", lwd = 2, xlab = "x", ylab = expression(hat(f)[X](x)), bty = "n")
lines(Xgrid, dnorm(Xgrid), lty = 2)
rug(X)
title(main = "Uniform kernel", font.main = 2, line = 3)
title(main = paste0("bandwidth = ", sprintf("%1.3f", opt.bw.dcv["uniform"])), font.main = 1, line = 0)

plot(Xgrid, muhat.opt[, "gaussian"], type = "l", xlab = "x", lwd = 2, ylab = expression(hat(mu)(x)), ylim = range(muhat.opt, Y, na.rm = TRUE), bty = "n")
points(X, Y, cex = 0.4, pch = 16, col = "#00000088") # Semi-transparent points
lines(Xgrid, Xgrid^2, lty = 2)
title(main = paste0("bandwidth = ", sprintf("%1.3f", opt.bw.lscv["gaussian"])), font.main = 1, line = 0)

plot(Xgrid, muhat.opt[, "uniform"], type = "l", xlab = "x", lwd = 2, ylab = expression(hat(mu)(x)), ylim = range(muhat.opt, Y, na.rm = TRUE), bty = "n")
points(X, Y, cex = 0.4, pch = 16, col = "#00000088")
lines(Xgrid, Xgrid^2, lty = 2)
title(main = paste0("bandwidth = ", sprintf("%1.3f", opt.bw.lscv["uniform"])), font.main = 1, line = 0)

top.plot.title <- c("Density and regression estimates", paste0("n = ", n), paste0("Time = ", sprintf("%1.3f", as.numeric(seconds.taken)), " s"))
mtext(top.plot.title, outer = TRUE, line = 1:-1)
dev.off()

