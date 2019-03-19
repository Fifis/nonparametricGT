# R code for multi-dimensional density estimation and regression.
# Author: Andreï V. Kostyrka
# v1.0: 2019-03-08
# v1.1: 2019-03-18 (made plotting optional)

rm(list = ls()) # Clear workspace.
source("smoothing-functions.R") # Load the functions.
start.time <- Sys.time() # Start clock.
set.seed(20190308) # Set seed for replication.
write.pdf <- TRUE # Do we want plots to be shown on the screen, or to be written in PDF files?

ker <- "gaussian" # In this simulation; we shall use the Gaussian kernel only; you can change this for any one you want provided that you do not run into the zero denominator problem

# Visualise the two-dimensional kernel in 3D.
x1 <- x2 <- seq(-4, 4, length.out = 51)
kernel.grid <- expand.grid(x1 = x1, x2 = x2)
product.kernel <- kernelFun(kernel.grid$x1, ker) * kernelFun(kernel.grid$x2, ker)
if (write.pdf) pdf(file = "20-kernel.pdf", width = 5, height = 5)
par(mar = c(1, 1, 3, 1))
persp(x1, x2, matrix(product.kernel, nrow = length(x1)), theta = 30, phi = 30, ticktype = "detailed", zlab = "K(X)", main = "2D Gaussian product kernel")
dev.off()

n <- 1000
f.fun <- function(x1, x2) dchisq(x1, 5) * dchisq(x2, 8) # This is the true density
mu.fun <- function(x1, x2) 1 + x1 + x2 + 3 * sin(x1) + 3 * cos(0.5 * x2) # This is the true function

X1 <- rchisq(n, 5)
X2 <- rchisq(n, 8)
X <- cbind(X1, X2)
X1max <- quantile(X1, 0.995) # End of grid for visualisation
X2max <- quantile(X2, 0.995)
mu <- mu.fun(X1, X2)
Y <- mu + rnorm(n, sd = 2) # These are the observed values

ngrid <- 41
X1grid <- seq(0, X1max, length.out = ngrid)
X2grid <- seq(0, X2max, length.out = ngrid)
mygrid <- as.matrix(expand.grid(X1grid, X2grid))

f.matrix <- matrix(f.fun(mygrid[, 1], mygrid[, 2]), nrow = ngrid) # True density
mu.matrix <- matrix(mu.fun(mygrid[, 1], mygrid[, 2]), nrow = ngrid) # True conditional expectation

if (!("plot3D" %in% rownames(installed.packages()))) install.packages("plot3D")
library(plot3D) # We can make a wireframe plot using this extra library

if (write.pdf) pdf(file = "21-theoretical-values-3d.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
par(mar = c(1, 2, 3, 1))
p <- persp3D(X1grid, X2grid, f.matrix, xlim = c(0, X1grid[ngrid]), ylim = c(0, X2grid[ngrid]), zlim = range(f.matrix), xlab = "X1", ylab = "X2", zlab = "Density", theta = 30, phi = 20, colvar = NULL, alpha = 0.5, border = "black", facets = NA, ticktype = "detailed", main = "True density")
points3D(X1, X2, rep(0, n), pch = 16, col = "red", cex = 0.5, alpha = 0.75, add = TRUE)

p <- persp3D(X1grid, X2grid, mu.matrix, xlim = c(0, X1grid[ngrid]), ylim = c(0, X2grid[ngrid]), zlim = range(Y), xlab = "X1", ylab = "X2", zlab = "mu", theta = 30, phi = 20, colvar = NULL, alpha = 0.5, border = "black", facets = NA, ticktype = "detailed", main = "True regression function")
points3D(X1, X2, Y, pch = 16, col = "red", cex = 0.5, alpha = 0.75, add = TRUE)
dev.off()

# Sometimes, a 3D plot is hard to read, so people provide 2D contour lines at a higher grid resolution. We make them pretty
nfine <- 101
X1grid.fine <- seq(0, X1max, length.out = nfine)
X2grid.fine <- seq(0, X2max, length.out = nfine)
mygrid.fine <- as.matrix(expand.grid(X1grid.fine, X2grid.fine))

nlev <- 25 # Using 25 contour levels
tcol <- rev(rainbow(nlev, start = 0, end = 0.7, v = 0.9))

if (write.pdf) pdf(file = "22-theoretical-values-2d.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
contour(X1grid.fine, X2grid.fine, matrix(f.fun(mygrid.fine[, 1], mygrid.fine[, 2]), nrow = nfine), xlab = expression(X[1]), ylab = expression(X[2]), bty = "n", levels = round(seq(0, max(f.matrix) * 0.99, length.out = nlev), 5), col = tcol, main = "True density levels")
contour(X1grid.fine, X2grid.fine, matrix(mu.fun(mygrid.fine[, 1], mygrid.fine[, 2]), nrow = nfine), xlab = expression(X[1]), ylab = expression(X[2]), bty = "n", levels = round(seq(min(mu.matrix) * 1.01, max(mu.matrix) * 0.99, length.out = nlev), 2), col = tcol, main = "True conditional expectation levels")
dev.off()

# Now, we can plot a preliminary estimator for some reasonable bandwidth
bnaive <- sd(c(X1, X2)) * n^(-1 / 6) # Getting the bandwidth for visualisation with Silverman's rule
fhat <- kernelDensity(x = cbind(X1, X2), xgrid = mygrid, bw = bnaive)
muhat <- kernelSmooth(x = cbind(X1, X2), y = Y, xgrid = mygrid, bw = bnaive)
if (write.pdf) pdf(file = "23-estimator.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
par(mar = c(1, 2, 3, 1))
persp3D(X1grid, X2grid, matrix(fhat, nrow = ngrid), xlab = "X1", ylab = "X2", zlab = "Density", colvar = NULL, alpha = 0.5, border = "black", facets = NA, main = paste0("Density estimator (b=", round(bnaive, 3), ")"))
persp3D(X1grid, X2grid, matrix(muhat, nrow = ngrid), xlab = "X1", ylab = "X2", zlab = "mu", colvar = NULL, alpha = 0.5, border = "black", facets = NA, main = paste0("Regression estimator (b=", round(bnaive, 3), ")"))
dev.off()


# The tricky part: density=cross-validating two bandwidths or one bandwidth
# Method 1, more computationally difficult: each bandwidth separately
# This section takes 60 seconds on 12 cores under Mac/Linux; on Windows, 60*12 seconds
nbgrid <- 25
b1grid <- b2grid <- exp(seq(log(bnaive / 4), log(bnaive * 4), length.out = nbgrid))
step.size <- log(b1grid[2] / b1grid[1])
bw.grid <- expand.grid(b1 = b1grid, b2 = b2grid)
system.time({ # Takes around 60 seconds on 12 cores
  if (.Platform$OS.type == "windows") { # On Windows, this must be computed single-threadedly
    DCV.values <- apply(bw.grid, 1, function(b) DCV(X, bw = b, kernel = ker))
  } else {
    library(parallel) # On everything else, for this plot, we avaluate DCV on a grid in parallel
    DCV.values <- unlist(mclapply(1:nrow(bw.grid), function(i) DCV(X, bw = as.numeric(bw.grid[i, ]), kernel = ker, rescale = TRUE), mc.cores = detectCores()))
  }
})

DCV.values <- matrix(DCV.values, ncol = nbgrid)
min.dcv.index <- which.min(DCV.values)
opt.bw.dcv <- as.numeric(bw.grid[min.dcv.index, ])
cat(opt.bw.dcv, "- DCV =", min(DCV.values), "\n") # We can improve this values a bit now that we see there are no multiple optima
opt.bw.dcv.fine <- optim(opt.bw.dcv, function(b) DCV(X, bw = b, kernel = ker), control = list(trace = 2))
cat(opt.bw.dcv.fine$par, "- DCV =", opt.bw.dcv.fine$value, "\n") # We were really close
plotcomplex <- function() {
  p <- persp(log(b1grid), log(b2grid), DCV.values, theta = 45, phi = 40, main = "DCV of two bandwidths", axes = FALSE)
  points(trans3D(log(opt.bw.dcv[1]), log(opt.bw.dcv[2]), min(DCV.values), p), pch = 16, col = "red", cex = 2)
  points(trans3D(log(opt.bw.dcv.fine$par[1]), log(opt.bw.dcv.fine$par[2]), opt.bw.dcv.fine$value, p), pch = 16, col = "orange", cex = 2)
  textpos1 <- trans3D(log(b1grid[seq(1, nbgrid, 3)]), rep(min(log(b2grid) - 2 * step.size), nbgrid)[seq(1, nbgrid, 3)], rep(min(DCV.values), nbgrid)[seq(1, nbgrid, 3)], p)
  text(textpos1$x, textpos1$y, labels = as.character(round(b1grid, 2))[seq(1, nbgrid, 3)], cex = 0.7)
  textpos2 <- trans3D(rep(max(log(b2grid)) + 2 * step.size, nbgrid)[seq(1, nbgrid, 3)], log(b2grid)[seq(1, nbgrid, 3)], rep(min(DCV.values), nbgrid)[seq(1, nbgrid, 3)], p)
  text(textpos2$x, textpos2$y, labels = as.character(round(b2grid, 2))[seq(1, nbgrid, 3)], cex = 0.7)
  textpos3 <- trans3D(c(mean(log(b1grid)), max(log(b1grid)) + 4 * step.size), c(min(log(b2grid)) - 4 * step.size, mean(log(b2grid))), rep(min(DCV.values), 2), p)
  text(textpos3$x, textpos3$y, labels = c("b1", "b2"))
  return(p)
}

# Method 2, assuming the bandwidth is the same---searching on the diagonal of the plot
DCV.values.same <- sapply(b1grid, function(b) DCV(X, bw = b, kernel = ker))
all(diag(DCV.values) == DCV.values.same)
opt.bw.dcv.same.fine <- optimise(function(b) DCV(X, bw = b, kernel = ker), interval = range(b1grid))
if (write.pdf) pdf("24-DCV-3D.pdf", 9, 5)
par(mfrow = c(1, 2))
plot(b1grid, DCV.values.same, type = "l", ylab = "DCV", xlab = "Unique bandwidth", bty = "n", log = "x", col = "blue", main = "DCV of one bandwidth")
points(b1grid[which.min(DCV.values.same)], min(DCV.values.same), col = "forestgreen", pch = 16, cex = 2)
points(opt.bw.dcv.same.fine$minimum, opt.bw.dcv.same.fine$objective, col = "green", pch = 16, cex = 2)
par(mar = c(1, 1, 3, 1))
p <- plotcomplex()
lines(trans3D(log(b1grid), log(b1grid), DCV.values.same, p), lwd = 2, col = "blue")
points(trans3D(log(opt.bw.dcv.same.fine$minimum), log(opt.bw.dcv.same.fine$minimum), opt.bw.dcv.same.fine$objective, p), pch = 16, col = "forestgreen", cex = 1.5)
dev.off()

# Cross-validating the bandwidths for Nadaraya---Watson estimator via Least Squares
# Method 1, more computationally difficult: each bandwidth separately
# This section takes 80 seconds on 12 cores under Mac/Linux; on Windows, 80*12 seconds
# Here, we are using the density-cross-validated bandwith as the initial values
nbgrid <- 31
b1grid <- opt.bw.dcv.fine$par[1] * exp(seq(-2, 0.3, length.out = nbgrid))
b2grid <- opt.bw.dcv.fine$par[2] * exp(seq(-2, 0.3, length.out = nbgrid))
step.size <- log(b1grid[2] / b1grid[1])
bw.grid <- expand.grid(b1 = b1grid, b2 = b2grid)

system.time({
  if (.Platform$OS.type == "windows") {
    LSCV.values <- apply(bw.grid, 1, function(b) LSCV(X, Y, bw = b, kernel = ker)) # Might take two minutes
  } else {
    LSCV.values <- unlist(mclapply(1:nrow(bw.grid), function(i) LSCV(X, Y, bw = as.numeric(bw.grid[i, ]), kernel = ker, rescale = TRUE), mc.cores = detectCores()))
  }
})

LSCV.values <- matrix(LSCV.values, ncol = nbgrid)
min.lscv.index <- which.min(LSCV.values)
opt.bw.lscv <- as.numeric(bw.grid[min.lscv.index, ])
cat(opt.bw.lscv, "- LSCV =", min(LSCV.values), "\n")
opt.bw.lscv.fine <- optim(opt.bw.lscv, function(b) LSCV(X, Y, bw = b, kernel = ker), control = list(trace = 2))
cat(opt.bw.lscv.fine$par, "- LSCV =", opt.bw.lscv.fine$value, "\n") # We were really close

plotcomplex <- function() {
  p <- persp(log(b1grid), log(b2grid), LSCV.values, theta = 45, phi = 40, main = "LSCV of two bandwidths", axes = FALSE)
  points(trans3D(log(opt.bw.lscv[1]), log(opt.bw.lscv[2]), min(LSCV.values), p), pch = 16, col = "red", cex = 2)
  points(trans3D(log(opt.bw.lscv.fine$par[1]), log(opt.bw.lscv.fine$par[2]), opt.bw.lscv.fine$value, p), pch = 16, col = "orange", cex = 2)
  textpos1 <- trans3D(log(b1grid[seq(1, nbgrid, 3)]), rep(min(log(b2grid) - 2 * step.size), nbgrid)[seq(1, nbgrid, 3)], rep(min(LSCV.values), nbgrid)[seq(1, nbgrid, 3)], p)
  text(textpos1$x, textpos1$y, labels = as.character(round(b1grid, 2))[seq(1, nbgrid, 3)], cex = 0.7)
  textpos2 <- trans3D(rep(max(log(b1grid)) + 2 * step.size, nbgrid)[seq(1, nbgrid, 3)], log(b2grid)[seq(1, nbgrid, 3)], rep(min(LSCV.values), nbgrid)[seq(1, nbgrid, 3)], p)
  text(textpos2$x, textpos2$y, labels = as.character(round(b2grid, 2))[seq(1, nbgrid, 3)], cex = 0.7)
  textpos3 <- trans3D(c(mean(log(b1grid)), max(log(b1grid)) + 5 * step.size), c(min(log(b2grid)) - 5 * step.size, mean(log(b2grid))), rep(min(LSCV.values), 2), p)
  text(textpos3$x, textpos3$y, labels = c("b1", "b2"))
  return(p)
}

# Method 2, assuming the bandwidth is the same---searching on the diagonal of the plot
LSCV.values.same <- sapply(b2grid, function(b) LSCV(X, Y, bw = b, kernel = ker))
opt.bw.lscv.same.fine <- optimise(function(b) LSCV(X, Y, bw = b, kernel = ker), interval = range(b1grid))
if (write.pdf) pdf("25-LSCV-3D.pdf", 9, 5)
par(mfrow = c(1, 2))
plot(b2grid, LSCV.values.same, type = "l", ylab = "LSCV", bty = "n", xlab = "Unique bandwidth", log = "x", col = "blue", main = "LSCV of one bandwidth")
points(b2grid[which.min(LSCV.values.same)], min(LSCV.values.same), col = "forestgreen", pch = 16, cex = 2)
points(opt.bw.lscv.same.fine$minimum, opt.bw.lscv.same.fine$objective, col = "green", pch = 16, cex = 2)
par(mar = c(1, 1, 3, 1))
p <- plotcomplex()
lines(trans3D(log(b2grid), log(b2grid), LSCV.values.same, p), lwd = 2, col = "blue")
points(trans3D(log(opt.bw.lscv.same.fine$minimum), log(opt.bw.lscv.same.fine$minimum), opt.bw.lscv.same.fine$objective, p), pch = 16, col = "forestgreen", cex = 1.5)
dev.off()

# We take the LSCV bandwidth
fhat.opt <- kernelDensity(X, mygrid, opt.bw.lscv.fine$par, kernel = ker)
muhat.opt <- kernelSmooth(X, Y, mygrid, opt.bw.lscv.fine$par, kernel = ker)

end.time <- Sys.time() # Start clock.
seconds.taken <- difftime(end.time, start.time, units = "s")

if (write.pdf) pdf(file = "26-density-and-regression-multi.pdf", width = 9, height = 5)
par(oma = c(0, 0, 3, 0))
par(mfrow = c(1, 2))
par(mar = c(1, 1, 3, 1))

p <- persp3D(X1grid, X2grid, matrix(fhat.opt, nrow = ngrid), xlim = c(0, X1grid[ngrid]), ylim = c(0, X2grid[ngrid]), zlim = range(f.matrix), xlab = "X1", ylab = "X2", zlab = "Density", theta = 30, phi = 20, colvar = NULL, alpha = 0.5, border = "black", facets = NA, ticktype = "detailed", main = "")
points3D(X1, X2, rep(0, n), pch = 16, col = "red", cex = 0.5, alpha = 0.75, add = TRUE)

p <- persp3D(X1grid, X2grid, matrix(muhat, nrow = ngrid), xlim = c(0, X1grid[ngrid]), ylim = c(0, X2grid[ngrid]), zlim = range(Y), xlab = "X1", ylab = "X2", zlab = expression(mu), theta = 30, phi = 20, colvar = NULL, alpha = 0.5, border = "black", facets = NA, ticktype = "detailed", main = "")
points3D(X1, X2, Y, pch = 16, col = "red", cex = 0.5, alpha = 0.75, add = TRUE)

top.plot.title <- c("Density and regression estimates", paste0("n = ", n, ", time(sec) = ", round(as.numeric(seconds.taken), 1)), paste0("LSCV bandwidth = (", paste(round(opt.bw.lscv.fine$par, 3), collapse = ", "), ")"))
mtext(top.plot.title, outer = TRUE, line = 1:-1)
dev.off()

if (write.pdf) pdf(file = "27-density-and-regression-multi-levels.pdf", width = 9, height = 5)
par(oma = c(0, 0, 3, 0))
par(mfrow = c(1, 2))
par(mar = c(4, 4, 2, 1))
contour(X1grid.fine, X2grid.fine, matrix(kernelDensity(X, mygrid.fine, opt.bw.lscv.fine$par, kernel = ker), nrow = nfine), xlab = expression(X[1]), ylab = expression(X[2]), bty = "n", levels = round(seq(0, max(fhat.opt) * 0.99, length.out = nlev), 5), col = tcol)
contour(X1grid.fine, X2grid.fine, matrix(kernelSmooth(X, Y, mygrid.fine, opt.bw.lscv.fine$par, kernel = ker), nrow = nfine), xlab = expression(X[1]), ylab = expression(X[2]), bty = "n", levels = round(seq(min(muhat.opt) * 1.01, max(muhat.opt) * 0.99, length.out = nlev), 2), col = tcol)
top.plot.title <- c("Density and regression estimates", paste0("n = ", n, ", time(sec) = ", round(as.numeric(seconds.taken), 1)), paste0("LSCV bandwidth = (", paste(round(opt.bw.lscv.fine$par, 3), collapse = ", "), ")"))
mtext(top.plot.title, outer = TRUE, line = 1:-1)
dev.off()
