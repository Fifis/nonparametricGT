# R code for multi-dimensional density estimation and regression.
# Author: Andreï V. Kostyrka
# v1.0: 2019-03-08
# v1.1: 2019-03-18 (made plotting optional)
# v1.2: 2021-03-10 (implemented a convenient interface for data-driven bandwidth selection)

rm(list = ls()) # Clear the workspace.
source("nonparametricGT/R/smoothing-functions.R") # Load the functions.
start.time <- Sys.time() # Start clock.
set.seed(1) # Set seed for replication.
write.pdf <- TRUE # Do we want plots to be shown on the screen, or to be written in PDF files?

ker <- "gaussian" # In this simulation; we shall use the Gaussian kernel only; you can change this for any one you want provided that you do not run into the zero denominator problem
print("Have you tried all the kernels? Try changing 'ker' to 'gaussian', 'uniform', or 'epanechnikov'!")

# Visualise the two-dimensional kernel in 3D.
x1 <- x2 <- seq(-4, 4, length.out = 51)
kernel.grid <- expand.grid(x1 = x1, x2 = x2)
product.kernel <- kernelFun(kernel.grid$x1, ker) * kernelFun(kernel.grid$x2, ker)
if (write.pdf) pdf(file = "20-kernel.pdf", width = 5, height = 5)
par(mar = c(1, 1, 3, 1))
persp(x1, x2, matrix(product.kernel, nrow = length(x1)), theta = 30, phi = 30, ticktype = "detailed", xlab = "X_1", ylab = "X_2", zlab = "K(X)", main = paste0("2D product kernel (", ker, ")"))
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
p <- persp3D(X1grid, X2grid, f.matrix, xlim = c(0, X1grid[ngrid]), ylim = c(0, X2grid[ngrid]), zlim = range(f.matrix), xlab = "X_1", ylab = "X_2", zlab = "Density", theta = 30, phi = 20, colvar = NULL, alpha = 0.5, border = "black", facets = NA, ticktype = "detailed", main = "True density")
points3D(X1, X2, rep(0, n), pch = 16, col = "red", cex = 0.5, alpha = 0.75, add = TRUE)

p <- persp3D(X1grid, X2grid, mu.matrix, xlim = c(0, X1grid[ngrid]), ylim = c(0, X2grid[ngrid]), zlim = range(Y), xlab = "X_1", ylab = "X_2", zlab = "Y", theta = 30, phi = 20, colvar = NULL, alpha = 0.5, border = "black", facets = NA, ticktype = "detailed", main = "True regression function")
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
contour(X1grid.fine, X2grid.fine, matrix(f.fun(mygrid.fine[, 1], mygrid.fine[, 2]), nrow = nfine), xlab = expression(X[1]), ylab = expression(X[2]), bty = "n", levels = round(seq(0, max(f.matrix) * 0.99, length.out = nlev), 4), col = tcol, main = "True density levels")
contour(X1grid.fine, X2grid.fine, matrix(mu.fun(mygrid.fine[, 1], mygrid.fine[, 2]), nrow = nfine), xlab = expression(X[1]), ylab = expression(X[2]), bty = "n", levels = round(seq(min(mu.matrix) * 1.01, max(mu.matrix) * 0.99, length.out = nlev), 2), col = tcol, main = "True regression function levels")
dev.off()

# Now, we can plot a preliminary estimator for some reasonable bandwidth
bnaive <- apply(cbind(X1, X2), 2, bw.rot) # Getting the bandwidth for visualisation with Silverman's rule in every dimension
fhat <- kernelDensity(x = cbind(X1, X2), xgrid = mygrid, bw = bnaive)
muhat <- kernelSmooth(x = cbind(X1, X2), y = Y, xgrid = mygrid, bw = bnaive)
if (write.pdf) pdf(file = "23-estimator.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
par(mar = c(1, 2, 3, 1))
persp3D(X1grid, X2grid, matrix(fhat, nrow = ngrid), xlab = "X_1", ylab = "X_2", zlab = "Density", colvar = NULL, alpha = 0.5, border = "black", facets = NA, main = paste0("Density estimate (b = [", paste(sprintf("%1.3f", bnaive), collapse = " "), "])"))
persp3D(X1grid, X2grid, matrix(muhat, nrow = ngrid), xlab = "X_1", ylab = "X_2", zlab = "E(Y | X)", colvar = NULL, alpha = 0.5, border = "black", facets = NA, main = paste0("Regression estimate (b = [", paste(sprintf("%1.3f", bnaive), collapse = " "), "])"))
dev.off()


# The tricky part: density=cross-validating two bandwidths or one bandwidth
# Method 1, more computationally difficult: each bandwidth separately
# This section takes 30 seconds on 4 cores under Linux or Mac OS; on Windows, multiply that by 4
nbgrid <- 15
b1grid <- b2grid <- exp(seq(log(mean(bnaive) / 4), log(mean(bnaive) * 4), length.out = nbgrid))
step.size <- log(b1grid[2] / b1grid[1])
bw.grid <- as.matrix(expand.grid(b1 = b1grid, b2 = b2grid))
library(parallel)
cores <- detectCores()
system.time({
  # On Windows, this must be computed single-threadedly
  # On everything else, for this plot, we evaluate the DCV on a grid in parallel
  DCV.values <- DCV(X, bw = bw.grid, kernel = ker, parallel = TRUE, ncores = cores)
})

DCV.values <- matrix(DCV.values, ncol = nbgrid)
min.dcv.index <- which.min(DCV.values)
opt.bw.dcv <- as.numeric(bw.grid[min.dcv.index, ])
cat(opt.bw.dcv, "--- DCV =", min(DCV.values), "\n") # We can improve this values a bit now that we see there are no multiple optima
opt.bw.dcv.fine <- bw.CV(X, start.bw = opt.bw.dcv, kernel = ker, opt.fun = "nlm", print.level = 2)
cat(opt.bw.dcv.fine, "--- DCV =", DCV(X, bw = opt.bw.dcv.fine, kernel = ker), "\n") # We were really close
plotcomplex <- function() {
  p <- persp(log(b1grid), log(b2grid), DCV.values, theta = 45, phi = 40, main = "DCV of two bandwidths", axes = FALSE)
  points(trans3D(log(opt.bw.dcv[1]), log(opt.bw.dcv[2]), min(DCV.values), p), pch = 16, col = "red", cex = 2)
  points(trans3D(log(opt.bw.dcv.fine[1]), log(opt.bw.dcv.fine[2]), DCV(X, opt.bw.dcv.fine, kernel = ker), p), pch = 16, col = "orange", cex = 2)
  textpos1 <- trans3D(log(b1grid[seq(1, nbgrid, 3)]), rep(min(log(b2grid) - 2 * step.size), nbgrid)[seq(1, nbgrid, 3)], rep(min(DCV.values), nbgrid)[seq(1, nbgrid, 3)], p)
  text(textpos1$x, textpos1$y, labels = sprintf("%1.2f", b1grid)[seq(1, nbgrid, 3)], cex = 0.7)
  textpos2 <- trans3D(rep(max(log(b2grid)) + 2 * step.size, nbgrid)[seq(1, nbgrid, 3)], log(b2grid)[seq(1, nbgrid, 3)], rep(min(DCV.values), nbgrid)[seq(1, nbgrid, 3)], p)
  text(textpos2$x, textpos2$y, labels = sprintf("%1.2f", b2grid)[seq(1, nbgrid, 3)], cex = 0.7)
  textpos3 <- trans3D(c(mean(log(b1grid)), max(log(b1grid)) + 4 * step.size), c(min(log(b2grid)) - 4 * step.size, mean(log(b2grid))), rep(min(DCV.values), 2), p)
  text(textpos3$x, textpos3$y, labels = c(expression(b[1]), expression(b[2])))
  return(p)
}

# Method 2, assuming the bandwidth is the same---searching on the diagonal of the plot
DCV.values.same <- DCV(X, bw = b1grid, kernel = ker, same = TRUE, parallel = TRUE, ncores = cores)
all(diag(DCV.values) == DCV.values.same)
opt.bw.dcv.same.fine <- bw.CV(X, kernel = ker, start.bw = b1grid[which.min(DCV.values.same)], same = TRUE, opt.fun = "nlm", print.level = 2)
if (write.pdf) pdf("24-DCV-3D.pdf", 9, 5)
par(mfrow = c(1, 2))
plot(b1grid, DCV.values.same, type = "l", ylab = "DCV", xlab = "Unique bandwidth", bty = "n", log = "x", col = "blue", main = "DCV of one bandwidth", lwd = 2)
points(b1grid[which.min(DCV.values.same)], min(DCV.values.same), col = "forestgreen", pch = 16, cex = 2)
points(opt.bw.dcv.same.fine, DCV(X, bw = opt.bw.dcv.same.fine, kernel = ker, same = TRUE), col = "green", pch = 16, cex = 2)
par(mar = c(1, 1, 3, 1))
p <- plotcomplex()
lines(trans3D(log(b1grid), log(b1grid), DCV.values.same, p), lwd = 2, col = "blue")
points(trans3D(log(opt.bw.dcv.same.fine), log(opt.bw.dcv.same.fine), DCV(X, bw = opt.bw.dcv.same.fine, kernel = ker, same = TRUE), p), pch = 16, col = "green", cex = 1.5)
dev.off()

# Cross-validating the bandwidths for Nadaraya---Watson estimator via Least Squares
# Method 1, more computationally difficult: each bandwidth separately
# This section takes 60 seconds on 4 cores under Linux or Mac OS; on Windows, 4 times longer!
# Here, we are using the density-cross-validated bandwidth as the initial values with some span
# For kernels with bounded support, we need to try much larger values!
nbgrid2 <- 19
mult.seq <- switch(ker, gaussian = exp(seq(-2, 0.3, length.out = nbgrid2)),
                   uniform = exp(seq(0.4, 1.7, length.out = nbgrid2)),
                   triangular = exp(seq(0, 1.6, length.out = nbgrid2)),
                   epanechnikov = exp(seq(0, 1.5, length.out = nbgrid2)),
                   quartic = exp(seq(0, 1.5, length.out = nbgrid2)))
b1grid2 <- opt.bw.dcv.fine[1] * mult.seq
b2grid2 <- opt.bw.dcv.fine[2] * mult.seq
step.size2 <- log(b1grid2[2] / b1grid2[1])
bw.grid2 <- expand.grid(b1 = b1grid2, b2 = b2grid2)

system.time({
  LSCV.values <- LSCV(X, Y, bw = bw.grid2, kernel = ker, parallel = TRUE, ncores = cores)
  LSCV.values[!is.finite(LSCV.values)] <- NA
})

LSCV.values <- matrix(LSCV.values, ncol = nbgrid2)
min.lscv.index <- which.min(LSCV.values)
opt.bw.lscv <- as.numeric(bw.grid2[min.lscv.index, ])
cat(opt.bw.lscv, "--- LSCV =", min(LSCV.values, na.rm = TRUE), "\n")
opt.bw.lscv.fine <- bw.CV(X, y = Y, start.bw = opt.bw.lscv, kernel = ker, CV = "LSCV", opt.fun = "nlm", print.level = 2)
cat(opt.bw.lscv.fine, "--- LSCV =", LSCV(X, Y, opt.bw.lscv.fine, kernel = ker), "\n") # We were really close

plotcomplex2 <- function() {
  p <- persp(log(b1grid2), log(b2grid2), LSCV.values, theta = 45, phi = 40, main = "LSCV of two bandwidths", axes = FALSE)
  points(trans3D(log(opt.bw.lscv[1]), log(opt.bw.lscv[2]), min(LSCV.values, na.rm = TRUE), p), pch = 16, col = "red", cex = 2)
  points(trans3D(log(opt.bw.lscv.fine[1]), log(opt.bw.lscv.fine[2]), LSCV(X, Y, opt.bw.lscv.fine, kernel = ker), p), pch = 16, col = "orange", cex = 2)
  textpos1 <- trans3D(log(b1grid2[seq(1, nbgrid2, 3)]), rep(min(log(b2grid2) - 2 * step.size2), nbgrid2)[seq(1, nbgrid2, 3)], rep(min(LSCV.values, na.rm = TRUE), nbgrid2)[seq(1, nbgrid2, 3)], p)
  text(textpos1$x, textpos1$y, labels = sprintf("%1.2f", b1grid2)[seq(1, nbgrid2, 3)], cex = 0.7)
  textpos2 <- trans3D(rep(max(log(b1grid2)) + 2 * step.size2, nbgrid2)[seq(1, nbgrid2, 3)], log(b2grid2)[seq(1, nbgrid2, 3)], rep(min(LSCV.values, na.rm = TRUE), nbgrid2)[seq(1, nbgrid2, 3)], p)
  text(textpos2$x, textpos2$y, labels = sprintf("%1.2f", b2grid2)[seq(1, nbgrid2, 3)], cex = 0.7)
  textpos3 <- trans3D(c(mean(log(b1grid2)), max(log(b1grid2)) + 5 * step.size2), c(min(log(b2grid2)) - 5 * step.size2, mean(log(b2grid2))), rep(min(LSCV.values, na.rm = TRUE), 2), p)
  text(textpos3$x, textpos3$y, labels = c(expression(b[1]), expression(b[2])))
  return(p)
}

# Method 2, assuming the bandwidth is the same---searching on the diagonal of the plot
LSCV.values.same <- LSCV(X, Y, bw = b2grid2, kernel = ker, same = TRUE, parallel = TRUE, ncores = cores)
LSCV.values.same[!is.finite(LSCV.values.same)] <- NA
opt.bw.lscv.same.fine <- bw.CV(X, Y, start.bw = b2grid2[which.min(LSCV.values.same)] * 1.01, kernel = ker, same = TRUE, CV = "LSCV", opt.fun = "nlm", print.level = 2)
if (write.pdf) pdf("25-LSCV-3D.pdf", 9, 5)
par(mfrow = c(1, 2))
plot(b2grid2, LSCV.values.same, type = "l", ylab = "LSCV", bty = "n", xlab = "Unique bandwidth", log = "x", col = "blue", main = "LSCV of one bandwidth")
points(b2grid2[which.min(LSCV.values.same)], min(LSCV.values.same, na.rm = TRUE), col = "forestgreen", pch = 16, cex = 2)
points(opt.bw.lscv.same.fine, LSCV(X, Y, opt.bw.lscv.same.fine, kernel = ker), col = "green", pch = 16, cex = 2)
par(mar = c(1, 1, 3, 1))
p <- plotcomplex2()
lines(trans3D(log(b2grid2), log(b2grid2), LSCV.values.same, p), lwd = 2, col = "blue")
points(trans3D(log(opt.bw.lscv.same.fine), log(opt.bw.lscv.same.fine), LSCV(X, Y, opt.bw.lscv.same.fine, kernel = ker), p), pch = 16, col = "green", cex = 1.5)
dev.off()

# We take the optimal bandwidths
fhat.opt <- kernelDensity(X, mygrid, opt.bw.dcv.fine, kernel = ker)
muhat.opt <- kernelSmooth(X, Y, mygrid, opt.bw.lscv.fine, kernel = ker)

end.time <- Sys.time() # Start clock.
seconds.taken <- difftime(end.time, start.time, units = "s")

if (write.pdf) pdf(file = "26-density-and-regression-multi.pdf", width = 9, height = 5)
par(oma = c(0, 0, 3, 0))
par(mfrow = c(1, 2))
par(mar = c(1, 1, 3, 1))

p <- persp3D(X1grid, X2grid, matrix(fhat.opt, nrow = ngrid), xlim = c(0, X1grid[ngrid]), ylim = c(0, X2grid[ngrid]), zlim = range(f.matrix), xlab = "X1", ylab = "X2", zlab = "Density", theta = 30, phi = 20, colvar = NULL, alpha = 0.5, border = "black", facets = NA, ticktype = "detailed", main = "")
points3D(X1, X2, rep(0, n), pch = 16, col = "red", cex = 0.5, alpha = 0.75, add = TRUE)

p <- persp3D(X1grid, X2grid, matrix(muhat.opt, nrow = ngrid), xlim = c(0, X1grid[ngrid]), ylim = c(0, X2grid[ngrid]), zlim = range(Y), xlab = "X1", ylab = "X2", zlab = expression(mu), theta = 30, phi = 20, colvar = NULL, alpha = 0.5, border = "black", facets = NA, ticktype = "detailed", main = "")
points3D(X1, X2, Y, pch = 16, col = "red", cex = 0.5, alpha = 0.75, add = TRUE)

top.plot.title <- c("Density and regression estimates",
                    paste0("n = ", n, ", time(sec) = ", round(as.numeric(seconds.taken), 1)),
                    paste0("DCV bandwidth = (", paste(sprintf("%1.3f", opt.bw.dcv.fine), collapse = ", "), ")"),
                    paste0("LSCV bandwidth = (", paste(sprintf("%1.3f", opt.bw.lscv.fine), collapse = ", "), ")"))
mtext(top.plot.title, outer = TRUE, line = 1:-2)
dev.off()

if (write.pdf) pdf(file = "27-density-and-regression-multi-levels.pdf", width = 9, height = 5)
par(oma = c(0, 0, 3, 0))
par(mfrow = c(1, 2))
par(mar = c(4, 4, 2, 1))
contour(X1grid.fine, X2grid.fine, matrix(kernelDensity(X, mygrid.fine, opt.bw.lscv.fine, kernel = ker), nrow = nfine), xlab = expression(X[1]), ylab = expression(X[2]), bty = "n", levels = round(seq(0, max(fhat.opt) * 0.99, length.out = nlev), 4), col = tcol)
contour(X1grid.fine, X2grid.fine, matrix(kernelSmooth(X, Y, mygrid.fine, opt.bw.lscv.fine, kernel = ker), nrow = nfine), xlab = expression(X[1]), ylab = expression(X[2]), bty = "n", levels = round(seq(min(muhat.opt, na.rm = TRUE) * 1.01, max(muhat.opt, na.rm = TRUE) * 0.99, length.out = nlev), 2), col = tcol)
mtext(top.plot.title, outer = TRUE, line = 1:-2)
dev.off()

# Why should we make our functions extensible? Because sometimes the minimised function is very non-linear,
# so in this afterword, I show how to choose the bandwidth via custom optimisers or even stochastic search
# For deterministic optimisation, I shall try Hooke-Jeeves algorithm and Mesh Adaptive Direct Search from `dfoptim`
# Then, I shall use one stochastic method: differential evolution
library(dfoptim)
dfoptim::hjk(par = c(2, 2), fn = function(b) LSCV(x = X, y = Y, bw = b, kernel = "uniform"), control = list(info = TRUE))
# It uses 'par' as the starting values and 'fn' as the function to minimise, and allows a custom 'control' list
# For more info, see ?hjk
bw.hjk <- bw.CV(x = X, y = Y, kernel = "uniform", CV = "LSCV", opt.fun =  "hjk", ret.fun = function(x) x[["par"]],
                start.bw = c(2, 2), par.name.in.opt = "par", fun.name.in.opt = "fn",
                control = list(info = TRUE))
bw.mesh <- bw.CV(x = X, y = Y, kernel = "uniform", CV = "LSCV", opt.fun =  "mads", ret.fun = function(x) x[["par"]],
                 start.bw = c(2, 2), par.name.in.opt = "par", fun.name.in.opt = "fn",
                 control = list(trace = TRUE))

# This one does not require any initial values
library(DEoptim)
library(foreach)
library(doParallel)
registerDoParallel(cores)
set.seed(1)
opt.evolution <- DEoptim(fn = function(b) LSCV(x = X, y = Y, bw = b, kernel = "uniform"), lower = c(0.2, 0.2), upper = c(2, 2), control = list(NP = 30, itermax = 20, trace = 1, parallelType = 2))
set.seed(1)
bw.evolution <- bw.CV(x = X, y = Y, kernel = "uniform", CV = "LSCV", opt.fun =  "DEoptim", ret.fun = function(x) x[["optim"]][["bestmem"]],
                  fun.name.in.opt = "fn",
                  lower = c(0.2, 0.2), upper = c(2, 2), control = list(NP = 30, itermax = 20, trace = 1, parallelType = 2))

bw.all <- rbind(Grid = opt.bw.lscv, Grid_Opt = opt.bw.lscv.fine, HJK = bw.hjk, Mesh = bw.mesh, DEoptim = bw.evolution)
colnames(bw.all) <- c("b1", "b2")
round(cbind(bw.all, CV = LSCV(x = X, y = Y, bw = bw.all, kernel = "uniform")), 3)

# However, some functions are buggy and cannot be used with this interface
library(hydroPSO)
set.seed(1)
opt.swarm <- hydroPSO(fn = function(b) LSCV(x = X, y = Y, bw = b, kernel = "uniform"),
                      lower = c(0.2, 0.2), upper = c(2, 2), control = list(npart = 30, maxit = 20, plot = TRUE, verbose = TRUE, REPORT = 1, parallel = "parallel"))
# ERROR!
bw.swarm <- bw.CV(x = X, y = Y, kernel = "uniform", CV = "LSCV", opt.fun =  "hydroPSO", ret.fun = function(x) x[["par"]],
                  fun.name.in.opt = "fn",
                  lower = c(0.2, 0.2), upper = c(2, 2), control = list(npart = 30, maxit = 20, plot = TRUE, verbose = TRUE, REPORT = 1, parallel = "parallel"))
