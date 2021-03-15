# R code for simulating the conditional distribution of the kernel estimator of pdf_X and E(Y|X).
# Code: Andreï V. Kostyrka
# March 15, 2021
# v1.0: 2019-04-01
# v1.2: 2021-03-15 (adapted to the existing package)

# Make sure you have ImageMagick installed on your computer, otherwise the animation magic will not happen!
# Go to https://www.imagemagick.org/script/download.php and install the latest version.

rm(list = ls()) # Clear the workspace.
source("nonparametricGT/R/smoothing-functions.R") # Load the functions.
write.img <- FALSE # We shall be using PNG and PDF image formats

# We have a design with two discrete Xs and continuous Z
# Now we have to condition on all possible combinations of values of the discrete variable
set.seed(1)
X1den <- function(x) {
  probs <- c(0.10, 0.20, 0.15, 0.10, 0.05, 0.20, 0.05, 0.15) # For X in 0:7
  inds <- match(x, 0:7)
  ret <- probs[inds]
  ret[!is.finite(ret)] <- 0
  return(ret)
}

if (write.img) pdf("40-truefreq", 6, 5)
plot(seq(-5, 15, 1), X1den(seq(-5, 15, 1)), pch = 16, xlab = "X", ylab = "Probability", bty = "n")
for (x in 0:7) lines(rep(x, 2), c(0, X1den(x)))
dev.off()

n <- 1000
educ <- sample(0:7, size = n, replace = TRUE, prob = X1den(0:7))
female <- rbinom(n, 1, 0.5)
age <- rchisq(n, 5)
max.age <- qchisq(0.995, 5)

plot(table(educ), bty = "n", xlab = "Education", ylab = "Frequency")

mu <- function(x1, x2, z) 1 + x1 + 2 * x2 + 2 * x1*x2 + z + 3 * sin(x1) + 3 * cos(0.5 * z)
Y <- mu(educ, female, age) + rnorm(n, sd = 2)
mysteps <- 51
X1vals <- 0:7
X2vals <- 0:1
Zvals <- seq(0, max.age, length.out = mysteps)
Xgrid0 <- expand.grid(X1 = X1vals, X2 = 0, Z = Zvals)
Xgrid1 <- expand.grid(X1 = X1vals, X2 = 1, Z = Zvals)
Ytrue0 <- mu(Xgrid0[, 1], Xgrid0[, 2], Xgrid0[, 3])
Ytrue1 <- mu(Xgrid0[, 1], Xgrid1[, 2], Xgrid0[, 3])

den.true <- dchisq(Xgrid0$Z, 5) * X1den(Xgrid0$X1) # X1, X2, and Z are independent,
# and the density given X2 = 1 is the same as X2 = 0

if (write.img) pdf("41-joint-density", 5, 5)
par(mar = c(0.5, 0.5, 0.5, 0.5))
persp(X1vals, Zvals, matrix(den.true, nrow = length(X1vals)), theta = 45, phi = 25, xlab = "Education", ylab = "Age", zlab = "Joint density")
dev.off()

if (write.img) pdf("42-cond-dist.pdf", 10, 5)
par(oma = c(0, 0, 0, 0))
par(mfrow = c(1, 2))
par(mar = c(1, 1, 3, 1))
persp(X1vals, Zvals, matrix(Ytrue0, nrow = length(X1vals)), zlim = range(Ytrue0, Ytrue1), theta = 45, phi = 25, xlab = "Education", ylab = "Age", zlab = "Y", ticktype = "detailed", main = "True function for males")
persp(X1vals, Zvals, matrix(Ytrue1, nrow = length(X1vals)), zlim = range(Ytrue0, Ytrue1), theta = 45, phi = 25, xlab = "Education", ylab = "Age", zlab = "Y", ticktype = "detailed", main = "True function for females")
dev.off()

X <- cbind(educ, female)
head(X)
X.vals <- apply(X, 1, function(x) paste0(x, collapse = "-"))
X.vals
X.vals.fac <- factor(X.vals)
plot(table(X.vals.fac), bty = "n")
X.vals.int <- as.integer(X.vals.fac)
length(unique(X.vals.int))
X.grid <- rep(Zvals, length(unique(X.vals.int)))
by.grid <- rep(sort(unique(X.vals.int)), each = length(Zvals))
names(by.grid) <- rep(levels(X.vals.fac), each = length(Zvals))

this.bw <- bw.CV(age) * length(unique(X1vals))^0.2

fhat <- kernelMixedDensity(x = age, xgrid = X.grid, by = X.vals.int, bygrid = by.grid, bw = this.bw)
muhat <- kernelMixedSmooth(x = age, y = Y, xgrid = X.grid, by = X.vals.int, bygrid = by.grid, bw = this.bw)

by.grid

if (write.img) pdf("43-mixed-density-estimate.pdf", 10, 5)
par(oma = c(0, 0, 0, 0))
par(mfrow = c(1, 2))
par(mar = c(1, 2, 3, 1))
p0 <- persp(x = c(1, 1.01), y = c(1, 1.1), z = matrix(c(1, 1, 1, 1), 2), xlim = c(-1, 8), ylim = c(0, max.age), zlim = c(0, max(den.true, fhat * 2)), xlab = "Education", ylab = "Age", zlab = "Density estimate", theta = 25, phi = 20, border = NA, ticktype = "detailed", main = "Density estimate for males")
for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, fhat[grep(paste0(i, "-0"), names(by.grid))]  / 0.5, pmat = p0), lwd = 2)
for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, dchisq(Zvals, 5) * X1den(i), pmat = p0), lwd = 1, col = "blue", lty = 2)
p1 <- persp(x = c(1, 1.01), y = c(1, 1.1), z = matrix(c(1, 1, 1, 1), 2), xlim = c(-1, 8), ylim = c(0, max.age), zlim = c(0, max(den.true, fhat * 2)), xlab = "Education", ylab = "Age", zlab = "Density estimate", theta = 25, phi = 20, border = NA, ticktype = "detailed", main = "Density estimate for females")
for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, fhat[grep(paste0(i, "-1"), names(by.grid))]  / 0.5, pmat = p1), lwd = 2)
for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, dchisq(Zvals, 5) * X1den(i), pmat = p1), lwd = 1, col = "red", lty = 2)
dev.off()

if (write.img) pdf("43-mixed-mu-estimate.pdf", 10, 5)
par(oma = c(0, 0, 0, 0))
par(mfrow = c(1, 2))
par(mar = c(1, 2, 3, 1))
p0 <- persp(x = c(1, 1.01), y = c(1, 1.1), z = matrix(c(1, 1, 1, 1), 2), xlim = c(-1, 8), ylim = c(0, max.age), zlim = c(0, max(Ytrue0, Ytrue1, muhat)), xlab = "Education", ylab = "Age", zlab = "Cond. exp. estimate", theta = 25, phi = 20, border = NA, ticktype = "detailed", main = "NW estimate for males")
for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, muhat[grep(paste0(i, "-0"), names(by.grid))], pmat = p0), lwd = 2)
for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, mu(rep(i, length(Zvals)), 0, Zvals), pmat = p0), lwd = 1, col = "blue", lty = 2)
p1 <- persp(x = c(1, 1.01), y = c(1, 1.1), z = matrix(c(1, 1, 1, 1), 2), xlim = c(-1, 8), ylim = c(0, max.age), zlim = c(0, max(Ytrue0, Ytrue1, muhat)), xlab = "Education", ylab = "Age", zlab = "Cond. exp. estimate", theta = 25, phi = 20, border = NA, ticktype = "detailed", main = "NW estimate for females")
for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, muhat[grep(paste0(i, "-1"), names(by.grid))], pmat = p1), lwd = 2)
for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, mu(rep(i, length(Zvals)), 1, Zvals), pmat = p1), lwd = 1, col = "red", lty = 2)
dev.off()

# Animation!
# I am using /tmp because it should be flushed anyway
j <- 1000
for (angle in seq(0, 359.5, 0.25)) {
  j <- j + 1
  png(file = paste0("/tmp/example", j, ".png"), width = 1280, height = 720, type = "cairo", pointsize = 16)
  par(oma = c(0, 0, 0, 0))
  par(mfrow = c(1, 2))
  par(mar = c(1, 2, 3, 1))
  p0 <- persp(x = c(1, 1.01), y = c(1, 1.1), z = matrix(c(1, 1, 1, 1), 2), xlim = c(-1, 8), ylim = c(0, max.age), zlim = c(0, max(den.true, fhat * 2)), xlab = "Education", ylab = "Age", zlab = "Density", theta = angle, phi = 20, border = NA, ticktype = "detailed", main = "Density estimate")
  for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, fhat[grep(paste0(i, "-0"), names(by.grid))]  / 0.5, pmat = p0), lwd = 4, col = "#0000AAAA")
  for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, dchisq(Zvals, 5) * X1den(i), pmat = p0), col = "blue", lwd = 2, lty = 2)
  for (i in 0:7) lines(trans3d(rep(i, length(Zvals)) + 0.1, Zvals, fhat[grep(paste0(i, "-1"), names(by.grid))]  / 0.5, pmat = p0), col = "#AA0000AA", lwd = 4)
  for (i in 0:7) lines(trans3d(rep(i, length(Zvals)) + 0.1, Zvals, dchisq(Zvals, 5) * X1den(i), pmat = p0), col = "red", lwd = 2, lty = 2)
  
  p1 <- persp(x = c(1, 1.01), y = c(1, 1.1), z = matrix(c(1, 1, 1, 1), 2), xlim = c(-1, 8), ylim = c(0, max.age), zlim = c(0, max(Ytrue0, Ytrue1, muhat)), xlab = "Education", ylab = "Age", zlab = "mu", theta = angle, phi = 20, border = NA, ticktype = "detailed", main = "NW estimate")
  for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, muhat[grep(paste0(i, "-0"), names(by.grid))], pmat = p1), lwd = 4, col = "#0000AAAA")
  for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, mu(rep(i, length(Zvals)), 0, Zvals), pmat = p1), col = "blue", lwd = 2, lty = 2)
  for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, muhat[grep(paste0(i, "-1"), names(by.grid))], pmat = p1), col = "#AA0000AA", lwd = 4)
  for (i in 0:7) lines(trans3d(rep(i, length(Zvals)), Zvals, mu(rep(i, length(Zvals)), 1, Zvals), pmat = p1), col = "red", lwd = 2, lty = 2)
  
  mtext("Red = females, blue = males", outer = TRUE, line = -4)
  dev.off()
  if (angle %% 5 == 0) cat(angle, "/", 360, "\n")
}
system("ffmpeg -y -framerate 30 -pattern_type glob -i '/tmp/example*.png' -c:v libx265 -crf 29 -preset slower -tune animation -pix_fmt yuv420p ./48-discrete-3d.mp4")
file.remove(list.files("/tmp/", pattern = ".png", full.names = TRUE))

