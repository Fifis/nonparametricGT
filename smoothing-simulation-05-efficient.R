# R code for simulating efficient estimators
# Code: Andreï V. Kostyrka
# April 1, 2019
# v1.0: 2019-04-01

rm(list = ls()) # Clear the workspace.
source("nonparametricGT/R/smoothing-functions.R") # Load the functions.
n <- 300 # Number of observations.
MC <- 1000 # Number of Monte-Carlo simulations
sked <- function(x) return(1 + abs(x) + sin(x)) # The true conditional variance function
ngrid <- 21
bw.grid <- exp(seq(log(0.1), log(4), length.out = ngrid))

# Substantial productivity gains can be achieved if parallslisation (unavailable on Windows) is used!
# Estimated run time: 300 s for 1000 simulations on 4 cores
library(parallel)
ncores <- if (.Platform$OS.type == "windows") 1 else detectCores()
system.time({
  a <- mclapply(1:MC, function(i) {
    set.seed(i)
    X <- rchisq(n, 3)
    U <- rnorm(n) * sked(X)
    Y <- 1 + X + U
    # plot(X, Y)
    # abline(1, 1)
    mod.ls <- lm(Y ~ X)
    b.ols <- mod.ls$coefficients
    hat.u2 <- mod.ls$residuals^2

    mod.log.aux <- lm(log(hat.u2) ~ X)
    hat.sigma2.log <- exp(predict(mod.log.aux)) # This is not a consistent estimator---only up to scale, under normality and a bunch of other assumptions
    mod.log <- lm(Y ~ X, weights = 1 / hat.sigma2.log)
    # plot(X, hat.u2, ylim = range(hat.u2, hat.sigma2.log))
    # points(X, hat.sigma2.log, col="red", pch=16)

    # We can estimate the unknown conditional variance function non-parametrically
    # But what bandwidth do we pick? We have several options:
    # --- Silverman's simplest rule of thumb (for one- and two-dimensional data)
    # --- Density cross-validation
    # --- Least-squares cross-validation
    # --- Trying many bandwidths and just looking at the behaviour of the estimator

    bw.rot1 <- 1.06 * sd(X) * n^(-1 / 5)
    bw.rot2 <- sd(X) * n^(-1 / 6)
    bw.dcv <- bw.CV(X, opt.fun = "optim", method = "BFGS")
    bw.lscv <- bw.CV(X, y = Y, CV = "LSCV", opt.fun = "optim", method = "BFGS")
    # plot(bw.grid, sapply(bw.grid, function(b) DCV(X, bw=b)))
    # plot(bw.grid, sapply(bw.grid, function(b) LSCV(x = X, y = hat.u2, bw=b)))
    hat.sigma2.all <- lapply(c(bw.rot1, bw.rot2, bw.dcv, bw.lscv, bw.grid), function(b) kernelSmooth(x = X, y = hat.u2, bw = b))
    mod.robinson.all <- lapply(hat.sigma2.all, function(hatvar) lm(Y ~ X, weights = 1 / hatvar))

    # However, it will not hurt if a second-stage estimator is used
    hat.u2.2step <- lapply(lapply(mod.robinson.all, "[[", "residuals"), function(x) x^2)
    bw.lscv.2step <- optimise(function(b) LSCV(x = X, y = hat.u2.2step[[4]], bw = b), interval = c(0.1, 4))$minimum
    bw.all.2step <- c(bw.rot1, bw.rot2, bw.dcv, bw.lscv.2step, bw.grid)
    hat.sigma2.all.2step <- lapply(1:length(bw.all.2step), function(j) kernelSmooth(x = X, y = hat.u2.2step[[j]], bw = bw.all.2step[j]))
    mod.robinson.all.2step <- lapply(hat.sigma2.all.2step, function(hatvar) lm(Y ~ X, weights = 1 / hatvar))

    mod.robinson.beta <- matrix(unlist(lapply(mod.robinson.all, "[[", "coefficients")), ncol = 2, byrow = TRUE)[, 2]
    mod.robinson.2step.beta <- matrix(unlist(lapply(mod.robinson.all.2step, "[[", "coefficients")), ncol = 2, byrow = TRUE)[, 2]
    betas <- c(mod.ls$coefficients[2], mod.log$coefficients[2], mod.robinson.beta, mod.robinson.2step.beta)
    if (i %% 10 == 0) cat(i, "\n")
    return(betas)
  }, mc.cores = ncores)
})

save(a, file = "Robinson.RData", compress = "xz")
load("Robinson.RData")

a <- data.frame(matrix(unlist(a), nrow = MC, byrow = TRUE))
colnames(a) <- c(
  "OLS", "logVar", "Rob.1s.rot1", "Rob.1s.rot2", "Rob.1s.DCV", "Rob.1s.LSCV", paste0("Rob.1s.", round(bw.grid, 2)),
  "Rob.2s.rot1", "Rob.2s.rot2", "Rob.2s.DCV", "Rob.2s.LSCV", paste0("Rob.2s.", round(bw.grid, 2))
)

plot(NULL, NULL, xlim = c(0.5, 1.5), bty = "n", ylim = c(0, 5), main = "Density of various slope estimators", xlab = "", ylab = "")
lines(density(a$OLS), lty = 1)
lines(density(a$logVar), lty = 2)
lines(density(a$Rob.1s.rot1), col = "blue")
lines(density(a$Rob.1s.rot2), col = "blue", lty = 2)
lines(density(a$Rob.1s.DCV), col = "blue", lty = 3)
lines(density(a$Rob.1s.LSCV), col = "blue", lty = 4)
mycols <- rainbow(ngrid, end = 0.3)
for (i in 7:(7+ngrid-1)) {
  lines(density(a[, i]), col = mycols[i - 6])
}

boxplot(a[, 7:(7+ngrid-1)], names = round(bw.grid, 2), frame = FALSE)

Bias1 <- apply(a[, 7:(7+ngrid-1)], 2, function(x) mean(x - 1))
SD1   <- apply(a[, 7:(7+ngrid-1)], 2, sd)
RMSE1 <- apply(a[, 7:(7+ngrid-1)], 2, function(x) sqrt(mean((x - 1)^2)))

Bias2 <- apply(a[, (ngrid+11):(10+ngrid*2)], 2, function(x) mean(x - 1))
SD2   <- apply(a[, (ngrid+11):(10+ngrid*2)], 2, sd)
RMSE2 <- apply(a[, (ngrid+11):(10+ngrid*2)], 2, function(x) sqrt(mean((x - 1)^2)))

cols <- c("#AA0000AA", "#0000AAAA")
plot(bw.grid, RMSE1, type = "l", bty = "n", ylim = range(Bias1, SD1, RMSE1), lwd = 2, main = "Properties of the Robinson estimator", col = cols[1], xlab = "Bandwidth", ylab = "")
points(bw.grid[which.min(RMSE1)], min(RMSE1), pch = 16)
lines(bw.grid, SD1, lty = 2, col = cols[1], lwd = 1.5)
lines(bw.grid, Bias1, col = cols[1], lwd = 1)

lines(bw.grid, RMSE2, col = cols[2], lwd = 2)
lines(bw.grid, SD2, lty = 2, col = cols[2], lwd = 1.5)
lines(bw.grid, Bias2, col = cols[2], lwd = 1)
abline(h = 0)
legend("right", c("Bias", "SD", "RMSE"), lty = c(1, 2, 1), lwd = c(1, 1.5, 2), bty = "n")
legend("left", c("1-step Robinson", "2-step Robinson"), lty = c(1, 2, 1), lwd = c(1, 1.5, 2), bty = "n")

plot(bw.grid, colMeans(a[, 7:(7+ngrid-1)]), ylim=c(0.5, 1.5), type="l", lwd=2, xlab="Bandwidth", ylab="Slope estimate", bty = "n", col = cols[1], main = "Monte-Carlo CI (50%, 90%, 95%, 99%)")
for (q in c(0.25, 0.75)) lines(bw.grid, apply(a[, 7:(7+ngrid-1)], 2, function(x) quantile(x, q)), lty=2, col = cols[1])
for (q in c(0.05, 0.95)) lines(bw.grid, apply(a[, 7:(7+ngrid-1)], 2, function(x) quantile(x, q)), lty=3, col = cols[1])
for (q in c(0.025, 0.975)) lines(bw.grid, apply(a[, 7:(7+ngrid-1)], 2, function(x) quantile(x, q)), lty=4, col = cols[1])
for (q in c(0.005, 0.995)) lines(bw.grid, apply(a[, 7:(7+ngrid-1)], 2, function(x) quantile(x, q)), lty=5, col = cols[1])
lines(bw.grid, colMeans(a[, (ngrid+11):(10+ngrid*2)]), lwd=2, col="blue")
for (q in c(0.25, 0.75)) lines(bw.grid, apply(a[, (ngrid+11):(10+ngrid*2)], 2, function(x) quantile(x, q)), lty=2, col = cols[2])
for (q in c(0.05, 0.95)) lines(bw.grid, apply(a[, (ngrid+11):(10+ngrid*2)], 2, function(x) quantile(x, q)), lty=3, col = cols[2])
for (q in c(0.025, 0.975)) lines(bw.grid, apply(a[, (ngrid+11):(10+ngrid*2)], 2, function(x) quantile(x, q)), lty=4, col = cols[2])
for (q in c(0.005, 0.995)) lines(bw.grid, apply(a[, (ngrid+11):(10+ngrid*2)], 2, function(x) quantile(x, q)), lty=5, col = cols[2])
points(rep(1, 8), as.numeric(apply(a[, c(3:6)], 2, function(x) quantile(x, c(0.25, 0.75)))), col = )
