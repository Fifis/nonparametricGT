# Original code and idea: Gautam Tripathi, 2017-03-08
# Rewrite: Andreï V. Kostyrka
# v0.0.1: 2019-03-07
# v0.0.2: 2019-03-18 (fixed an error in DCV that caused severe over-smoothing)
# v0.1: 2021-09-10 Made it into an R package, uploaded to GitHub

#' Basic univatiate kernel functions
#'
#' Computes 5 most popular kernel functions with the potential of returning an analytical convolution kernel for density cross-validation.
#'
#' @param x A numeric vector of values at which to compute the kernel function.
#' @param kernel Kernel type: uniform, Epanechnikov, triangular, quartic, or Gaussian.
#' @param rescale Logical: rescale to unit variance? If \code{TRUE}, ensures
#' \eqn{\int_{-\infty}^{+\infty} x^2 k(x) = \sigma^2_k = 1}{\int_{-Inf}^{+Inf} x^2 k(x) = \sigma^2_k = 1}.
#' This is useful because in this case, the constant \eqn{k_2}{k_2} in formulæ 3.12 and 3.21
#' from \insertRef{silverman1986density} is equal to 1.
#' @param convolution Logical: return the convolution kernel? (Useful for density cross-validation.)
#'
#' @return A numeric vector of the same length as input.
#' @importFrom Rdpack reprompt
#' @export
#'
#' @examples
#' all.kernels <- c("uniform", "triangular", "epanechnikov", "quartic", "gaussian")
#' my.colours <- c("#000000CC", "#0000CCCC", "#CC0000CC", "#00AA00CC", "#BB8800CC")
#' curve(kernelFun(x, kernel = "uniform", rescale = FALSE), -2, 2, n = 501, ylim = c(0, 1.1),
#' col = my.colours[1], lwd = 2, ylab = "Kernel", main = "Kernels used in smoothing", bty = "n")
#' for (i in 2:5) curve(kernelFun(x, kernel = all.kernels[i], rescale = FALSE), -2, 2,
#'   add = TRUE, col = my.colours[i], lwd = 2)
#' legend("topright", legend = all.kernels, lwd = 2, col = my.colours, bty = "n")
#'
#' # All kernels integrate to one
#' for (k in all.kernels) print(integrate(function(x) kernelFun(x, kernel = k, rescale = FALSE),
#'   lower = -Inf, upper = Inf, rel.tol = 1e-8))
#' for (k in all.kernels) print(integrate(function(x) kernelFun(x, kernel = k, rescale = TRUE),
#'   lower = -Inf, upper = Inf, rel.tol = 1e-8))
#' # Without scaling, the variance is different for different kernels
#' for (k in all.kernels) print(integrate(function(x) x^2*kernelFun(x, kernel=k, rescale=FALSE),
#'   lower = -Inf, upper = Inf, rel.tol = 1e-8))
#' for (k in all.kernels) print(integrate(function(x) x^2*kernelFun(x, kernel=k, rescale=TRUE),
#'   lower = -Inf, upper = Inf, rel.tol = 1e-8))
kernelFun <- function(x,
                      kernel = c("gaussian", "uniform", "triangular", "epanechnikov", "quartic"),
                      rescale = TRUE,
                      convolution = FALSE
) {
  kernel <- kernel[1]
  adj.factor <- switch(kernel,
    uniform = sqrt(3),
    triangular = sqrt(6),
    epanechnikov = sqrt(5),
    quartic = sqrt(7),
    gaussian = 1
  )
  if (!rescale) adj.factor <- 1
  x <- x / adj.factor
  if (!convolution) {
    k <- switch(kernel,
      uniform = 1 / 2 * (abs(x) < 1),
      triangular = (1 - abs(x)) * (abs(x) < 1),
      epanechnikov = 3 / 4 * (1 - x^2) * (abs(x) < 1),
      quartic = 15 / 16 * (1 - x^2)^2 * (abs(x) < 1),
      gaussian = stats::dnorm(x)
    )
  } else {
    k <- switch(kernel,
      uniform = 1 / 4 * (2 - abs(x)) * (abs(x) < 2),
      triangular = 1 / 6 * ((3 * abs(x)^3 - 6 * x^2 + 4) * (abs(x) <= 1) + (8 - 12 * abs(x) + 6 * x^2 - abs(x)^3) * (abs(x) > 1 & abs(x) < 2)),
      epanechnikov = 3 / 160 * (32 - 40 * x^2 + 20 * abs(x)^3 - abs(x)^5) * (abs(x) < 2),
      quartic = 225 / 256 * (-128 / 105 * x^2 + 16 / 15 * x^4 + 256 / 315 + 4 / 105 * abs(x)^7 - 1 / 630 * abs(x)^9 - 8 / 15 * abs(x)^5) * (abs(x) < 2),
      gaussian = stats::dnorm(x, sd = sqrt(2))
    )
  }
  k <- k / adj.factor
  return(k)
}

#' Kernel weight matrix computation
#'
#' @param x A numeric vector or numeric matrix of predictors.
#' @param xgrid A numeric vector or numeric matrix with \code{ncol(xgrid) = ncol(x)}.
#' @param bw Bandwidth: a scalar or a vector of the same length as ncol(x).
#' @param kernel Passed to \code{kernelFun}.
#' @param rescale Passed to \code{kernelFun}.
#' @param convolution Passed to \code{kernelFun}.
#'
#' @seealso \code{\link{kernelFun}} for bare-bones kernels.
#' @return A numeric matrix of size \code{nrow(x) x nrow(xgrid)} with kernel weights.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- sort(rnorm(30))
#' kw <- kernelWeights(x, xgrid = seq(-3, 3, 0.05), bw = bw.rot(x))
#' image(kw)
#' mycols <- rainbow(length(x), end = 0.7, v = 0.8)
#' plot(NULL, NULL, xlim = c(-3, 3), ylim = c(0, max(kw)), bty = "n", xlab = "x", ylab = "Weights")
#' for (i in 1:length(x)) lines(seq(-3, 3, 0.05), kw[i, ], col = mycols[i])
#' abline(v = x, col = mycols, lty = 3)
kernelWeights <- function(x,
                          xgrid = NULL,
                          bw,
                          kernel = "gaussian",
                          rescale = TRUE,
                          convolution = FALSE
) {
  if (is.null(xgrid)) xgrid <- x # If no grid was passed, use existing data points as the grid
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (one.dim) {
    if (length(bw) > 1) stop("For one-dimensional kernel weights, the bandwidth must be a scalar.")
    diffs <- outer(x, xgrid, "-")
    PK <- kernelFun(diffs / bw, kernel = kernel, rescale = rescale, convolution = convolution)
  } else {
    if (!is.matrix(x)) stop("x should be either a numeric vector or a numeric matrix.")
    s <- ncol(x) # The dimension of data
    if (ncol(x) != ncol(xgrid)) stop("x and xgrid must have the same number of columns.")
    if (length(bw) == 1) bw <- rep(bw, s)
    if (length(bw) != s) stop("For multi-dimensional kernel weights, the bandwidth must have the same length as ncol(x).")
    nx <- nrow(x)
    nxgrid <- nrow(xgrid)
    PK <- matrix(1, nrow = nx, ncol = nxgrid) # Initialising the product kernel
    for (i in 1:s) {
      K <- kernelFun(outer(x[, i], xgrid[, i], "-") / bw[i], kernel = kernel, rescale = rescale, convolution = convolution)
      PK <- PK * K
    }
  }
  return(PK)
}

#' Kernel density estimation
#'
#' @param x A numeric vector or numeric matrix of predictors.
#' @param xgrid A numeric vector or numeric matrix with \code{ncol(xgrid) = ncol(x)} of points at which the density is estimated.
#' @param bw Bandwidth: a scalar or a vector of the same length as \code{ncol(x)}.
#' @param kernel Which kernel to use? Passed to \code{kernelWeights}.
#' @param rescale Passed to \code{kernelWeights}.
#'
#' @return \code{kernelDensity}: a numeric vector of kernel density estimator of \code{x} evaluated at \code{xgrid}.
#' \code{kernelSmooth}: a numeric vector of kernel regression estimator values (best predictor of \code{y} given \code{x}) evaluated at \code{xgrid}.
#' @export
#' @seealso \code{\link{kernelWeights}} for kernel weight matrix construction.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(110)
#' y <- x^2 + rnorm(110)
#' xgrid <- seq(-5, 5, length.out = 301)
#' opt.dcv <- bw.CV(x = x)
#' opt.lscv <- bw.CV(x = x, y = y, CV = "LSCV")
#' myfhat       <- kernelDensity(x = x, xgrid = xgrid, bw = opt.dcv, kernel = "gaussian")
#' myfhat.under <- kernelDensity(x = x, xgrid = xgrid, bw = opt.dcv * 0.2, kernel = "gaussian")
#' myfhat.over  <- kernelDensity(x = x, xgrid = xgrid, bw = opt.dcv * 5, kernel = "gaussian")
#' mymuhat       <- kernelSmooth(x = x, y = y, xgrid = xgrid, bw = opt.lscv)
#' mymuhat.under <- kernelSmooth(x = x, y = y, xgrid = xgrid, bw = opt.lscv * 0.2)
#' mymuhat.over  <- kernelSmooth(x = x, y = y, xgrid = xgrid, bw = opt.lscv * 5)
#' plot(xgrid, myfhat, type = "l", ylim = c(0, max(myfhat.under)), ylab = "Density", xlab = "x",
#'   bty = "n", lwd = 2, main = "Effect of the bandwidth on the density estimator")
#' lines(xgrid, dnorm(xgrid), lty = 2, lwd = 2)
#' rug(x)
#' lines(xgrid, myfhat.under, col = "red", lwd = 2)
#' lines(xgrid, myfhat.over, col = "blue", lwd = 2)
#' legend("topright", paste0(c("True density", "Optimal (b=", "Under-smth (b=", "Over-smth (b="),
#'     c("", sprintf("%1.3f", bw.rot(x) * c(0.2, 1, 5))), c("", rep(")", 3))), bty = "n",
#'   col = c("black", "black", "red", "blue"), lty = c(2, 1, 1, 1), lwd = 2)
#' plot(x, y, bty = "n")
#' lines(xgrid, mymuhat, col = "black", lwd = 2)
#' lines(xgrid, xgrid^2, col = "black", lwd = 2, lty = 2)
#' lines(xgrid, mymuhat.under, col = "red", lwd = 2)
#' lines(xgrid, mymuhat.over, col = "blue", lwd = 2)
#' legend("top", paste0(c("True law of y | x", "Optimal (b=", "Under-smth (b=", "Over-smth (b="),
#'     c("", sprintf("%1.3f", bw.rot(x) * c(0.2, 1, 5))), c("", rep(")", 3))), bty = "n",
#'   col = c("black", "black", "red", "blue"), lty = c(2, 1, 1, 1), lwd = 2)
kernelDensity <- function(x,
                          xgrid = NULL,
                          bw = NULL,
                          kernel = "gaussian",
                          rescale = TRUE
) {
  if (is.null(bw)) {
    bw <- bw.rot(x)
    warning(paste0("No bandwidth supplied, using Silverman's one-dimensional rule of thumb: bw = ", round(bw, 5), "."))
  }
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (!one.dim & length(bw) == 1) bw <- rep(bw, ncol(x))
  K <- kernelWeights(x = x, xgrid = xgrid, bw = bw, kernel = kernel, rescale = rescale)
  dens <- colSums(K) / (nrow(K) * prod(bw))
  return(dens)
}

#' @param y A numeric vector of responses (dependent variable).
#' @param LOO Logical: return the leave-one-out estimator?
#' @rdname kernelDensity
#' @export
kernelSmooth <- function(x, y,
                         xgrid = NULL,
                         bw = NULL,
                         kernel = "gaussian",
                         rescale = TRUE,
                         LOO = FALSE
) {
  if (is.null(bw)) {
    bw <- bw.rot(x)
    warning(paste0("No bandwidth supplied, using Silverman's one-dimensional rule of thumb: bw = ", round(bw, 5), "."))
  }
  if (is.numeric(xgrid)) {
    if (LOO & !isTRUE(all.equal(x, xgrid))) stop("The Leave-one-out estimator must use the same xgrid as x or NULL.")
  }
  K <- kernelWeights(x = x, xgrid = xgrid, bw = bw, kernel = kernel, rescale = rescale)
  if (LOO) diag(K) <- 0
  num <- colSums(sweep(K, 1, y, "*"))
  den <- colSums(K)
  muhat <- num / den
  muhat[!is.finite(muhat)] <- NA
  return(muhat)
}

#' Silverman's rule-of-thumb bandwidth
#'
#' A fail-safe function that would return a nice Silverman-like bandwidth suggestion for data for which
#' the standard deviation might be NA or 0.
#'
#' It is obtained under the assumption that the true density is multivariate normal with zero covariances
#' (i.e. a diagonal variance-covariance matrix
#' \eqn{\Sigma = \mathrm{\mathop{diag}}(\sigma^2_k)}{\Sigma = diag(\sigma^2_k)} with
#' \eqn{\det \Sigma = \prod_k \sigma^2_k}{det \Sigma = prod(\sigma^2_k)} and \eqn{\Sigma^{-1} = diag(\sigma^{-2}_k)}{\Sigma = diag(1/\sigma^2_k)}).
#' Then, the formula 4.12 in Silverman (1986) depends only on \eqn{\alpha}{\alpha}, \eqn{\beta}{\beta}
#' (which depend only on the kernel and are fixed for a multivariate normal), and on the L2-norm of the
#' second derivative of the density. The (i, i)th element of the Hessian of multi-variate normal
#' (\eqn{\phi(x_1, \ldots, x_d) = \phi(X)}{\phi(x_1, ..., x_d) = \phi(X)}) is
#' \eqn{\phi(X)(x_i^2 - \sigma^2_i)/\sigma_i^4}{\phi(X)(x_i^2 - \sigma^2_i)/\sigma_i^4}.
#'
#' @param x A numeric vector without non-finite values.
#' @return A bandwidth that might be optimal for non-parametric density estimation of \code{x}.
#' @examples
#' set.seed(1)
#' bw.rot(rnorm(100)) # Should be 0.3787568 in R version 4.0.4
#' @export
bw.rot <- function(x) {
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (one.dim) x <- matrix(x, ncol = 1)
  d <- ncol(x)
  s <- apply(x, 2, stats::sd)
  AK <- (4 / (2*d + 1))^(1 / (d + 4)) # (4.15) from Silverman (1986)
  if (!is.finite(s)) {
    stop("bw.rot: Could not compute the bandwidth, check your data, most likely it has length 1.")
  } else if (s > 0) {
    return((4 / 3)^0.2 * s * length(x)^(-1 / 5))
  } else {
    return(1)
  }
}


#' Density and Least-Squares cross-validation criteria
#'
#' Density cross-validation is performed based purely on the values of \code{x} for the kernel density estimator.
#' Least-squares cross-validation is performed for the Nadaraya-Watson estimator of E(Y|X).
#'
#' @param x A numeric vector or numeric matrix.
#' @param bw A numeric scalar, vector, or matrix of candidate bandwidths.
#' @param kernel Which kernel to use? Passed to \code{kernelWeights}.
#' @param rescale Passed to \code{kernelWeights}.
#' @param same Logical: use the same bandwidth in all dimensions for multi-variate x?
#' @param parallel Logical: use \code{parallel::mclapply} to speed up computation?
#' @param ncores The number of cores to use for parallel computation. Reset to 1 on Windows machines.
#'
#' @return For one-dimensional \code{x}, returns a numeric vector of the same length as \code{length(bw)}.
#' For matrix \code{x}, returns a numeric vector of the same length as \code{nrow(bw).}
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' y <- x^2 + rnorm(100)
#' plot(x, y, bty = "n")
#' bw.grid <- exp(seq(log(bw.rot(x) / 4), log(bw.rot(x) * 2), length.out = 101))
#' DCV.values <- DCV(x, bw = bw.grid)
#' DCV.values <- (DCV.values - min(DCV.values)) / diff(range(DCV.values))
#' LSCV.values <- LSCV(x = x, y = y, bw = bw.grid)
#' LSCV.values <- (LSCV.values - min(LSCV.values)) / diff(range(LSCV.values))
#' plot(bw.grid, DCV.values, bty = "n", xlab = "Bandwidth", type = "l")
#' lines(bw.grid, LSCV.values, col = "red")
DCV <- function(x, bw,
                kernel = "gaussian",
                same = FALSE,
                rescale = TRUE,
                parallel = FALSE,
                ncores = 2) {
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (is.data.frame(bw)) bw <- as.matrix(bw)
  if (one.dim) {
    n <- length(x)
    many.bw <- (length(bw) > 1)
  } else {
    n <- nrow(x)
    many.bw <- (!is.null(dim(bw))) | (is.vector(bw) & length(bw) > 1 & same)
    if (many.bw) {
      if (!is.vector(bw)) bw <- lapply(seq_len(nrow(bw)), function(i) bw[i, ]) # If the input is a matrix, split it into a list
      # Otherwise, [mc]lapply will happily eat a vector
    } else bw <- list(bw) # If there is only one bw, make it a list
  }
  CV <- function(b) { # A sub-function to compute the CV for one BW, parallelisable
    if (any(b <= 0)) return(Inf)
    if (!one.dim & length(b) == 1) b <- rep(b, ncol(x))
    KK <- kernelWeights(x, x, bw = b, kernel = kernel, rescale = rescale, convolution = TRUE)
    term1 <- sum(KK) / (n^2 * prod(b))
    # Computing the LOO estimator efficiently: fhat_i(x) = n/(n-1) * fhat(x) - 1/((n-1)*b^s) * K((X[i] - x)/b)
    fhat <- kernelDensity(x, bw = b, kernel = kernel, rescale = rescale)
    fhat.LOO <- (n * fhat - kernelFun(0, kernel = kernel, rescale = rescale)^length(b) / (prod(b))) / (n - 1)
    term2 <- -2 * mean(fhat.LOO)
    return(term1 + term2)
  }
  if (.Platform$OS.type == "windows" & parallel & ncores > 1) {
    parallel <- FALSE
    warning("On windows, the multi-core functionality via 'mclapply' is not available, using 1 core only. Try Linux or Mac!")
  }
  if (parallel) CV.values <- parallel::mclapply(bw, CV, mc.cores = ncores) else CV.values <- lapply(bw, CV)
  return(unlist(CV.values))
}

#' @param y A numeric vector of responses (dependent variable).
#' @rdname DCV
#' @export
LSCV <- function(x, y, bw,
                 kernel = "gaussian",
                 same = FALSE,
                 rescale = TRUE,
                 parallel = FALSE,
                 ncores = 2) {
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (is.data.frame(bw)) bw <- as.matrix(bw)
  if (one.dim) {
    n <- length(x)
    many.bw <- (length(bw) > 1)
  } else {
    n <- nrow(x)
    many.bw <- (!is.null(dim(bw))) | (is.vector(bw) & length(bw) > 1 & same)
    if (many.bw) {
      if (!is.vector(bw)) bw <- lapply(seq_len(nrow(bw)), function(i) bw[i, ]) # If the input is a matrix, split it into a list
      # Otherwise, [mc]lapply will happily eat a vector
    } else bw <- list(bw) # If there is only one bw, make it a list
  }
  ASE <- function(b) {
    if (any(b <= 0)) return(Inf)
    muhat_i <- kernelSmooth(x = x, y = y, bw = b, kernel = kernel, rescale = rescale, LOO = TRUE)
    m <- mean((y - muhat_i)^2)
    if (!is.finite(m)) m <- Inf
    return(m)
  }
  if (.Platform$OS.type == "windows" & parallel & ncores > 1) parallel <- FALSE
  if (parallel) ASE.values <- parallel::mclapply(bw, ASE, mc.cores = ncores) else ASE.values <- lapply(bw, ASE)
  return(unlist(ASE.values))
}

#' Bandwidth Selectors for Kernel Density Estimation
#'
#' Finds the optimal bandwidth by minimising the density cross-valication or least-squares criteria.
#' Remember that since usually, the CV function is highly non-linear, the return value should be taken with a grain of salt.
#' With non-smooth kernels (such as uniform), it will oftern return the local minimum after starting from a reasonable value.
#' The user might want to standardise the input matrix \code{x} by column (divide by some estimator of scale, like \code{sd}
#' or \code{IQR}) and examine the behaviour of the CV criterion as a function of unique bandwidth (\code{same} argument).
#' If it seems that the optimum is unique, then they may proceed by multiplying the bandwidth by the scale measure,
#' and start the search for the optimal bandwidth in multiple dimensions.
#'
#' @param x A numeric vector or numeric matrix.
#' @param y A numeric vector of responses (dependent variable) if \code{CV == "LSCV"}.
#' @param kernel Which kernel to use? Passed to \code{kernelWeights}.
#' @param start.bw A starting value for the optimiser. If \code{NULL}, then Silverman’s rule of thumb is applied to each column of \code{x} via \code{bw.rot}.
#' @param same Logical: use a unique single bandwidth for all columns of \code{x}? Best used if the variability of data across the columns of \code{x} is of the same order of magnitude.
#' @param CV Density or least-squares cross-validation?
#' @param opt.fun Optimiser to be used. If 'optimise', then \code{...} must contain an interval.
#' @param ret.fun If a custom minimiser is used and returns a list, this function should extract the argmin.
#' @param par.name.in.opt If a custom minimiser is used and returns a list, this is the name of the argument corresponding to the initial parameter value.
#' @param fun.name.in.opt If a custom minimiser is used and returns a list, this is the name of the argument corresponding to the objective function to be minimised.
#' @param ... Passed to \code{opt.fun}.
#'
#' @return An estimate of the optimal bandwidth.
#' @importFrom stats nlm optim optimise nlminb
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' y <- x^2 + rnorm(100)
#' plot(x, y, bty = "n")
#' bw.grid <- seq(0.12, 1, length.out = 101)
#' DCV.vals1 <- DCV(x = x, bw = bw.grid)
#' DCV.vals2 <- DCV(x = x, bw = bw.grid, kernel = "epanechnikov")
#' DCV.vals3 <- LSCV(x = x, y = y, bw = bw.grid, kernel = "epanechnikov")
#' DCV.vals1 <- (DCV.vals1 - min(DCV.vals1)) / diff(range(DCV.vals1))
#' DCV.vals2 <- (DCV.vals2 - min(DCV.vals2)) / diff(range(DCV.vals2))
#' DCV.vals3 <- (DCV.vals3 - min(DCV.vals3)) / diff(range(DCV.vals3))
#' plot(bw.grid, DCV.vals1 - min(DCV.vals1), type = "l", bty = "n",
#'   xlab = "Bandwidth", ylab = "CV criterion")
#' lines(bw.grid, DCV.vals2, lty = 2)
#' lines(bw.grid, DCV.vals3, lwd = 2, col = "red")
#' bw.CV(x)
#' bw.CV(x, kernel = "epanechnikov", start.bw = 0.5, opt.fun = "optim", method = "BFGS")
#' bw.CV(x, y = y, kernel = "epanechnikov", CV = "LSCV")
bw.CV <- function(x, y = NULL, kernel = "gaussian", start.bw = NULL, same = FALSE, CV = c("DCV", "LSCV"),
                  opt.fun = c("nlm", "optim", "nlminb", "optimise"),
                  ret.fun = NULL,
                  par.name.in.opt = NULL,
                  fun.name.in.opt = NULL, ...) {
  opt.fun <- opt.fun[1]
  CV <- CV[1]
  one.dim <- is.vector(x) # Are our data one-dimensional?
  opt <- get(opt.fun)
  f.to.min <- if (CV == "DCV") function(b) DCV(x = x, bw = b, kernel = kernel, same = same) else if (CV == "LSCV") function(b) LSCV(x = x, y = y, bw = b, kernel = kernel, same = same) else stop("bw.CV: 'CV' should be either 'DCV' or 'LSCV'!")
  if (is.null(ret.fun)) ret.fun <- switch(opt.fun, nlm = function(x) x[["estimate"]],
                                          optim = function(x) x[["par"]],
                                          optimise = function(x) x[["minimum"]],
                                          nlminb = function(x) x[["par"]])
  if (one.dim) {
    if (is.null(start.bw)) start.bw <- bw.rot(x)
  } else {
    if (is.null(start.bw)) start.bw <- apply(x, 2, bw.rot)
  }
  opt.result <- switch(opt.fun,
                       nlm = stats::nlm(f = f.to.min, p = start.bw, ...),
                       optim = stats::optim(par = start.bw, fn = f.to.min, ...),
                       optimise = stats::optimise(f = f.to.min, ...),
                       nlminb = stats::nlminb(start = start.bw, objective = f.to.min, ...))
  if (is.null(opt.result)) {
    arg.list <- list()
    if (!is.null(par.name.in.opt)) arg.list[[par.name.in.opt]] <- start.bw
    if (!is.null(fun.name.in.opt)) arg.list[[fun.name.in.opt]] <- f.to.min
    arg.list <- c(arg.list, list(...))
    opt.result <- do.call(opt, args = arg.list)
  }
  return(ret.fun(opt.result))
}

#' Density with conditioning on discrete and continuous variables
#'
#' @param x A numeric vector or a numeric matrix of continuous predictors
#' @param by An integer defining the grouping (all possible unique combinations of discrete predictor values)
#' @param xgrid A numeric vector or numeric matrix with \code{ncol(xgrid) = ncol(x)} of points at which the density is estimated.
#' @param bygrid An integer defining the grouping for the grid (all possible unique combinations of discrete predictor values)
#' @param bw Bandwidth: a scalar or a vector of the same length as \code{ncol(x)}.
#' @param ... Passed to kernelDensity.
#'
#' @return A numeric vector of kernel density estimator of \code{x} evaluated at \code{xgrid}.
#' @export
#' @seealso \code{\link{kernelDensity}} for the underlying function.
#'
#' @examples
kernelMixedDensity <- function(x,
                               by,
                               xgrid = NULL,
                               bygrid = NULL,
                               bw,
                               ...
) {
  by <- as.integer(by)
  one.dim <- is.vector(x)
  if (is.null(xgrid)) xgrid <- x
  if (is.null(bygrid)) bygrid <- by
  ngrid <- if (one.dim) length(xgrid) else nrow(xgrid)
  x.tab <- table(by)
  n <- sum(x.tab)
  res <- numeric(ngrid)
  for (v in sort(unique(by))) {
    s <- by == v
    g <- bygrid == v
    prob <- sum(s) / n
    x.sub <- if (one.dim) x[s] else x[s, ]
    xgrid.sub <- if (one.dim) xgrid[g] else xgrid[g, ]
    res[g] <- kernelDensity(x = x.sub, xgrid = xgrid.sub, bw = bw, ...) * prob
  }
  return(res)
}

#' Kernel smoothing with conditioning on discrete and continuous variables
#'
#' @param x A numeric vector or a numeric matrix of continuous predictors
#' @param y A numeric vector of responses (dependent variable).
#' @param by An integer defining the grouping (all possible unique combinations of discrete predictor values)
#' @param xgrid A numeric vector or numeric matrix with \code{ncol(xgrid) = ncol(x)} of points at which the density is estimated.
#' @param bygrid An integer defining the grouping for the grid (all possible unique combinations of discrete predictor values)
#' @param bw Bandwidth: a scalar or a vector of the same length as \code{ncol(x)}.
#' @param ... Passed to kernelSmooth.
#'
#' @return A numeric vector of kernel density estimator of \code{x} evaluated at \code{xgrid}.
#' @export
#' @seealso \code{\link{kernelSmooth}} for the underlying function.
#'
#' @examples
kernelMixedSmooth <- function(x,
                              y,
                              by,
                              xgrid = NULL,
                              bygrid = NULL,
                              bw,
                              ...
) {
  by <- as.integer(by)
  one.dim <- is.vector(x)
  if (is.null(xgrid)) xgrid <- x
  if (is.null(bygrid)) bygrid <- by
  ngrid <- if (one.dim) length(xgrid) else nrow(xgrid)
  res <- numeric(ngrid)
  for (v in sort(unique(by))) {
    s <- by == v
    g <- bygrid == v
    x.sub <- if (one.dim) x[s] else x[s, ]
    y.sub <- y[s]
    xgrid.sub <- if (one.dim) xgrid[g] else xgrid[g, ]
    res[g] <- kernelSmooth(x = x.sub, y = y.sub, xgrid = xgrid.sub, bw = bw, ...)
  }
  return(res)
}

