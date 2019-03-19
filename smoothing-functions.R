# R code for kernel estimator of pdf_X(x) and E(Y|X=x).
# Original code and idea: Gautam Tripathi, 2017-03-08
# Rewrite: Andreï V. Kostyrka
# v1.0: 2019-03-07
# v1.1: 2019-03-18 (fixed an error in DCV that caused severe over-smoothing)

kernelFun <- function(x, # The values to compute the kernel
                      kernel = "gaussian", # Kernel type: uniform, epanechnikov, triangular, quartic, gaussian
                      rescale = TRUE, # Rescale to unit variance: int_-Inf^+Inf x^2 k(x) = sigma^2_K = 1
                      convolution = FALSE # Return the convolution kernel? Used for CV
) {
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
      gaussian = dnorm(x)
    )
  } else {
    k <- switch(kernel,
      uniform = 1 / 4 * (2 - abs(x)) * (abs(x) < 2),
      triangular = 1 / 6 * ((3 * abs(x)^3 - 6 * x^2 + 4) * (abs(x) <= 1) + (8 - 12 * abs(x) + 6 * x^2 - abs(x)^3) * (abs(x) > 1 & abs(x) < 2)),
      epanechnikov = 3 / 160 * (32 - 40 * x^2 + 20 * abs(x)^3 - abs(x)^5) * (abs(x) < 2),
      quartic = 225 / 256 * (-128 / 105 * x^2 + 16 / 15 * x^4 + 256 / 315 + 4 / 105 * abs(x)^7 - 1 / 630 * abs(x)^9 - 8 / 15 * abs(x)^5) * (abs(x) < 2),
      gaussian = dnorm(x, sd = sqrt(2))
    )
  }
  k <- k / adj.factor
  return(k)
}

kernelWeights <- function(x, # A numeric vector or numeric matrix
                          xgrid = NULL, # A numeric vector or numeric matrix with ncol(xgrid) = ncol(x)
                          bw, # Bandwidth: a scalar or a vector of the same length as ncol(x)
                          kernel = "gaussian", # Passed to kernelFun
                          rescale = TRUE, # Passed to kernelFun
                          convolution = FALSE # Passed to kernelFun
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

# Function for estimating pdf_X(x).
kernelDensity <- function(x, xgrid = NULL, bw, kernel = "gaussian", rescale = TRUE # Arguments passed to kernelWeights(...)
) {
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (!one.dim & length(bw) == 1) bw <- rep(bw, ncol(x))
  K <- kernelWeights(x = x, xgrid = xgrid, bw = bw, kernel = kernel, rescale = rescale)
  dens <- colSums(K) / (nrow(K) * prod(bw))
  return(dens)
}

# Function for estimating E(Y|X=x) (possibly leaving one observation out).
kernelSmooth <- function(x, # Passed to kernelWeights(...)
                         y, # A vector of observations
                         xgrid = NULL, bw, kernel = "gaussian", rescale = TRUE, # Passed to kernelWeights(...)
                         LOO = FALSE # Return the leave-one-out estimator?
) {
  if (is.numeric(xgrid)) {
    if (LOO & !isTRUE(all.equal(x, xgrid))) stop("The Leave-one-out estimator must use the same xgrid as x or NULL.")
  }
  K <- kernelWeights(x = x, xgrid = xgrid, bw = bw, kernel = kernel, rescale = rescale)
  if (LOO) diag(K) <- 0
  num <- colSums(sweep(K, 1, y, "*"))
  den <- colSums(K)
  muhat <- num / den
  return(muhat)
}

# Density cross-validation for the estimator of f_X(x).
DCV <- function(x, bw, kernel = "gaussian", rescale = TRUE) {
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (one.dim) {
    n <- length(x)
  } else {
    n <- nrow(x)
    if (length(bw) == 1) bw <- rep(bw, ncol(x))
  }
  KK <- kernelWeights(x, x, bw = bw, kernel = kernel, rescale = rescale, convolution = TRUE)
  term1 <- sum(KK) / (n^2 * prod(bw))
  # Computing the LOO estimator efficiently: fhat_i(x) = n/(n-1) * fhat(x) - 1/((n-1)*b^s) * K((X[i] - x)/b)
  fhat <- kernelDensity(x, bw = bw, kernel = kernel, rescale = rescale)
  fhat.LOO <- (n * fhat - kernelFun(0, kernel = kernel, rescale = rescale) / (prod(bw))) / (n - 1)
  term2 <- -2 * mean(fhat.LOO)
  return(term1 + term2)
}

# Least-squares cross-validation function for Nadaraya-Watson estimator of E(Y|X).
LSCV <- function(x, y, bw, kernel = "gaussian", rescale = TRUE) {
  muhat_i <- kernelSmooth(x = x, y = y, bw = bw, kernel = kernel, rescale = rescale, LOO = TRUE)
  ASE <- mean((y - muhat_i)^2)
  return(ASE)
}
