# Original code and idea: Gautam Tripathi, 2017-03-08
# Rewrite: Andreï V. Kostyrka
# v0.0.1: 2019-03-07
# v0.0.2: 2019-03-18 (fixed an error in DCV that caused severe over-smoothing)
# v0.1: 2021-03-10 Made it into an R package, uploaded to GitHub
# v0.2: 2021-03-16 Added higher-order kernels, multi-dimensional rule-of-thumb bandwidth, semi-parametric regression

#' Basic univatiate kernel functions
#'
#' Computes 5 most
#'  popular kernel functions of orders 2, 4, and 6 with the potential of returning an analytical convolution kernel for density cross-validation.
#'
#' @param x A numeric vector of values at which to compute the kernel function.
#' @param kernel Kernel type: uniform, Epanechnikov, triangular, quartic, or Gaussian.
#' @param order Kernel order. 2nd-order kernels are always non-negative. kth-order kernels have all moments from 1 to (k-1) equal to zero, which is achieved by having some negative values.
#' @param rescale Logical: rescale to unit variance? If \code{TRUE}, ensures that, for the chosen kernel
#' name, the second-order kernel integrates to 1:
#' \eqn{\int_{-\infty}^{+\infty} x^2 k(x) = \sigma^2_k = 1}.
#' This is useful because in this case, the constant \code{k_2} in formulæ 3.12 and 3.21
#' from \insertCite{silverman1986density;textual}{nonparametricGT} is equal to 1.
#' @param convolution Logical: return the convolution kernel? (Useful for density cross-validation.)
#'
#' @return A numeric vector of the same length as input.
#' @importFrom Rdpack reprompt
#' @export
#' @references
#' \insertRef{silverman1986density}{nonparametricGT}
#'
#'
#'
#' @examples
#' ks <- c("uniform", "triangular", "epanechnikov", "quartic", "gaussian"); names(ks) <- ks
#' os <- c(2, 4, 6); names(os) <- paste0("o", os)
#' cols <- c("#000000CC", "#0000CCCC", "#CC0000CC", "#00AA00CC", "#BB8800CC")
#' put.legend <- function() legend("topright", legend = ks, lty = 1, col = cols, bty = "n")
#' xgrid <- seq(-4, 4, length.out = 301)
#' plot(NULL, NULL, xlim = range(xgrid), ylim = c(0, 1.1),
#'   xlab = "", ylab = "", main = "Unscaled kernels", bty = "n"); put.legend()
#' for (i in 1:5) lines(xgrid, kernelFun(xgrid, kernel = ks[i], rescale = FALSE), col = cols[i])
#' par(mfrow = c(1, 2))
#' plot(NULL, NULL, xlim = range(xgrid), ylim = c(-0.1, 0.8), xlab = "", ylab = "",
#'   main = "4th-order scaled kernels", bty = "n"); put.legend()
#' for (i in 1:5) lines(xgrid, kernelFun(xgrid, kernel = ks[i], order = 4), col = cols[i])
#' plot(NULL, NULL, xlim = range(xgrid), ylim = c(-0.25, 1.2), xlab = "", ylab = "",
#'   main = "6th-order scaled kernels", bty = "n"); put.legend()
#' for (i in 1:5) lines(xgrid, kernelFun(xgrid, kernel = ks[i], order = 6), col = cols[i])
#' par(mfrow = c(1, 1))
#' plot(NULL, NULL, xlim = range(xgrid), ylim = c(-0.25, 1.4), xlab = "", ylab = "",
#'   main = "Convolution kernels", bty = "n"); put.legend()
#' for (i in 1:5) {
#'   for (j in 1:3) lines(xgrid, kernelFun(xgrid, kernel = ks[i], order = os[j],
#'   convolution = TRUE), col = cols[i], lty = j)
#' }; legend("topleft", c("2nd order", "4th order", "6th order"), lty = 1:3, bty = "n")
#'
#' # All kernels integrate to correct values; we compute the moments
#' mom <- Vectorize(function(k, o, m, c) integrate(function(x) x^m * kernelFun(x, k, o,
#'   rescale = FALSE, convolution = c), lower = -Inf, upper = Inf)$value)
#' for (m in 0:6) {
#'   cat("\nComputing integrals of x^", m, " * f(x). \nSimple unscaled kernel:\n", sep = "")
#'   print(round(outer(os, ks, function(o, k) mom(k, o, m = m, c = FALSE)), 4))
#'   cat("Convolution kernel:\n")
#'   print(round(outer(os, ks, function(o, k) mom(k, o, m = m, c = TRUE)), 4))
#' }
#'
kernelFun <- function(x,
                      kernel = c("gaussian", "uniform", "triangular", "epanechnikov", "quartic"),
                      order = c(2, 4, 6),
                      rescale = TRUE,
                      convolution = FALSE
) {
  order <- order[1]
  if (!(order %in% c(2, 4, 6))) stop("Only kernels of orders 2, 4, 6 have been implemented.")
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
  x <- abs(x)
  if (!convolution) {
    abc <- if (order == 2) c(1, 0, 0) else if (order == 4) {
      switch(kernel,
             uniform = c(NA, NA, NA),
             triangular = c(12, -30, 0) / 7,
             epanechnikov = c(15, -35, 0) / 8,
             quartic = c(7, -21, 0) / 4,
             gaussian = c(3, -1, 0) / 2
      )
    } else if (order == 6) {
      switch(kernel,
             uniform = c(NA, NA, NA),
             triangular = c(1635, -10500, 11970) / 683,
             epanechnikov = c(175, -1050, 1155) / 64,
             quartic = c(315, -2310, 3003) / 128,
             gaussian = c(15, -10, 1) / 8
      )
    }
    k <- switch(kernel,
      uniform = 1/2 * (x < 1) * (order == 2) + (-1/6 * (x < 1) + 4/3 * (x < 0.5)) * (order == 4) + (1/20 * (x < 1) - 9/20 * (x < 2/3) + 9/4 * (x < 1/3)) * (order == 6),
      triangular = (1 - x) * (x < 1),
      epanechnikov = 3/4 * (1 - x^2) * (x < 1),
      quartic = 15/16 * (1 - x^2)^2 * (x < 1),
      gaussian = stats::dnorm(x)
    )
    if (order > 2 & kernel != "uniform") {
      polynomial <- abc[1] + abc[2] * x^2 + abc[3] * x^4
      k <- k * polynomial
    }
  } else { # Convolution kernel
    if (order == 2) {
      k <- switch(kernel,
                  uniform = 1 / 4 * (2 - x) * (x < 2),
                  triangular = 1 / 6 * ((3 * x^3 - 6 * x^2 + 4) * (x <= 1) + (8 - 12 * x + 6 * x^2 - x^3) * (x > 1 & x < 2)),
                  epanechnikov = 3 / 160 * (2 - x)^3 * (x^2 + 6*x + 4) * (x < 2),
                  quartic = 5 / 3584 * (2 - x)^5 * (16 + (2*x + x^2)*(20 + 8*x + x^2)) * (x < 2),
                  gaussian = stats::dnorm(x, sd = sqrt(2))
      )
    } else if (order == 4) {
    k <- switch(kernel,
                uniform = 1/36*(64*(1 - x)*(x < 1) - 16*(1*(2*x < 1) + (3/2 - x)*(1/2 <= x & x < 3/2)) + (2 - x)*(x<2)),
                triangular = 3/343 * ((152 + x^2*(-616 + x*(238 + x*(560 + x*(-322 + 5*x*(-14 + 9*x))))))*(x<1) + (2-x)^3 * (30 + x*(-4 + x*(-74 + 5*x*(4 + 3*x))))*(x>=1 & x<2)),
                epanechnikov = 1/2048 * (5*(2-x)^3*(64 + x*(96 + x*(-144 + x*(-160 + x*(48 + 7*x*(6 + x))))))) * (x < 2),
                quartic = 1/2342912 * (35*(2 - x)^5*(2944 + x*(7360 + x*(-1440 + x*(-18320 + x*(-9896 + 9*x*(560 + x*(536 + 15*x*(10 + x))))))))) * (x < 2),
                gaussian = 1/64 * stats::dnorm(x, sd = sqrt(2)) * (108 - 28*x^2 + x^4)
    )
  } else if (order == 6) {
    k <- switch(kernel,
                uniform = 1/400*(956 - 2107*x)*(x<1/3) + 1/400*(680 - 1279*x)*(x>=1/3 & x<2/3) + (-61/40 + 41/25*x)*(x>=2/3 & x<1) + (1/2 - 77/200*x)*(x>=1 & x<4/3) + 1/400*(-28 + 17*x)*(x>=4/3 & x<5/3) + (2-x)/400*(x>=5/3 & x<2),
                triangular = 1/466489*((15/22)*(1353180 + x^2*(-10733690 + x*(3663605 + 2*x*(11535524 + x*(-6217288 + x*(-8435812 + 3*x*(1846130 + 133*x*(4400 + x*(-3168 + 19*x*(-22 + 15*x))))))))))*(x<1) + (-(15/22))*(-2 + x)^3*(240311 + 2*x*(-24496 + x*(-770780 + x*(257608 + x*(1018338 + 133*x*(-3072 + x*(-2180 + 57*x*(8 + 5*x))))))))*(x>=1 & x<2)),
                epanechnikov = -((105*(-61440 + x^2*(465920 + x*(-232960 + x*(-838656 + x*(640640 + x*(439296 - 486720*x + 84448*x^3 - 9828*x^5 + 495*x^7)))))))/3407872) * (x < 2),
                quartic = (-((315*(-2 + x)^5*(380928 + x*(952320 + x*(-1739776 + x*(-6254080 + x*(478464 + x*(9024512 + x*(2918912 + x*(-3982464 + x*(-2272576 + 143*x*(1000 + x*(3012 + 91*x*(10 + x)))))))))))))/1853882368)) * (x < 2),
                gaussian = stats::dnorm(x, sd = sqrt(2)) * (36240 - 19360*x^2 + 2312*x^4 - 88*x^6 + x^8)/16384
    )
  }
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
#' @param order Passed to \code{kernelFun}.
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
                          order = 2,
                          rescale = TRUE,
                          convolution = FALSE
) {
  if (is.null(xgrid)) xgrid <- x # If no grid was passed, use existing data points as the grid
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (one.dim) {
    if (length(bw) > 1) stop("For one-dimensional kernel weights, the bandwidth must be a scalar.")
    diffs <- outer(x, xgrid, "-")
    PK <- kernelFun(diffs / bw, kernel = kernel, order = order, rescale = rescale, convolution = convolution)
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
      K <- kernelFun(outer(x[, i], xgrid[, i], "-") / bw[i], kernel = kernel, order = order, rescale = rescale, convolution = convolution)
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
#' @param kernel Passed to \code{kernelWeights}.
#' @param order Passed to \code{kernelWeights}.
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
                          order = 2,
                          rescale = TRUE
) {
  if (is.null(bw)) {
    bw <- bw.rot(x)
    warning(paste0("No bandwidth supplied, using Silverman's rule of thumb: bw = ", paste0(round(bw, 5), collapse = ", "), "."))
  }
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (!one.dim & length(bw) == 1) bw <- rep(bw, ncol(x))
  K <- kernelWeights(x = x, xgrid = xgrid, bw = bw, kernel = kernel, order = order, rescale = rescale)
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
                         order = 2,
                         rescale = TRUE,
                         LOO = FALSE
) {
  if (is.null(bw)) {
    bw <- bw.rot(x)
    warning(paste0("No bandwidth supplied, using Silverman's rule of thumb: bw = ", paste0(round(bw, 5), collapse = ", "), "."))
  }
  if (is.numeric(xgrid)) {
    if (LOO & !isTRUE(all.equal(x, xgrid))) stop("The Leave-one-out estimator must use the same xgrid as x or NULL.")
  }
  K <- kernelWeights(x = x, xgrid = xgrid, bw = bw, kernel = kernel, order = order, rescale = rescale)
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
#' @param na.rm Logical: should missing values be removed? Setting it to TRUE may cause issued because variable-wise removal of NAs may return a bandwidth that is inappropriate for the final data set to which it is applied.
#' @return A bandwidth that might be optimal for non-parametric density estimation of \code{x}.
#' @examples
#' set.seed(1); bw.rot(rnorm(100)) # Should be 0.3787568 in R version 4.0.4
#' set.seed(1); bw.rot(matrix(rnorm(500), ncol = 10)) # 0.4737872 ... 0.7089850
#' @export
bw.rot <- function(x, na.rm = FALSE) {
  if (any(is.na(x))) {
    if (na.rm) warning("bw.rot: There are missing values in the data, and you should do something about it because proper analysis is impossible with NA, and your results might be unreliable with these bandwidths.") else
      stop("bw.rot: There are missing values in the data, but non-parametric methods rely on data with finite numeric values only.")
  }
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (one.dim) x <- matrix(x, ncol = 1)
  d <- ncol(x)
  n <- nrow(x)
  s <- apply(x, 2, function(x) stats::sd(x, na.rm = na.rm))
  AK <- (4 / (2*d + 1))^(1 / (d + 4)) # (4.15) from Silverman (1986)
  if (any(!is.finite(s))) {
    stop("bw.rot: Could not compute the bandwidth, check your data, most likely it has length 1.")
  } else if (all(s > 0)) {
    return(AK * s * length(x)^(-1/(d+4)))
  } else {
    return(rep(1, d))
  }
}


#' Density and Least-Squares cross-validation criteria
#'
#' Density cross-validation is performed based purely on the values of \code{x} for the kernel density estimator.
#' Least-squares cross-validation is performed for the Nadaraya-Watson estimator of E(Y|X).
#'
#' @param x A numeric vector or numeric matrix.
#' @param bw A numeric scalar, vector, or matrix of candidate bandwidths.
#' @param same Logical: use the same bandwidth in all dimensions for multi-variate x?
#' @param kernel Which kernel to use? Passed to \code{kernelWeights}.
#' @param order Passed to \code{kernelWeights}.
#' @param rescale Passed to \code{kernelWeights}.
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
                same = FALSE,
                kernel = "gaussian",
                order = 2,
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
    KK <- kernelWeights(x, x, bw = b, kernel = kernel, order = order, rescale = rescale, convolution = TRUE)
    term1 <- sum(KK) / (n^2 * prod(b))
    # Computing the LOO estimator efficiently: fhat_i(x) = n/(n-1) * fhat(x) - 1/((n-1)*b^s) * K((X[i] - x)/b)
    fhat <- kernelDensity(x, bw = b, kernel = kernel, rescale = rescale)
    fhat.LOO <- (n * fhat - kernelFun(0, kernel = kernel, order = order, rescale = rescale)^length(b) / (prod(b))) / (n - 1)
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
                 same = FALSE,
                 kernel = "gaussian",
                 order = 2,
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
    muhat_i <- kernelSmooth(x = x, y = y, bw = b, kernel = kernel, order = order, rescale = rescale, LOO = TRUE)
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
#' @param order Passed to \code{kernelWeights}.
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
#' pars <- expand.grid(k = c("triangular", "gaussian"), o = c(2, 4, 6), stringsAsFactors = FALSE)
#' DCV.vals <- lapply(1:6, function(i) DCV(x, bw.grid, kernel = pars$k[i], order = pars$o[i]))
#' LSCV.vals <- lapply(1:6, function(i) LSCV(x, y, bw.grid, kernel = pars$k[i], order = pars$o[i]))
#' DCV.vals <- lapply(DCV.vals, function(x) (x - min(x))/diff(quantile(x, c(0, 0.9))))
#' LSCV.vals <- lapply(LSCV.vals, function(x) (x - min(x)) / diff(quantile(x, c(0, 0.80))))
#'
#' par(mfrow = c(1, 2))
#' plot(NULL, NULL, bty = "n", xlab = "Bandwidth", ylab = "DCV", xlim=c(0.13, 1), ylim=c(0, 1))
#' lapply(seq_along(DCV.vals), function(i) {
#'   lines(bw.grid, DCV.vals[[i]], col = round(i/2 + 0.1), lty = 2 - i %% 2)
#'   points(bw.grid[which.min(DCV.vals[[i]])], min(DCV.vals[[i]]), col = ceiling(i/2), pch = 16)
#' })
#' abline(v = bw.rot(x), lty = 3)
#' legend("top", c("triangular", "gaussian"), lty = 1:2, bty = "n")
#' plot(NULL, NULL, bty = "n", xlab = "Bandwidth", ylab = "LSCV", xlim=c(0.1, 1), ylim=c(0, 1))
#' lapply(seq_along(LSCV.vals), function(i) {
#'   lines(bw.grid, LSCV.vals[[i]], col = ceiling(i/2), lty = 2 - i %% 2)
#'   points(bw.grid[which.min(LSCV.vals[[i]])], min(LSCV.vals[[i]]), col = ceiling(i/2), pch = 16)
#' })
#' abline(v = bw.rot(x), lty = 3)
#' legend("top", paste0("Order ", c(2, 4, 6)), col = 1:3, lwd = 1, bty = "n")

#' bw.CV(x) # The easiest invocation
#' bw.CV(x, opt.fun = "optimise", interval = c(0.1, 1))
#' bw.CV(x, opt.fun = "nlminb")
#' bw.CV(x, kernel = "triangular", order = 4, start.bw = 0.5, opt.fun = "optim", method = "BFGS")
#' bw.CV(x, y = y, kernel = "triangular", order = 4, CV = "LSCV")
#' bw.CV(x, y = y, kernel = "triangular", order = 4, CV = "LSCV", start.bw = 0.2) # Unlucky start
#' bw.CV(x, y = y, kernel = "triangular", order = 4, CV = "LSCV", start.bw = 0.82) # Unlucky start
#' bw.CV(x, y = y, kernel = "triangular", order = 4, CV = "LSCV", start.bw = 0.9) # Unlucky start
bw.CV <- function(x, y = NULL, kernel = "gaussian", order = 2, start.bw = NULL, same = FALSE, CV = c("DCV", "LSCV"),
                  opt.fun = c("nlm", "optim", "nlminb", "optimise"),
                  ret.fun = NULL,
                  par.name.in.opt = NULL,
                  fun.name.in.opt = NULL, ...) {
  opt.fun <- opt.fun[1]
  CV <- CV[1]
  one.dim <- is.vector(x) # Are our data one-dimensional?
  opt <- get(opt.fun)
  f.to.min <- if (CV == "DCV") function(b) DCV(x = x, bw = b, kernel = kernel, order = order, same = same) else if (CV == "LSCV") function(b) LSCV(x = x, y = y, bw = b, kernel = kernel, order = order, same = same) else stop("bw.CV: 'CV' should be either 'DCV' or 'LSCV'!")
  if (is.null(ret.fun)) ret.fun <- switch(opt.fun, nlm = function(x) x[["estimate"]],
                                          optim = function(x) x[["par"]],
                                          optimise = function(x) x[["minimum"]],
                                          nlminb = function(x) x[["par"]])
  if (is.null(start.bw)) start.bw <- bw.rot(x)
  if (same) start.bw <- mean(start.bw)
  opt.result <- switch(opt.fun,
                       nlm = suppressWarnings(stats::nlm(f = f.to.min, p = start.bw, ...)),
                       optim = stats::optim(par = start.bw, fn = f.to.min, ...),
                       optimise = stats::optimise(f = f.to.min, ...),
                       nlminb = stats::nlminb(start = start.bw, objective = f.to.min, lower = 1e-8, ...))
  switch(opt.fun,
         nlm = message(paste0("nlm exit code ", opt.result$code, ", ||gradient||=", round(sqrt(sum(opt.result$gradient^2)), 6), ", done in ", opt.result$iterations, " iterations.")),
         optim = message(paste0("optim exit code ", opt.result$convergence, ", done in (", paste0(opt.result$counts, collapse = ", "), ") iterations.")),
         optimise = message(paste0("optimise does not return useful convergence information.")),
         nlminb = message(paste0("nlminb exit code ", opt.result$convergence, ", done in ", paste0(opt.result$iterations, collapse = ", "), " iterations.")))
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
#' @param bw Bandwidth: a scalar, a vector of the same length as \code{ncol(x)} or a list of the same length as \code{unique(sort(by))}.
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
                               bw = NULL,
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
  by.vals <- sort(unique(by))
  if (!is.list(bw)) bw <- lapply(seq_along(by.vals), function(x) bw)
  for (i in seq_along(by.vals)) {
    s <- by == by.vals[i]
    g <- bygrid == by.vals[i]
    prob <- sum(s) / n
    x.sub <- if (one.dim) x[s] else x[s, ]
    xgrid.sub <- if (one.dim) xgrid[g] else xgrid[g, ]
    res[g] <- kernelDensity(x = x.sub, xgrid = xgrid.sub, bw = bw[[i]], ...) * prob
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

#' Robinson’s semi-parametric regression
#'
#' @param xnames A character vector of exogenous regressors entering the linear part of the model.
#' @param znames A character vector of exogenous regressors entering the non-parametric part of the model. Must be continuous variables (an assumption made by Robinson)
#' @param yname A character scalar: dependent variable name.
#' @param data A data frame containing the variables specified in \code{xnames}, \code{znames}, \code{yname}.
#' @param bw.x A list containing bandwidth vectors or scalars for smoothing X1(Z), X2(Z), ...
#' @param bw.y A numeric vector or scalar for smoothing Y(Z).
#' @param order Non-parametric smoothing order. Robinson (1988) inststs that \code{order >= length(znames)}.
#' @param ... Passed to \code{kernelSmooth}.
#'
#' @return Returns an object of class "lm", like an ordinary \code{lm()}, would, with an extra component: \code{non.parametric} for the estimated mu(Z).
#' @importFrom stats lm
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 1000
#' age <- rchisq(n, 5)
#' female <- rbinom(n, 1, 0.5)
#' muscle <- rnorm(n, sd = 3) - female
#' knowledge <- runif(n, 1, 8) + age/2 + female - muscle/4
#' mu <- function(x1, x2) 0.5*x1 + 0.5*x2 + 0.02*x1*x2 + 2*sin(x1) + 2*cos(0.5*x2)
#' wage <- 1 + age + female + age*female + mu(muscle, knowledge) + rnorm(n, sd = 2)
#' d <- data.frame(cbind(wage, age, female, age.female = age*female, muscle, knowledge))
#' diag(var(d))
#' plot(d)
#' x1grid <- quantile(muscle, 1:49/50)
#' x2grid <- quantile(knowledge, 1:49/50)
#' persp(x1grid, x2grid, outer(x1grid, x2grid, mu), theta = 45, phi = 30, ticktype = "detailed")
#'
#' mod.linear <- lm(wage ~ age + female + age*female + muscle + I(muscle^2) + knowledge +
#'   I(knowledge^2), data = d)
#' mod.semipar <- semiParLM(xnames = c("age", "female", "age.female"),
#'   znames = c("muscle", "knowledge"), yname = "wage", data = d)
#'
#' plot(d$age - mod.semipar$predicted.X[, 1], d$wage - mod.semipar$predicted.Y,
#'   col = c("blue", "red")[d$female + 1])
#' plot(scale(d$muscle, scale = FALSE), scale(mu(muscle, knowledge), scale = FALSE))
#' points(scale(d$muscle, scale = FALSE), scale(mod.semipar$non.parametric, scale = FALSE),
#'   col = "red", pch = 16, cex = 0.5)
#' plot(scale(d$knowledge, scale = FALSE),   scale(mu(muscle, knowledge),      scale = FALSE))
#' points(scale(d$knowledge, scale = FALSE), scale(mod.semipar$non.parametric, scale = FALSE),
#'   col = "red", pch = 16, cex = 0.5)
semiParLM <- function(xnames, znames, yname, data, bw.x = NULL, bw.y = NULL, order = 4, ...) {
  Yhat <- kernelSmooth(x = as.matrix(data[, znames]), y = as.numeric(data[, yname]), bw = bw.y, order = order, ...)
  Xhat <- do.call(cbind, lapply(seq_along(xnames), function(i) kernelSmooth(x = as.matrix(data[, znames]), y = as.numeric(data[, xnames[i]]), bw = bw.x[[i]], order = order, ...)))
  Yc <- data[, yname] - Yhat
  one.dim.x <- length(xnames) == 1
  if (one.dim.x) Xhat <- as.numeric(Xhat)
  Xc <- as.matrix(data[, xnames]) - Xhat
  if (one.dim.x) Xc <- as.numeric(Xc)
  d <- data.frame(Xc, y = Yc)
  colnames(d)[1:length(xnames)] <- xnames
  fl <- paste0("y ~ ", paste0(xnames, collapse = " + "))
  m <- lm(fl, data = d)
  m$predicted.Y <- Yhat
  m$predicted.X <- Xhat
  colnames(m$predicted.X) <- xnames
  if (!one.dim.x) m$predicted.X <- as.data.frame(m$predicted.X)
  m$non.parametric <- Yhat - as.numeric(cbind(1, Xhat) %*% m$coefficients)
  m$bw.x <- if (!is.null(bw.x)) bw.x else bw.rot(data[, znames])
  m$bw.y <- if (!is.null(bw.y)) bw.x else lapply(1:length(xnames), function(i) bw.rot(data[, znames]))
  return(m)
}

