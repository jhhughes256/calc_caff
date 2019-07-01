# Functions Utility for Nivolumab Simulation Study
# -----------------------------------------------------------------------------
# CKD-EPI Equation for determining eGFR
  ckdepi.fn <- function(secr, age, sex, black, mgdl = F) {
  # secr: units for equation are mg/dL
    if (mgdl) {
      scr <- secr
    } else {
      scr <- secr/88.4
    }
  # sex: 0 - female, 1 - male
  # black: 0 - non-african american, 1 - african american
  # Determine constants
    if (sex) {
      alpha <- -0.411
      k <- 0.9
    } else {
      alpha <- -0.329
      k <- 0.7
    }
  # Determine eGFR
    141 * min(c(scr/k, 1))^alpha * max(c(scr/k, 1))^-1.209 *   # continues
      0.993^age * 1.018^(1-sex) * 1.159^black
  }
  ckdepi_fn <- ckdepi.fn
  
# BMI equation for use with multivariate sampling output
  bmi.fn <- function(mat, cm = T) {
  # Determine wt and ht objects from matrix
    wt <- mat[1]
    ht <- mat[2]
  # Determine conversion value for height units
    if (cm) { conv <- 100 }
    else { conv <- 1 }
  # Determine BMI
    wt/(ht/conv)^2
  }
  
# BSA equation
  bsa.fn <- function(wt, ht) {
    0.007184*ht^0.725*wt^0.425
  }

# Truncated distribution sampling functions
# Create resampling for distribution if outside the desired range
# Normal Distribution
  trunc.rnorm <- function(n, mean, sd, range, log = F) {
    lower <- min(range)
    upper <- max(range)
    # Convert mean and sd if necessary
    if (log) {  # transform mu and sigma to work in log domain
      mu <- log(mean^2/sqrt((sd^2 + mean^2)))
      sigma <- sqrt(log(sd^2/mean^2 + 1))
      # https://au.mathworks.com/help/stats/lognstat.html
    } else {
      mu <- mean
      sigma <- sd
    }
    # Sample from normal distribution
    out <- rnorm(n, mu, sigma)
    if (log) { out <- exp(out) }
    # Determine values to be resampled
    repeat {
      resample <- out < lower | out > upper
      if (!any(resample)) break  # if none to resample exit loop
      replace <- which(resample)
      new <- rnorm(length(replace), mu, sigma)
      if (log) { new <- exp(new) }
      out[replace] <- new
    }
    out
  }
  trunc_rnorm <- trunc.rnorm

# Multivariate Normal Distribution
  trunc.mvrnorm <- function(n, mean, sd, corr_mat, lower, upper, log = F) {
    require(MASS)
    require(MBESS)
    mat_dim <- length(mean)  # frequently used term
    # Check for errors
    if (length(sd) != mat_dim | length(lower) != mat_dim | length(upper) != mat_dim) {
      stop("mean, sd, lower and upper all must have the same length")
    } else if (!is.matrix(corr_mat)) {
      stop("correlation matrix must be of class matrix")
    } else if (any(!dim(corr_mat) %in% mat_dim)) {
      stop("correlation matrix dimensions must match length of mean")
    }
    # Convert mean and sd if necessary
    if (log) {  # transform mu and sigma to work in log domain
      mu <- log(mean^2/sqrt((sd^2 + mean^2)))
      sigma <- sqrt(log(sd^2/mean^2 + 1))
      # https://au.mathworks.com/help/stats/lognstat.html
    } else {
      mu <- mean
      sigma <- sd
    }
    # Sample from multivariate normal distribution
    cov_mat <- cor2cov(corr_mat, sigma)
    if (log) {  # take normally distributed error terms for transformation
      mvr_mat <- mvrnorm(n, double(mat_dim), cov_mat)
      out <- apply(mvr_mat, 2, exp) %*% diag(mean)
    } else {  # take the normal distribution
      out <- mvrnorm(n, mean, cov_mat)
    }
    # Determine values to be resampled
    repeat {
      # mapply - check if values are outside of their respective range
      # matrix - to coerce from vector to matrix
      # apply - to see which rows have any values that need resampling
      resample <- apply(MARGIN = 1, FUN = any,
        matrix(ncol = mat_dim,
          mapply(function(val, low, upp) {
            val < low | val > upp
          }, out, rep(lower, each = n), rep(upper, each = n))
        )  # matrix
      )  # apply
      if (!any(resample)) break  # if none to resample exit loop
      if (log) {
        mvr_mat <- matrix(ncol = mat_dim,
          mvrnorm(length(which(resample)), double(mat_dim), cov_mat)
        )  # matrix used as length(which(resample)) == 1 returns a vector
        out[resample, ] <- apply(mvr_mat, 2, exp) %*% diag(mean)
      } else {
        out[resample, ] <- mvrnorm(length(which(resample)), mean, cov_mat)
      }
    }
    out
  }
  trunc_mvrnorm <- trunc.mvrnorm
  
# Confidence interval functions
  CI90lo <- function(x) quantile(x, probs = 0.05)
  CI90hi <- function(x) quantile(x, probs = 0.95)

  CI95lo <- function(x) quantile(x, probs = 0.025)
  CI95hi <- function(x) quantile(x, probs = 0.975)