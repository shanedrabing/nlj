# GENERAL


# :-)


#' Sum of Log-Cosh Loss
#'
#' Computes the sum of the log-cosh of residuals, a smooth approximation to absolute error loss.
#'
#' @param m A linear model object, which should contain residuals.
#' @return A single numeric value representing the sum of the log-cosh of the residuals.
#' @export
sum_log_cosh <- function(m) {
    sum(log(cosh(m$residuals)))
}

#' Total Absolute Deviation from Theoretical CDF
#'
#' Calculates the total absolute deviation from a theoretical cumulative distribution function (CDF).
#'
#' @param p A numeric vector representing sorted probabilities.
#' @return A single numeric value representing the total absolute deviation from the theoretical CDF.
#' @export
qqtad <- function(p) {
    sum(abs(p - seq(0, 1, length.out = length(p))))
}

#' Generalized Inverse Hyperbolic Sine Transformation
#'
#' Computes the generalized inverse hyperbolic sine transformation.
#'
#' @param x A numeric vector to transform.
#' @param xi Numeric, location parameter (default = 0).
#' @param lambda Numeric, scale parameter (default = 1).
#' @return A numeric vector of transformed values.
#' @export
gasinh <- function(x, xi = 0, lambda = 1) {
    asinh((x - xi) / lambda)
}

#' Generalized Hyperbolic Sine Transformation
#'
#' Computes the generalized hyperbolic sine transformation (inverse of the asinh transformation).
#'
#' @param x A numeric vector to transform back.
#' @param xi Numeric, location parameter (default = 0).
#' @param lambda Numeric, scale parameter (default = 1).
#' @return A numeric vector of back-transformed values.
#' @export
gsinh <- function(x, xi = 0, lambda = 1) {
    sinh(x) * lambda + xi
}


# JOHNSON SU DISTRIBUTION


#' Density of Johnson SU Distribution
#'
#' Computes the density of the Johnson SU distribution.
#'
#' @param x A numeric vector of values at which to evaluate the density.
#' @param gamma Numeric, shape parameter (default = 0).
#' @param delta Numeric, shape parameter (default = 1).
#' @param xi Numeric, location parameter (default = 0).
#' @param lambda Numeric, scale parameter (default = 1).
#' @return A numeric vector of density values for the Johnson SU distribution.
#' @export
djohnson <- function(x, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    z <- (x - xi) / lambda
    scale <- delta / (lambda * sqrt(2 * pi))
    shape <- 1 / sqrt(1 + z ^ 2)
    tail <- exp(-0.5 * (gamma + delta * asinh(z)) ^ 2)
    scale * shape * tail
}

#' Cumulative Distribution of Johnson SU
#'
#' Computes the cumulative distribution function of the Johnson SU distribution.
#'
#' @param q A numeric vector of quantiles.
#' @param gamma Numeric, shape parameter (default = 0).
#' @param delta Numeric, shape parameter (default = 1).
#' @param xi Numeric, location parameter (default = 0).
#' @param lambda Numeric, scale parameter (default = 1).
#' @return A numeric vector of cumulative probabilities for the Johnson SU distribution.
#' @export
pjohnson <- function(q, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    stats::pnorm(gamma + delta * asinh((q - xi) / lambda))
}

#' Inverse Cumulative Distribution of Johnson SU
#'
#' Computes the inverse cumulative distribution (quantile function) of the Johnson SU distribution.
#'
#' @param p A numeric vector of probabilities.
#' @param gamma Numeric, shape parameter (default = 0).
#' @param delta Numeric, shape parameter (default = 1).
#' @param xi Numeric, location parameter (default = 0).
#' @param lambda Numeric, scale parameter (default = 1).
#' @return A numeric vector of quantiles for the Johnson SU distribution.
#' @export
qjohnson <- function(p, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    xi + lambda * sinh((stats::qnorm(p) - gamma) / delta)
}

#' Random Deviates from Johnson SU Distribution
#'
#' Generates random deviates from the Johnson SU distribution.
#'
#' @param n Integer, the number of random values to generate.
#' @param gamma Numeric, shape parameter (default = 0).
#' @param delta Numeric, shape parameter (default = 1).
#' @param xi Numeric, location parameter (default = 0).
#' @param lambda Numeric, scale parameter (default = 1).
#' @return A numeric vector of random deviates from the Johnson SU distribution.
#' @export
rjohnson <- function(n, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    xi + lambda * sinh((stats::qnorm(stats::runif(n)) - gamma) / delta)
}


# NORMALIZATION


#' Z-Score Normalization
#'
#' Normalizes data using Z-score normalization.
#'
#' @param x A numeric vector to normalize.
#' @return A list containing:
#'   \describe{
#'     \item{normalize}{A function that applies the normalization to a new vector.}
#'     \item{denormalize}{A function that reverses the normalization.}
#'     \item{par}{A numeric vector of parameters (mean and standard deviation).}
#'     \item{fd}{A function to compute the density under Z-normalization.}
#'     \item{fp}{A function to compute cumulative probability under Z-normalization.}
#'     \item{fq}{A function to compute quantiles under Z-normalization.}
#'     \item{fr}{A function to generate random deviates under Z-normalization.}
#'     \item{x}{The original input data.}
#'     \item{z}{The normalized data.}
#'   }
#' @export
znorm <- function(x) {
    mu <- mean(x)
    sigma <- stats::sd(x)

    normalize <- function(q) (q - mu) / sigma
    denormalize <- function(z) z * sigma + mu
    fd <- function(x) stats::dnorm(x, mu, sigma)
    fp <- function(q) stats::pnorm(q, mu, sigma)
    fq <- function(p) stats::qnorm(p, mu, sigma)
    fr <- function(n) stats::rnorm(n, mu, sigma)

    list(normalize = normalize,
         denormalize = denormalize,
         par = c(mu, sigma),
         fd = fd,
         fp = fp,
         fq = fq,
         fr = fr,
         x = x,
         z = normalize(x))
}

#' Johnson SU Normalization
#'
#' Normalizes data using the Johnson SU distribution.
#'
#' @param x A numeric vector to normalize.
#' @return A list containing:
#'   \describe{
#'     \item{normalize}{A function that applies the Johnson SU normalization.}
#'     \item{denormalize}{A function that reverses the normalization.}
#'     \item{par}{A numeric vector of optimized Johnson SU parameters.}
#'     \item{fd}{A function to compute density under Johnson SU normalization.}
#'     \item{fp}{A function to compute cumulative probability under Johnson SU normalization.}
#'     \item{fq}{A function to compute quantiles under Johnson SU normalization.}
#'     \item{fr}{A function to generate random deviates under Johnson SU normalization.}
#'     \item{x}{The original input data.}
#'     \item{z}{The normalized data.}
#'   }
#' @export
zjohnson <- function(x) {
    gdxl <- c(gamma = 0, delta = 1, xi = 0, lambda = 1)
    opt <- stats::optim(gdxl, function(gdxl) {
        p <- pjohnson(x, gdxl[1], gdxl[2], gdxl[3], gdxl[4])
        qqtad(p[order(x)])
    })

    gdxl <- opt$par
    opt <- stats::optim(gdxl, function(gdxl) {
        p <- pjohnson(x, gdxl[1], gdxl[2], gdxl[3], gdxl[4])
        ks <- suppressWarnings(stats::ks.test(stats::qnorm(p), "pnorm"))
        -ks$p.value
    })

    gdxl <- opt$par
    normalize <- function(q)
        stats::qnorm(pjohnson(q, gdxl[1], gdxl[2], gdxl[3], gdxl[4]))
    denormalize <- function(z)
        qjohnson(stats::pnorm(z), gdxl[1], gdxl[2], gdxl[3], gdxl[4])
    fd <- function(x) djohnson(x, gdxl[1], gdxl[2], gdxl[3], gdxl[4])
    fp <- function(q) pjohnson(q, gdxl[1], gdxl[2], gdxl[3], gdxl[4])
    fq <- function(p) qjohnson(p, gdxl[1], gdxl[2], gdxl[3], gdxl[4])
    fr <- function(n) rjohnson(n, gdxl[1], gdxl[2], gdxl[3], gdxl[4])

    list(normalize = normalize,
         denormalize = denormalize,
         par = gdxl,
         fd = fd,
         fp = fp,
         fq = fq,
         fr = fr,
         x = x,
         z = normalize(x))
}


# GENERALIZED ASINH TRANSFORMATION MODEL (GATM)


#' Generalized Asinh Transformation Model (GATM)
#'
#' Fits a linear model using the generalized asinh transformation.
#'
#' @param formula An object of class `formula` for the model.
#' @param data An optional data frame for variable data.
#' @param loss A function that defines the loss to be minimized (default = `sum_log_cosh`).
#' @param iterations Integer, number of iterations for optimization (default = 1).
#' @param penalty Numeric, regularization penalty (default = 1e-9).
#' @param verbose Logical, if TRUE, prints optimization progress (default = FALSE).
#' @return A list containing:
#'   \describe{
#'     \item{detransform}{A function to reverse the generalized asinh transformation.}
#'     \item{par}{Optimized transformation parameters for each variable.}
#'     \item{fit}{Fitted linear model object.}
#'     \item{opt}{Optimization result object.}
#'     \item{df.sub}{Subset of data used for modeling.}
#'     \item{df.mut}{Transformed data frame used in modeling.}
#'     \item{z}{Detransformed fitted values of the model.}
#'   }
#' @export
lm.gat <- function(formula, data = NULL,
                   loss = sum_log_cosh,
                   iterations = 1, penalty = 1e-9,
                   verbose = FALSE) {

    as.df <- as.data.frame
    use.gasinh <- function(i, df, prm) {
        gasinh(df[, i], prm[2 * (i - 1) + 1], prm[2 * (i - 1) + 2])
    }

    env <- if (is.null(data)) environment(formula) else as.environment(data)
    sub <- as.df(sapply(all.vars(formula), get, envir = env))

    n <- ncol(sub)
    prm <- rep(c(0, 1), n)

    evaluate <- function(prm) {
        mut <- stats::setNames(as.df(sapply(1:n, use.gasinh, sub, prm)), names(sub))

        tryCatch({
            m <- stats::lm(formula, mut)
            cost <- penalty * mean(prm ^ 2)
            loss(m) + cost
        }, error = function(e) {
            Inf
        })
    }

    for (k in 1:iterations) {
        opt <- stats::optim(prm, evaluate)
        prm <- opt$par
        if (verbose) {
            message("Loss: ", format(opt$value, scientific = TRUE))
        }
    }

    prm <- opt$par
    mut <- stats::setNames(as.df(sapply(1:n, use.gasinh, sub, prm)), names(sub))
    m <- stats::lm(formula, mut)

    detransform <- function(y) gsinh(y, prm[1], prm[2])

    list(detransform = detransform,
         par = prm,
         fit = m,
         opt = opt,
         df.sub = sub,
         df.mut = mut,
         z = detransform(m$fitted.values))
}
