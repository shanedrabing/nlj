# GENERAL


#' Total Absolute Deviation from Theoretical CDF
#'
#' Calculates the total absolute deviation from a theoretical cumulative distribution function (CDF).
#'
#' @param p A numeric vector representing probabilities.
#' @return A single numeric value representing the total absolute deviation.
#' @export
qqtad <- function(p) {
    sum(abs(sort(p) - seq(0, 1, length.out = length(p))))
}

#' Generalized Inverse Hyperbolic Sine Transformation
#'
#' Computes the generalized inverse hyperbolic sine transformation.
#'
#' @param x A numeric vector to transform.
#' @param xi Numeric, location parameter.
#' @param lambda Numeric, scale parameter.
#' @return A numeric vector of transformed values.
#' @export
gasinh <- function(x, xi = 0, lambda = 1) {
    asinh((x - xi) / lambda)
}

#' Inverse of Generalized Inverse Hyperbolic Sine
#'
#' Computes the inverse of the generalized inverse hyperbolic sine transformation.
#'
#' @param y A numeric vector to transform back.
#' @param xi Numeric, location parameter.
#' @param lambda Numeric, scale parameter.
#' @return A numeric vector of back-transformed values.
#' @export
inverse.gasinh <- function(y, xi = 0, lambda = 1) {
    sinh(y) * lambda + xi
}


# JOHNSON SU DISTRIBUTION


#' Density of Johnson SU Distribution
#'
#' Computes the density of the Johnson SU distribution.
#'
#' @param x A numeric vector of values at which to evaluate the density.
#' @param gamma Numeric, shape parameter.
#' @param delta Numeric, shape parameter.
#' @param xi Numeric, location parameter.
#' @param lambda Numeric, scale parameter.
#' @return A numeric vector of density values.
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
#' @param gamma Numeric, shape parameter.
#' @param delta Numeric, shape parameter.
#' @param xi Numeric, location parameter.
#' @param lambda Numeric, scale parameter.
#' @return A numeric vector of cumulative probabilities.
#' @export
pjohnson <- function(q, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    pnorm(gamma + delta * asinh((q - xi) / lambda))
}

#' Inverse Cumulative Distribution of Johnson SU
#'
#' Computes the inverse cumulative distribution (quantile function) of the Johnson SU distribution.
#'
#' @param p A numeric vector of probabilities.
#' @param gamma Numeric, shape parameter.
#' @param delta Numeric, shape parameter.
#' @param xi Numeric, location parameter.
#' @param lambda Numeric, scale parameter.
#' @return A numeric vector of quantiles.
#' @export
qjohnson <- function(p, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    xi + lambda * sinh((qnorm(p) - gamma) / delta)
}

#' Random Deviates from Johnson SU Distribution
#'
#' Generates random deviates from the Johnson SU distribution.
#'
#' @param n Number of random values to generate.
#' @param gamma Numeric, shape parameter.
#' @param delta Numeric, shape parameter.
#' @param xi Numeric, location parameter.
#' @param lambda Numeric, scale parameter.
#' @return A numeric vector of random deviates.
#' @export
rjohnson <- function(n, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    xi + lambda * sinh((qnorm(runif(n)) - gamma) / delta)
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
#'     \item{x}{The original input data.}
#'     \item{z}{The normalized data.}
#'   }
#' @export
znorm <- function(x) {
    # parameterization
    mu <- mean(x)
    sigma <- sd(x)

    # binding
    normalize <- function(q)
        (q - mu) / sigma
    denormalize <- function(z)
        z * sigma + mu

    # return
    list(normalize = normalize,
         denormalize = denormalize,
         par = c(mu, sigma),
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
#'     \item{x}{The original input data.}
#'     \item{z}{The normalized data.}
#'   }
#' @export
zjohnson <- function(x) {
    # parameterization (Q-Q)
    gdxl <- c(gamma = 0, delta = 1, xi = 0, lambda = 1)
    opt <- optim(gdxl, function(gdxl) {
        p <- pjohnson(x, gdxl[1], gdxl[2], gdxl[3], gdxl[4])
        qqtad(p)
    })

    # parameterization (K-S test)
    gdxl <- opt$par
    opt <- optim(gdxl, function(gdxl) {
        p <- pjohnson(x, gdxl[1], gdxl[2], gdxl[3], gdxl[4])
        ks <- suppressWarnings(ks.test(qnorm(p), "pnorm"))
        ks$statistic
    })

    # binding
    gdxl <- opt$par
    normalize <- function(q)
        qnorm(pjohnson(q, gdxl[1], gdxl[2], gdxl[3], gdxl[4]))
    denormalize <- function(z)
        qjohnson(pnorm(z), gdxl[1], gdxl[2], gdxl[3], gdxl[4])

    # return
    list(normalize = normalize,
         denormalize = denormalize,
         par = gdxl,
         x = x,
         z = normalize(x))
}


# GENERALIZED ASINH TRANSFORMATION (GAT) MODEL


#' Generalized Asinh Transformation (GAT) Model
#'
#' Fits a linear model using the generalized asinh transformation.
#'
#' @param formula An object of class `formula`.
#' @param data An optional data frame.
#' @param iterations Number of iterations for optimization.
#' @param penalty Penalty for regularization.
#' @param verbose Logical, if TRUE, prints progress.
#' @return A list containing:
#'   \describe{
#'     \item{predict}{A function to predict new data.}
#'     \item{detransform}{A function to reverse the transformation.}
#'     \item{par}{Optimized transformation parameters.}
#'     \item{fit}{Fitted linear model object.}
#'     \item{opt}{Optimization object.}
#'     \item{df.sub}{Subset of data frame used in modeling.}
#'     \item{df.mut}{Transformed data frame.}
#'     \item{z}{Detransformed fitted values.}
#'   }
#' @export
lm.gat <- function(formula, data = NULL,
                   iterations = 1, penalty = 1e-9,
                   verbose = FALSE) {

    # helper
    as.df <- as.data.frame
    use.gasinh <- function(i, df, prm) {
        gasinh(df[, i], prm[2 * (i - 1) + 1], prm[2 * (i - 1) + 2])
    }

    # subset
    env <- if (is.null(data)) environment(formula) else as.environment(data)
    sub <- as.df(sapply(all.vars(formula), get, envir = env))

    # initialization
    n <- ncol(sub)
    prm <- rep(c(0, 1), n)

    # binding
    evaluate <- function(prm) {
        mut <- setNames(as.df(sapply(1:n, use.gasinh, sub, prm)), names(sub))

        tryCatch({
            m <- lm(formula, mut)
            loss <- sum(log(cosh(m$residuals)))
            cost <- penalty * mean(prm ^ 2)
            loss + cost
        }, error = function(e) {
            Inf
        })
    }

    # optimization
    for (k in 1:iterations) {
        opt <- optim(prm, evaluate)
        prm <- opt$par
        if (verbose) {
            message("Loss: ", format(opt$value, scientific = TRUE))
        }
    }

    # evaluation
    prm <- opt$par
    mut <- setNames(as.df(sapply(1:n, use.gasinh, sub, prm)), names(sub))
    m <- lm(formula, mut)

    # binding
    detransform <- function(y)
        inverse.gasinh(y, prm[1], prm[2])
    predict <- function(data) {
        rhs <- all.names(formula)[-(1:3)]
        sub <- data[, rhs]
        predict(m, sub)
    }

    # return
    list(predict = predict,
         detransform = detransform,
         par = prm,
         fit = m,
         opt = opt,
         df.sub = sub,
         df.mut = mut,
         z = detransform(m$fitted.values))
}
