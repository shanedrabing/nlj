# GENERAL


capply <- function(X, FUN, ...) {
    do.call(cbind, lapply(X, FUN, ...))
}

# total absolute deviation from theoretical CDF
qqtad <- function(p) {
    sum(abs(sort(p) - seq(0, 1, length.out = length(p))))
}

# generalized inverse hyperbolic sine
gasinh <- function(x, xi = 0, lambda = 1) {
    asinh((x - xi) / lambda)
}

inverse.gasinh <- function(y, xi = 0, lambda = 1) {
    sinh(y) * lambda + xi
}


# JOHNSON SU DISTRIBUTION


# density
djohnson <- function(x, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    z <- (x - xi) / lambda
    scale <- delta / (lambda * sqrt(2 * pi))
    shape <- 1 / sqrt(1 + z ^ 2)
    tail <- exp(-0.5 * (gamma + delta * asinh(z)) ^ 2)
    scale * shape * tail
}

# cumulative distribution
pjohnson <- function(q, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    pnorm(gamma + delta * asinh((q - xi) / lambda))
}

# inverse cumulative distribution
qjohnson <- function(p, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    xi + lambda * sinh((qnorm(p) - gamma) / delta)
}

# random deviates
rjohnson <- function(n, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    xi + lambda * sinh((qnorm(runif(n)) - gamma) / delta)
}


# NORMALIZATION


# Z score normalization
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

# normalization using the Johnson SU
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


# GENERALIZED ASINH TRANSFORMATION MODEL (GATM)


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
