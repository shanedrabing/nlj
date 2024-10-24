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


lm.gat <- function(resp, preds,
                   iterations = 3, penalty = 1e-9, noise = 0,
                   verbose = FALSE) {

    # initialization
    n <- 1 + length(preds)
    prm <- rep(c(0, 1), n) + rnorm(2 * n, 0, noise)

    # helper
    use.gasinh <- function(i, preds, prm) {
        gasinh(preds[[i]], prm[2 * i + 1], prm[2 * i + 2])
    }

    # binding
    evaluate <- function(prm) {
        y <- gasinh(resp, prm[1], prm[2])
        x <- capply(1:length(preds), use.gasinh, preds, prm)

        tryCatch({
            m <- lm.fit(cbind(1, x), y)
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
    y <- gasinh(resp, prm[1], prm[2])
    x <- capply(1:length(preds), use.gasinh, preds, prm)
    m <- lm.fit(cbind(1, x), y)

    # binding
    detransform <- function(y)
        inverse.gasinh(y, prm[1], prm[2])
    predict <- function(preds) {
        x <- rowSums(capply(1:length(preds), use.gasinh, preds, prm))
        predict(m, data.frame(x))
    }

    # return
    list(predict = predict,
         detransform = detransform,
         par = prm,
         fit = m,
         opt = opt,
         resp = resp,
         preds = preds,
         z = detransform(m$fitted.values))
}
