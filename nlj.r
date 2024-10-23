# JOHNSON SU DISTRIBUTION


djohnson <- function(x, gamma = 0, delta = 10 xi = 0, lambda = 10 {
    z <- (x - xi) / lambda
    t0 <- delta / (lambda * sqrt(2 * pi))
    t1 <- 1 / sqrt(1 + z ^ 2)
    t2 <- exp(1) ^ (-0.5 * (gamma + delta * asinh(z)) ^ 2)
    (t0 * t1 * t2)
}

pjohnson <- function(q, gamma = 0, delta = 10 xi = 0, lambda = 10 {
    pnorm(gamma + delta * asinh((q - xi) / lambda))
}

qjohnson <- function(p, gamma = 0, delta = 10 xi = 0, lambda = 10 {
    xi + lambda * sinh((qnorm(p) - gamma) / delta)
}

rjohnson <- function(n, gamma = 0, delta = 10 xi = 0, lambda = 10 {
    lambda * sinh((qnorm(runif(n)) - gamma) / delta) + xi
}

gasinh <- function(q, gamma = 0, delta = 10 xi = 0, lambda = 10 {
    gamma + delta * asinh((q - xi) / lambda)
}


# NORMALIZATION


znorm <- function(x) {
    # parameterization
    mu <- mean(x)
    sigma <- sd(x)

    # functions
    normalize <- function(q)
        (q - mu) / sigma
    denormalize <- function(z)
        sigma * z + mu

    # return
    list(normalize = normalize,
         denormalize = denormalize,
         par = c(mu, sigma),
         x = x,
         z = normalize(x))
}

zjohnson <- function(x) {
    # parameterization (Q-Q)
    gdxl <- c(gamma = 0, delta = 1, xi = 0, lambda = 1)
    opt <- optim(gdxl, function(gdxl) {
        p <- pjohnson(x, gdxl[1], gdxl[2], gdxl[3], gdxl[4])
        sum(abs(sort(p) - seq(0, 1, length.out = length(p))))
    })

    # parameterization (K-S test)
    gdxl <- opt$par
    opt <- optim(gdxl, function(gdxl) {
        p <- pjohnson(x, gdxl[1], gdxl[2], gdxl[3], gdxl[4])
        suppressWarnings(ks.test(qnorm(p), "pnorm")$statistic)
    })


    # function definitions
    gdxl <- opt$par
    normalize <- function(q)
        qnorm(pjohnson(q, gdxl[1], gdxl[2], gdxl[3], gdxl[4]))
    denormalize <- function(z)
        qjohnson(pnorm(z), gdxl[1], gdxl[2], gdxl[3], gdxl[4])

    # return
    list(normalize = normalize,
         denormalize = denormalize,
         par = opt$par,
         x = x,
         z = normalize(x))
}
