# JOHNSON SU DISTRIBUTION


djohnson <- function(x, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    z <- (x - xi) / lambda
    t0 <- delta / (lambda * sqrt(2 * pi))
    t1 <- 1 / sqrt(1 + z ^ 2)
    t2 <- exp(1) ^ (-0.5 * (gamma + delta * asinh(z)) ^ 2)
    (t0 * t1 * t2)
}

pjohnson <- function(q, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    pnorm(gamma + delta * asinh((q - xi) / lambda))
}

qjohnson <- function(p, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    xi + lambda * sinh((qnorm(p) - gamma) / delta)
}

rjohnson <- function(n, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    lambda * sinh((qnorm(runif(n)) - gamma) / delta) + xi
}

gasinh <- function(q, gamma = 0, delta = 1, xi = 0, lambda = 1) {
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


# JOINT NORMALIZATION


# This process optimizes the linear relationship in a bivariate setting through joint normalization with the Johnson SU distribution.
# lm.johnson <- function(x, y, iterations = 3) {
#     prm <- c(0, 1, 0, 1,
#              0, 1, 0, 1)

#     for (k in 1:iterations) {
#         opt <- optim(prm, function(prm) {
#             gdxl_x <- prm[1:4]
#             gdxl_y <- prm[5:8]

#             xn <- qnorm(pjohnson(x, gdxl_x[1], gdxl_x[2], gdxl_x[3], gdxl_x[4]))
#             yn <- qnorm(pjohnson(y, gdxl_y[1], gdxl_y[2], gdxl_y[3], gdxl_y[4]))

#             tryCatch({
#                 m <- lm(yn ~ xn)
#                 suppressWarnings(ks.test(m$residuals, "pnorm")$p.value)
#             }, error = function(e) {
#                 Inf
#             })
#         })

#         prm <- opt$par
#     }

#     gdxl_x <- prm[1:4]
#     gdxl_y <- prm[5:8]

#     xn <- qnorm(pjohnson(x, gdxl_x[1], gdxl_x[2], gdxl_x[3], gdxl_x[4]))
#     yn <- qnorm(pjohnson(y, gdxl_y[1], gdxl_y[2], gdxl_y[3], gdxl_y[4]))

#     # xi <- seq(min(xn), max(xn), length.out = 300)
#     xi <- xn
#     m <- lm(yn ~ xn)
#     print(summary(m))
#     yi <- predict(m, data.frame(xn = xi))

#     xq <- qjohnson(pnorm(xi), gdxl_x[1], gdxl_x[2], gdxl_x[3], gdxl_x[4])
#     yq <- qjohnson(pnorm(yi), gdxl_y[1], gdxl_y[2], gdxl_y[3], gdxl_y[4])

#     list(x = xq, y = yq)
# }


# TESTING


df <- mtcars[order(mtcars$hp), ]
str(df)

{
    attach(df)

    jhp <- zjohnson(hp)
    jwt <- zjohnson(wt)
    jmpg <- zjohnson(mpg)

    summary(lm(mpg ~ hp))
    summary(lm(log(mpg) ~ log(hp)))
    summary(lm(jmpg$z ~ jhp$z))

    summary(lm(mpg ~ hp + wt))
    summary(lm(log(mpg) ~ log(hp) + log(wt)))
    summary(lm(jmpg$z ~ jhp$z + jwt$z))

    summary(lm(mpg ~ hp + wt + disp + drat + qsec))

    # setup
    response <- qsec
    predictors <- list(mpg, hp, wt)

    # initialization
    init <- FALSE
    if (init) {
        prm <- unname(c(zjohnson(response)$par, unlist(lapply(predictors, function(x) zjohnson(x)$par))))
    } else {
        prm <- rep(c(0, 1, 0, 1), 1 + length(predictors))
    }

    # optimization
    last <- NULL
    better <- TRUE
    while (better) {
        opt <- optim(prm, function(prm) {
            y <- gasinh(response, prm[1], prm[2], prm[3], prm[4])
            x <- rowSums(sapply(1:length(predictors), function(k) {
                gasinh(predictors[[k]], prm[4 * k + 1], prm[4 * k + 2], prm[4 * k + 3], prm[4 * k + 4])
            }))
            tryCatch({
                m <- lm(y ~ x)
                -summary(m)$fstatistic[1]
            }, error = function(e) {
                Inf
            })
        })

        prm <- opt$par
        print(c(last / opt$value, opt$value, last))
        if (is.null(last)) {
            last <- opt$value
        } else if (0.999 < (last / opt$value)) {
            better <- FALSE
        }
        last <- opt$value
    }

    y <- m$fitted.values
    # y <- gamma + delta * asinh((x - xi) / lambda)
    y_hat <- sinh((m$fitted.values - prm[1]) / prm[2]) * prm[4] + prm[3]
    ratio <- response / y_hat
    round(cbind(mtcars[, c(1, 4, 6, 7)], y_hat = y_hat, ratio = ratio), 2)[order(ratio), ]

}

with(mtcars, {
    pdf("nlj.pdf", 9, 6)

    plot(mpg, qsec, pch = 16)
    segments(mpg, qsec, mpg, y_hat, pch = 16, col = "red")
    text(mpg, qsec, 1:nrow(mtcars), pos = 3)

    plot(hp, qsec, pch = 16)
    segments(hp, qsec, hp, y_hat, pch = 16, col = "red")
    text(hp, qsec, 1:nrow(mtcars), pos = 3)

    plot(wt, qsec, pch = 16)
    segments(wt, qsec, wt, y_hat, pch = 16, col = "red")
    text(wt, qsec, 1:nrow(mtcars), pos = 3)

    dev.off()
})

with(df, {
    pdf("nlj.pdf", 9, 6)

    plot(hp, mpg, pch = 16)
    lines(hp, lm(mpg ~ hp)$fitted.values, col = "red")

    plot(hp, mpg, pch = 16)
    lines(hp, exp(lm(log(mpg) ~ log(hp))$fitted.values), col = "red")

    jhp <- zjohnson(hp)
    jmpg <- zjohnson(mpg)

    plot(hp, mpg, pch = 16)
    lines(hp, jmpg$denormalize(lm(jmpg$z ~ jhp$z)$fitted.values), col = "red")

    dev.off()
})

{
    pdf("nlj.pdf", 9, 9)

    y <- qnorm(pjohnson(response, prm[1], prm[2], prm[3], prm[4]))
    x <- sapply(1:length(predictors), function(k) {
        qnorm(pjohnson(predictors[[k]], prm[4 * k + 1], prm[4 * k + 2], prm[4 * k + 3], prm[4 * k + 4]))
    })
    mat <- as.data.frame(cbind(y, x))
    m <- lm(y ~ rowSums(x))
    summary(m)

    plot(mat, pch = 16)
    pnorm(y)
    qjohnson(pnorm(y), prm[1], prm[2], prm[3], prm[4])
    mpg_hat <- qjohnson(pnorm(m$fitted.values), prm[1], prm[2], prm[3], prm[4])
    rownames(df)[order(mpg / mpg_hat)]

    dev.off()
}
