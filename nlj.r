# JOHNSON


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
         mean = mu,
         sd = sigma,
         x = x,
         z = normalize(x))
}

zjohnson <- function(x) {
    # default parameters
    gdxl <- c(gamma = 0, delta = 1, xi = 0, lambda = 1)

    # parameterization
    opt <- optim(gdxl, function(gdxl) {
        z <- qnorm(pjohnson(x, gdxl[1], gdxl[2], gdxl[3], gdxl[4]))
        suppressWarnings(ks.test(z, "pnorm")$p.value)
    })

    # binding
    gdxl <- opt$par
    zn <- znorm(qnorm(pjohnson(x, gdxl[1], gdxl[2], gdxl[3], gdxl[4])))

    # functions
    normalize <- function(q)
        zn$normalize(qnorm(pjohnson(q, gdxl[1], gdxl[2], gdxl[3], gdxl[4])))
    denormalize <- function(z)
        qjohnson(pnorm(zn$denormalize(z)), gdxl[1], gdxl[2], gdxl[3], gdxl[4])

    # return
    list(normalize = normalize,
         denormalize = denormalize,
         gamma = gdxl[1],
         delta = gdxl[2],
         xi = gdxl[3],
         lambda = gdxl[4],
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

    response <- mpg
    predictors <- list(hp, wt, disp, drat, qsec)

    prm <- rep(c(0, 1), 2 * (1 + length(predictors)))
    last <- 1
    better <- TRUE
    while (better) {
        opt <- optim(prm, function(prm) {
            y <- qnorm(pjohnson(response, prm[1], prm[2], prm[3], prm[4]))
            x <- rowSums(sapply(1:length(predictors), function(k) {
                qnorm(pjohnson(predictors[[k]], prm[4 * k + 1], prm[4 * k + 2], prm[4 * k + 3], prm[4 * k + 4]))
            }))
            tryCatch({
                m <- lm(y ~ x)
                summary(m)$coefficients[2, 4]
            }, error = function(e) {
                Inf
            })
        })

        prm <- opt$par
        print(opt$value / last)
        if (0.999 < (opt$value / last)) {
            better <- FALSE
        }
        last <- opt$value
    }
}

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
