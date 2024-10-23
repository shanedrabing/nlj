# GENERAL


# generalized inverse hyperbolic sine
gasinh <- function(x, xi = 0, lambda = 1) {
    asinh((x - xi) / lambda)
}

inverse.gasinh <- function(y, xi = 0, lambda = 1) {
    sinh(y) * lambda + xi
}

# total absolute deviation from theoretical CDF
qqtad <- function(p) {
    sum(abs(sort(p) - seq(0, 1, length.out = length(p))))
}


# JOHNSON SU DISTRIBUTION


djohnson <- function(x, gamma = 0, delta = 1, xi = 0, lambda = 1) {
    z <- (x - xi) / lambda
    scale <- delta / (lambda * sqrt(2 * pi))
    shape <- 1 / sqrt(1 + z ^ 2)
    tail <- exp(-0.5 * (gamma + delta * asinh(z)) ^ 2)
    scale * shape * tail
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

    # binding
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
         par = opt$par,
         x = x,
         z = normalize(x))
}


# GENERALIZED ASINH TRANSFORMATION MODEL (GATM)


use.gasinh <- function(i, preds, prm) {
    gasinh(preds[[i]], prm[2 * i + 1], prm[2 * i + 2])
}

capply <- function(X, FUN, ...) {
    do.call(cbind, lapply(X, FUN, ...))
}

lm.johnson <- function(resp, preds, iterations = 3, lambda = 1e-9, verbose = FALSE) {
    # initialization
    prm <- rep(c(0, 1), 1 + length(preds))

    # binding
    evaluate <- function(prm) {
        y <- gasinh(resp, prm[1], prm[2])
        x <- capply(1:length(preds), use.gasinh, preds, prm)

        tryCatch({
            m <- lm.fit(cbind(1, x), y)
            loss <- sum(log(cosh(m$residuals)))
            penalty <- lambda * mean(prm ^ 2)
            loss + penalty
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


# TESTING


resp <- mtcars$qsec
preds <- list(mtcars$mpg, mtcars$hp, mtcars$wt)

raw <- with(mtcars, lm(qsec ~ mpg + hp + wt))
ari <- with(mtcars, lm(log(qsec) ~ log(mpg) + log(hp) + log(wt)))
geo <- with(mtcars, lm(log(qsec) ~ log(mpg) * log(hp) * log(wt)))
lmj <- lm.johnson(resp, preds, iterations = 5, lambda = 1e-15, verbose = TRUE)

round(lmj$par, 3)

sum((mtcars$qsec - raw$fitted.values) ^ 2)
sum((mtcars$qsec - exp(ari$fitted.values)) ^ 2)
sum((mtcars$qsec - exp(geo$fitted.values)) ^ 2)
sum((mtcars$qsec - lmj$z) ^ 2)

df <- cbind(mtcars[, c(1, 4, 6, 7)],
            qsec_raw = raw$fitted.values,
            qsec_ari = exp(ari$fitted.values),
            qsec_geo = exp(geo$fitted.values),
            qsec_lmj = lmj$z)

{
    pdf("nlj.pdf", 7, 6)

    plot(df$qsec_geo / df$qsec, pch = 16, log = "y")
    points(df$qsec_raw / df$qsec, pch = 16, col = "orange")
    points(df$qsec_lmj / df$qsec, pch = 16, col = "green")
    abline(h = 1, col = "red")

    dev.off()
}

round(df, 2)
