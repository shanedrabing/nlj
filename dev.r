# PACKAGING

# 0) Set working directory
setwd("/home/spooks/Code/github/nlj")

# 1) Create framework
# devtools::create("/Users/Shane/Code/packaging/nlj")

# 2) usethis (run once)
# usethis::use_readme_rmd()
# usethis::use_mit_license()
# usethis::use_vignette("nlj_guide")
# usethis::use_testthat()

# 3) Documentation
devtools::document()
devtools::build_readme()
# devtools::build_vignettes()

# 4) check code
devtools::check()
# devtools::check_man()
# devtools::spell_check()
# devtools::run_examples()

# test code
# devtools::unload()
# devtools::load_all()
# devtools::test()

# build
# devtools::build()


# TESTING


{
    # reset    
    tryCatch(devtools::unload(), error = function(e) {})
    devtools::load_all()
}

{
    detach(mtcars)
    attach(mtcars)
}

formula <- qsec ~ hp * mpg * wt

m <- lm(qsec ~ hp * mpg * wt, mtcars)
gat <- nlj::lm.gat(qsec ~ hp * mpg * wt, mtcars,
                   iterations = 6, penalty = 1e-9, verbose = TRUE)

summary(m)
summary(gat$fit)

df <- data.frame(qsec = mtcars$qsec,
                 simple = m$fitted.values,
                 nonlin = gat$z)

round(df, 2)
sum(abs(log(df$qsec / df$simple)))
sum(abs(log(df$qsec / df$nonlin)))


# NUMERICALLY STABLE SUM(LOG(COSH(X)))


sum_log_cosh <- function(m) {
    sum(log(cosh(m$residuals)))
}

sum_log_cosh_stable <- function(m) {
    r <- m$residuals
    sum(r + log1p(exp(-2 * r)) - log(2))
}

x <- NULL
x$residuals <- rnorm(100)
100 * sum_log_cosh(x) / sum_log_cosh_stable(x)


# MORE EFFICIENT QQTAD


qqtad <- function(p) {
    sum(abs(p - seq(0, 1, length.out = length(p))))
}

qqtad_fast <- function(p) {
    n <- length(p)
    ref <- (seq_len(n) - 1) / (n - 1)
    sum(abs(p - ref))
}

n <- 1e4
system.time(replicate(n, qqtad(rnorm(n))))
system.time(replicate(n, qqtad_fast(rnorm(n))))


# BOUNDED OPTIM (ZJOHNSON)


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

zjohnson_bounded <- function(x) {
    # Initial parameters
    gdxl <- c(gamma = 0, delta = 1, xi = 0, lambda = 1)

    # Bounds for constrained optimization
    lower_bd <- c(gamma = -Inf, delta = 1e-6, xi = -Inf, lambda = 1e-6)
    upper_bd <- c(gamma =  Inf,  delta = Inf,  xi =  Inf,  lambda = Inf)

    # Stage 1: Optimize using total absolute deviation from uniform (qqtad)
    opt <- stats::optim(par = gdxl,
                        fn = function(gdxl) {
                            p <- pjohnson(x, gdxl[1], gdxl[2], gdxl[3], gdxl[4])
                            qqtad(p[order(x)])
                        },
                        method = "L-BFGS-B",
                        lower = lower_bd,
                        upper = upper_bd)

    gdxl <- opt$par

    # Stage 2: Refine using KS test on transformed probabilities
    opt <- stats::optim(par = gdxl,
                        fn = function(gdxl) {
                            p <- pjohnson(x, gdxl[1], gdxl[2], gdxl[3], gdxl[4])
                            ks <- suppressWarnings(stats::ks.test(stats::qnorm(p), "pnorm"))
                            -ks$p.value
                        },
                        method = "L-BFGS-B",
                        lower = lower_bd,
                        upper = upper_bd)

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

# Example data: Wind speed
x <- sort(airquality$Wind)
d <- density(x)

# Apply Johnson-SU normalization
zj <- zjohnson(x)
zjb <- zjohnson_bounded(x)

zj$par
zjb$par

# Extract the density estimates
xd <- d$x
yd <- d$y
yj <- zj$fd(xd)
yjb <- zjb$fd(xd)

# Plot the density comparison
{
    plot(range(xd), range(c(yd, yj)), type = "n",
         main = "Wind Speed Density Estimation",
         xlab = "Wind Speed", ylab = "Density")
    lines(xd, yd, col = "black")  # Original density
    lines(xd, yj, col = "red")   # Johnson-SU normalized density
    lines(xd, yjb, col = "blue")   # Johnson-SU normalized density
    legend(15, max(c(yd, yj)),
           c("Kernel", "Johnson-SU", "Johnson-SU (Bounded)"),
           col = c("black", "red", "blue"),
           lwd = 1, bg = "white")
}


# MODEL PERFORMANCE COMPARISON


m <- lm(qsec ~ hp * mpg * wt, mtcars)
summary(m)

gat <- lm.gat(qsec ~ hp * mpg * wt, mtcars, iterations = 6, penalty = 1e-9)
summary(gat$fit)


# COMPARISONS


# USArrests dataset
summary(scale(summary(lm(Murder ~ Assault * UrbanPop * Rape, data = USArrests))$residuals))
summary(scale(summary(lm.gat(Murder ~ Assault * UrbanPop * Rape, data = USArrests, iterations = 6, penalty = 1e-9)$fit)$residuals))

# swiss dataset
summary(scale(summary(lm(Fertility ~ Agriculture * Examination * Education, data = swiss))$residuals))
summary(scale(summary(lm.gat(Fertility ~ Agriculture * Examination * Education, data = swiss, iterations = 6, penalty = 1e-9)$fit)$residuals))

# trees dataset
summary(scale(summary(lm(Volume ~ Girth * Height, data = trees))$residuals))
summary(scale(summary(lm.gat(Volume ~ Girth * Height, data = trees, iterations = 6, penalty = 1e-9)$fit)$residuals))

# airquality dataset (omit NAs)
summary(scale(summary(lm(Ozone ~ Solar.R * Wind * Temp, data = na.omit(airquality)))$residuals))
summary(scale(summary(lm.gat(Ozone ~ Solar.R * Wind * Temp, data = na.omit(airquality), iterations = 6, penalty = 1e-9)$fit)$residuals))

# faithful dataset (simple nonlinear model)
summary(scale(summary(lm(eruptions ~ waiting, data = faithful))$residuals))
summary(scale(summary(lm.gat(eruptions ~ waiting, data = faithful, iterations = 6, penalty = 1e-9)$fit)$residuals))

# CO2 dataset (numeric subset only)
summary(scale(summary(lm(uptake ~ conc, data = subset(CO2, select = c(conc, uptake))))$residuals))
summary(scale(summary(lm.gat(uptake ~ conc, data = subset(CO2, select = c(conc, uptake)), iterations = 6, penalty = 1e-9)$fit)$residuals))
