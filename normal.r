evaluate_point <- function(i, f, a, w) {
    f(a + i * w)
}

integrate_trapezoidal <- function(f, a, b, n) {
    w <- (b - a) / n
    u <- (f(a) + f(b)) / 2
    v <- sum(sapply(1:(n - 1), evaluate_point, f, a, w))
    w * (u + v)
}

erf <- function(x) {
    (2 / sqrt(pi)) * integrate(function(t) exp(-t ^ 2), 0, x)$value
}

invert <- function(f) {
    function(y) {
        suppressWarnings(optim(0, function(x) {
            abs(y - f(x))
        })$par)
    }
}

dnormie <- function(x, mean = 0, sd = 1) {
    (1 / (sqrt(2 * pi * sd ^ 2))) * exp(-((x - mean) ^ 2) / (2 * sd ^ 2))
}

pnormie <- function(q, mean = 0, sd = 1) {
    (1 + erf((q - mean) / (sd * sqrt(2)))) / 2
}

qnormie <- function(p, mean = 0, sd = 1) {
    mean + sd * sqrt(2) * invert(erf)(2 * p - 1)
}

rnormie <- function(n, mean = 0, sd = 1) {
    (1 + erf((q - mean) / (sd * sqrt(2)))) / 2
}

dnorm(-2)
dnormie(-2)
pnorm(-2)
pnormie(-2)
qnorm(0.05)
qnormie(0.05)
rnorm(1)
rnorm(1)
dnorm

{
    pdf("nlj.pdf", 9, 6)

    i <- seq(-3, 3, length.out = 100)
    x <- sapply(i, erf)
    plot(i, x, type = "l")

    dev.off()
}
