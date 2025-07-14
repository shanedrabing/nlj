
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NLJ

The **nlj** package—named in honor of British statistician Norman Lloyd
Johnson—provides specialized tools for working with the unbounded
Johnson SU distribution. In addition to functions for density
evaluation, distribution and quantile computation, and random variate
generation, the package offers automatic normalization routines and a
generalized transformation-based regression model that leverages the
inverse hyperbolic sine transformation.

------------------------------------------------------------------------

## Installation

Install the development version of **nlj** from GitHub using
**devtools**:

``` r
# install.packages("devtools")
devtools::install_github("shanedrabing/nlj")
```

------------------------------------------------------------------------

## Overview of Functions and Examples

**nlj** provides a suite of tools for the Johnson SU distribution,
automatic data normalization, and advanced regression modeling. The
following examples demonstrate the package’s core functionality.

------------------------------------------------------------------------

### Johnson-SU Distribution Functions

The Johnson-SU distribution functions—**djohnson**, **pjohnson**,
**qjohnson**, and **rjohnson**—provide density evaluation, cumulative
distribution computation, quantile inversion, and random sample
generation. They support the full Johnson-SU parameterization for
flexible control of skewness and kurtosis.

#### Density

**djohnson** computes the Johnson-SU density for a specified set of
input values.

``` r
x <- seq(-pi, pi, length.out = 300)
d <- nlj::djohnson(x, gamma = -1.7, delta = 2.1, xi = -1.45, lambda = 1.45)
plot(x, d, type = "l",
     main = "Johnson-SU Density",
     xlab = "X", ylab = "Density")
```

<img src="man/figures/djohnson-1.png" width="70%" />

#### Cumulative Distribution

**pjohnson** returns the cumulative distribution function (CDF) values
for specified quantiles under the Johnson-SU model.

``` r
q <- seq(-pi, pi, length.out = 300)
p <- nlj::pjohnson(q, gamma = -1.7, delta = 2.1, xi = -1.45, lambda = 1.45)
plot(q, p, type = "l",
     main = "Johnson-SU Cumulative Distribution",
     xlab = "Quantile", ylab = "Probability")
```

<img src="man/figures/pjohnson-1.png" width="70%" />

#### Quantile Function

**qjohnson** obtains the quantile (inverse CDF) corresponding to
specified probability levels, facilitating probability-based analyses
within the Johnson-SU framework.

``` r
p <- seq(0, 1, length.out = 300)
q <- nlj::qjohnson(p, gamma = -1.7, delta = 2.1, xi = -1.45, lambda = 1.45)
plot(p, q, type = "l",
     main = "Johnson-SU Quantile",
     xlab = "Probability", ylab = "Quantile")
```

<img src="man/figures/qjohnson-1.png" width="70%" />

#### Random Deviate Generation

**rjohnson** generates random samples from the Johnson-SU distribution,
useful for simulation, bootstrapping, or Monte Carlo studies.

``` r
nlj::rjohnson(5, gamma = -1.7, delta = 2.1, xi = -1.45, lambda = 1.45)
#> [1]  0.1163819  2.4153221 -0.9112831 -0.1996155  0.4708389
```

------------------------------------------------------------------------

### Automatic Normalization Functions

The package provides **znorm** and **zjohnson**, automated normalization
functions that compute and apply data transformations, track
transformation parameters, and support denormalization. These utilities
streamline the process of preparing data in standardized or
distribution-adapted form.

#### Johnson-SU Normalization

**zjohnson** performs Johnson-SU normalization by estimating
distribution parameters that align the input data to a standard normal
reference. It adapts to skewness and heavy tails, producing a
transformation that accommodates non-Gaussian data characteristics.

``` r
# Example data: Wind speed
x <- sort(airquality$Wind)
d <- density(x)

# Apply Johnson-SU normalization
zj <- nlj::zjohnson(x)

# Extract the density estimates
xd <- d$x
yd <- d$y
yj <- zj$fd(xd)

# Plot the density comparison
plot(range(xd), range(c(yd, yj)), type = "n",
     main = "Wind Speed Density Estimation",
     xlab = "Wind Speed", ylab = "Density")
lines(xd, yd, col = "black")  # Original density
lines(xd, yj, col = "red")   # Johnson-SU normalized density
legend(min(xd), max(c(yd, yj)),
       c("Kernel", "Johnson-SU"),
       col = c("black", "red"),
       lwd = 1, bg = "white")
```

<img src="man/figures/johnson-density-1.png" width="70%" />

------------------------------------------------------------------------

### Generalized Asinh Transformation Model (GATM)

**lm.gat** implements the Generalized Asinh Transformation Model,
fitting regression models with a parameterized inverse hyperbolic sine
transform on predictors. This approach can enhance model fit in the
presence of nonlinearity, particularly when linear models are
insufficient.

#### Automatic Relationship Optimization

The following example illustrates **lm.gat** in action. We predict miles
per gallon (**mpg**) from engine horsepower (**hp**) using both a
standard linear model and a GATM model. The comparison highlights how
the GATM transformation captures non-linear patterns that a simple
linear fit does not.

``` r
# Example data: Horsepower and miles per gallon (MPG)
i <- order(mtcars$hp)
x <- mtcars$hp[i]
y <- mtcars$mpg[i]

# Fit simple regression
m <- lm(y ~ x)

# Fit non-linear model
gat <- nlj::lm.gat(y ~ x, iterations = 3)

# Plot
plot(x, y, pch = 16,
     main = "Horsepower (HP) vs Miles Per Gallon (MPG)",
     xlab = "HP", ylab = "MPG")
lines(x, m$fitted.values, col = "red")
lines(x, gat$z, col = "blue")
legend(min(x), max(y),
       c("Simple", "GATM"),
       col = c("red", "blue"),
       lwd = 1, bg = "white")
```

<img src="man/figures/gatm-1.png" width="70%" />

#### Complex Nonlinear Interaction Optimization

The next example extends **lm.gat** to a model with multiple interaction
terms, demonstrating its ability to manage high-dimensional, non-linear
relationships. We illustrate this by modeling quarter-mile time
(**qsec**) as a function of **hp**, **mpg**, and **wt**, including all
two-way and three-way interactions.

First, we fit a standard linear model with these terms and examine its
summary. Then, we apply **lm.gat** to the same formula, showcasing how
the parameterized asinh transformation can yield a superior fit by
accommodating complex variable interactions.

``` r
# Fit simple regression
m <- lm(qsec ~ hp * mpg * wt, mtcars)

summary(m)
#> 
#> Call:
#> lm(formula = qsec ~ hp * mpg * wt, data = mtcars)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.5219 -0.5882 -0.0978  0.4013  3.1959 
#> 
#> Coefficients:
#>               Estimate Std. Error t value Pr(>|t|)
#> (Intercept) 14.0426771 10.9005197   1.288    0.210
#> hp           0.0017784  0.0590665   0.030    0.976
#> mpg          0.0249053  0.3921269   0.064    0.950
#> wt           0.2427188  3.4981471   0.069    0.945
#> hp:mpg      -0.0005924  0.0026322  -0.225    0.824
#> hp:wt        0.0016601  0.0179826   0.092    0.927
#> mpg:wt       0.1050975  0.1451373   0.724    0.476
#> hp:mpg:wt   -0.0003813  0.0008658  -0.440    0.664
#> 
#> Residual standard error: 1.102 on 24 degrees of freedom
#> Multiple R-squared:  0.7055, Adjusted R-squared:  0.6196 
#> F-statistic: 8.215 on 7 and 24 DF,  p-value: 4.052e-05
```

``` r
# Fit complex regression
gat <- nlj::lm.gat(qsec ~ hp * mpg * wt, mtcars,
                   iterations = 6, penalty = 1e-9)

summary(gat$fit)
#> 
#> Call:
#> stats::lm(formula = formula, data = mut)
#> 
#> Residuals:
#>        Min         1Q     Median         3Q        Max 
#> -2.506e-03 -7.267e-04 -2.282e-05  6.394e-04  2.601e-03 
#> 
#> Coefficients:
#>               Estimate Std. Error  t value Pr(>|t|)    
#> (Intercept)  2.130e+00  1.880e-03 1133.138  < 2e-16 ***
#> hp          -4.703e-04  2.501e-04   -1.881  0.07220 .  
#> mpg         -2.400e-03  2.681e-04   -8.952 4.08e-09 ***
#> wt          -7.029e-03  6.350e-04  -11.069 6.52e-11 ***
#> hp:mpg       1.682e-06  4.890e-05    0.034  0.97284    
#> hp:wt       -3.599e-04  1.176e-04   -3.061  0.00537 ** 
#> mpg:wt      -7.182e-04  9.360e-05   -7.673 6.55e-08 ***
#> hp:mpg:wt   -8.171e-05  1.108e-05   -7.371 1.30e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.001326 on 24 degrees of freedom
#> Multiple R-squared:  0.9224, Adjusted R-squared:  0.8997 
#> F-statistic: 40.74 on 7 and 24 DF,  p-value: 8.352e-12
```
