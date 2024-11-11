
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nlj

<!-- badges: start -->
<!-- badges: end -->

The goal of nlj is to provide functions for statistical transformations,
distributions, and normalization methods, particularly leveraging the
Johnson SU distribution and the generalized inverse hyperbolic sine
transformation.

## Installation

You can install the development version of nlj from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("shanedrabing/nlj")
```

## Example

Below are examples showcasing some key functions in the `nlj` package.

### Total Absolute Deviation from Theoretical CDF

The `qqtad` function calculates the total absolute deviation from a
theoretical cumulative distribution function (CDF):

``` r
library(nlj)

# Example probabilities
p <- c(0.1, 0.3, 0.6, 0.8)

# Calculate the total absolute deviation
qqtad(p)
#> [1] 0.4
```

### Generalized Inverse Hyperbolic Sine Transformation

The `gasinh` and `inverse.gasinh` functions apply and reverse the
generalized inverse hyperbolic sine transformation:

``` r
x <- c(-1, 0, 1, 2, 3)

# Apply generalized inverse hyperbolic sine transformation
transformed <- gasinh(x)
transformed
#> [1] -0.8813736  0.0000000  0.8813736  1.4436355  1.8184465

# Reverse the transformation
inverse_transformed <- inverse.gasinh(transformed)
inverse_transformed
#> [1] -1  0  1  2  3
```

### Johnson SU Distribution Functions

The Johnson SU distribution functions `djohnson`, `pjohnson`,
`qjohnson`, and `rjohnson` provide density, cumulative distribution,
quantile, and random variate generation capabilities.

#### Density

``` r
# Density at x = 0
djohnson(0, gamma = 0, delta = 1, xi = 0, lambda = 1)
#> [1] 0.3989423
```

#### Cumulative Distribution

``` r
# Cumulative probability at q = 0.5
pjohnson(0.5, gamma = 0, delta = 1, xi = 0, lambda = 1)
#> [1] 0.684817
```

#### Quantile Function

``` r
# Quantile for probability p = 0.5
qjohnson(0.5, gamma = 0, delta = 1, xi = 0, lambda = 1)
#> [1] 0
```

#### Random Deviate Generation

``` r
# Generate 5 random deviates
rjohnson(5, gamma = 0, delta = 1, xi = 0, lambda = 1)
#> [1] -0.5492611  0.7361620  0.3126695 -0.4571788  0.7614111
```

### Z-Score Normalization

The `znorm` function normalizes a numeric vector using Z-score
normalization:

``` r
# Sample data
x <- c(1, 2, 3, 4, 5)

# Apply Z-score normalization
z_norm <- znorm(x)
z_norm$z  # Normalized values
#> [1] -1.2649111 -0.6324555  0.0000000  0.6324555  1.2649111

# Reverse the normalization
denormalized <- z_norm$denormalize(z_norm$z)
denormalized
#> [1] 1 2 3 4 5
```

### Johnson SU Normalization

The `zjohnson` function normalizes data using the Johnson SU
distribution:

``` r
# Normalize using Johnson SU distribution
johnson_norm <- zjohnson(x)
johnson_norm$z  # Normalized values
#> [1] -1.214545721 -0.592801987  0.009346948  0.592802011  1.158435547

# Reverse the normalization
denormalized_johnson <- johnson_norm$denormalize(johnson_norm$z)
denormalized_johnson
#> [1] 1 2 3 4 5
```

### Generalized Asinh Transformation Model (GATM)

The `lm.gat` function fits a linear model using the generalized asinh
transformation. Here’s an example of its usage:

``` r
# Generate some example data
set.seed(123)
example_data <- data.frame(
  x = rnorm(100),
  y = 2 * rnorm(100) + 3
)

# Fit the GATM model
model <- lm.gat(y ~ x, data = example_data, iterations = 10, penalty = 1e-6, verbose = TRUE)
#> Loss: 1.341005e-02
#> Loss: 1.319487e-02
#> Loss: 1.31893e-02
#> Loss: 1.31893e-02
#> Loss: 1.31893e-02
#> Loss: 1.31893e-02
#> Loss: 1.31893e-02
#> Loss: 1.31893e-02
#> Loss: 1.31893e-02
#> Loss: 1.31893e-02

# View fitted model parameters
model$par
#> [1] -1.332197e+02  9.233658e+01  1.136046e+00 -1.027282e-29
```
