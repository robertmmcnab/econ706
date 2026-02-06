---
title: "Matrix Operations"
format:
  html:
    toc: true
    toc-depth: 3
    toc-location: right
    self-contained: false     # must be false when using webr
urlcolor: blue
filters:
  - webr
execute:
  webR: true
---


::: {.cell}

:::


## Matrix Operations

This module extends our discussion of matrix operations. We begin with simulating least squares estimation before discussing matrix properties. We examine derivative properties before covering the projection and residual maker matrices. Understanding these operations is critical as they form the backbone of the Ordinary Least Squares (OLS) estimator, which we will derive and implement both manually and using built-in R functions.

You will find these examples useful as you derive the least squares estimator and examine the properties of the least squares estimator in later modules.

## Simulating Data

To illustrate how to use the least squares estimator in R, we can use "simulated" data. The use of simulated data allows us to specify the properties of the observations that are used in the estimation of the least squares coefficients.

For replication, we use the **set.seed** function.

By specifying **set.seed(1234)**, we will obtain the same results each time we run the code.

We first need to create the matrix of independent variables.

Assume that we want to 6,000 observations with 5 independent variables. We will add an additional column for the intercept, so our matrix will be $(6000 \times 6)$.

Our first step is to create the $(6000 \times 6)$ matrix. The matrix is filled with zeros.

We then fill the first column of the **X** matrix with 1's for the intercept term.

The next step is to fill columns 2 thorough 6 of the **X** matrix. We use a loop to fill these columns.

Each observation is filled with a random value from a standard normal distribution.


::: {.cell}

```{.r .cell-code}
# Clear memory
rm(list = ls())

# Set Seed to 1234
set.seed(1234)

# Create a (6000 x 6) matrix of 0s 
x <- matrix(0,6000,6)

# Fill the first column with 1s

x[,1] <- 1

# Loop over remaining columns
# Fill with values from standard normal
for (l in 2:6) {
    x[,l] <- rnorm(6000)
}

# Output the first three rows of the X matrix
kable(x[1:3,],
      digits = 3,
      col.names = c('X1','X2','X3','X4','X5','X6')) %>%
kable_styling(font_size = 12,
              position  = "center",
              full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> X1 </th>
   <th style="text-align:right;"> X2 </th>
   <th style="text-align:right;"> X3 </th>
   <th style="text-align:right;"> X4 </th>
   <th style="text-align:right;"> X5 </th>
   <th style="text-align:right;"> X6 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -1.207 </td>
   <td style="text-align:right;"> 0.023 </td>
   <td style="text-align:right;"> 0.249 </td>
   <td style="text-align:right;"> 1.500 </td>
   <td style="text-align:right;"> -0.069 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.277 </td>
   <td style="text-align:right;"> -0.205 </td>
   <td style="text-align:right;"> -1.030 </td>
   <td style="text-align:right;"> -0.519 </td>
   <td style="text-align:right;"> -0.238 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1.084 </td>
   <td style="text-align:right;"> -1.049 </td>
   <td style="text-align:right;"> -0.341 </td>
   <td style="text-align:right;"> -2.724 </td>
   <td style="text-align:right;"> -1.111 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Estimating Coefficients

We can now specify the population coefficients, that is, as we are simulating population data, we can define the values of the population parameters or:

$$
\beta = \begin{bmatrix}{1 \\ 2 \\ 3 \\ 4 \\ 5 \\ 6} \end{bmatrix}
$$ Let the disturbance term, $\epsilon$ be distributed normally with 0 mean and constant variance equal to 1 or:

$$
\epsilon \sim N(0,1)
$$ We specify that $y$ is a linear equation equal to:

$$
y = X \beta + \epsilon
$$

The least squares estimator is defined as:

$$
b = (X^T X)^{-1} X^T y
$$

We estimate the least squares coefficients and produce a table of coefficients.

You should note that even though we have population data and know the population parameters, the least squares estimates are not equal to the population parameters.

The stochastic disturbance captures unobserved variation in the dependent variable, that is, we can not fully explain the variation in the dependent variable.


::: {.cell}

```{.r .cell-code}
set.seed(1234)
 
x <- matrix(0,6000,6)
      x[,1] <- 1
      for (l in 2:6) {
        x[,l] <- rnorm(6000)
      }

# Population coefficients

beta <- c(1,2,3,4,5,6)

# Disturbance term
epsilon <- rnorm(6000, 0, 1)

# y is equal to X*beta + epsilon

y <- x %*% beta + epsilon

# Solve for b =  (X^T X)^(-1) X^T y

bhat <- solve(t(x) %*% x) %*% t(x) %*% y

b_names <- c("b1", "b2", "b3", "b4", "b5", "b6")

bhat_named <- cbind(b_names, 
                    beta,
                    round(bhat,3))

kable(bhat_named,
      align     = 'c',
      col.names = c("Coefficient",
                    "Population Parameter",
                    "Estimate")) %>%
kable_styling(font_size = 12,
              position  = "center",
              full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:center;"> Coefficient </th>
   <th style="text-align:center;"> Population Parameter </th>
   <th style="text-align:center;"> Estimate </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> b1 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 1.015 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> b2 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 1.995 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> b3 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 3.003 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> b4 </td>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 3.986 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> b5 </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 4.988 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> b6 </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 5.999 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Estimating Using LM Function

We can also compare the least squares estimates using matrix algebra with those obtained using the **lm** function (which uses matrix algebra behind the scenes).

In the example below, we replicate the methodology above to obtain the least squares estimates using matrix algebra. Note that we have changed the number of observations to 10,000.

To use the **lm** function and obtain results comparable to those using matrix algebra, we only need to include the independent variables. This means we must exclude the intercept variable (**X1**) from the data that will be used by the **lm** function.

Why?

The **lm** function, by default, estimates an intercept. If we included the intercept variable (which is a column of 1's in the **X** matrix), then our results would not be comparable.

We find that the estimates are the same, regardless if we use matrix algebra or the **lm** function.


::: {.cell}

```{.r .cell-code}
rm(list = ls())

set.seed(1234)

x <- matrix(0,10000,6)
      x[,1] <- 1
      for (l in 2:6) {
        x[,l] <- rnorm(10000)
        gc()
      }

beta <- c(7,8,9,10,11,12)

epsilon <- rnorm(10000, 0, 1)

# y is equal to X*beta + random disturbance

y <- x %*% beta + epsilon

# Solve for b =  (X^T X)^(-1) X^T y

bhat <- solve(t(x) %*% x) %*% t(x) %*% y

b_names <- c("b1", "b2", "b3", "b4", "b5", "b6")

# Create dataframe for LM function
# Include Y and independent variables 

data <- tibble(y = y[,1],
               x2 = x[,2],
               x3 = x[,3],
               x4 = x[,4],
               x5 = x[,5],
               x6 = x[,6])

# Generate least squares estimates using LM function

lm_mod <- lm(y ~ x2 + x3 + x4 + x5 + x6, data = data)$coefficients

# Bind population parameters and estimates together
bhat_named <- cbind(b_names, 
                     beta, 
                     round(bhat,3),
                     round(lm_mod,3))
kable(bhat_named,
      align     = 'c',
      col.names = c("Coefficient",
                    "Population Parameter",
                    "LS Estimates",
                    "LM Estimates")) %>%
kable_styling(font_size = 12,
              position  = "center",
              full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">  </th>
   <th style="text-align:center;"> Coefficient </th>
   <th style="text-align:center;"> Population Parameter </th>
   <th style="text-align:center;"> LS Estimates </th>
   <th style="text-align:center;"> LM Estimates </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> (Intercept) </td>
   <td style="text-align:center;"> b1 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 6.992 </td>
   <td style="text-align:center;"> 6.992 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x2 </td>
   <td style="text-align:center;"> b2 </td>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 7.991 </td>
   <td style="text-align:center;"> 7.991 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x3 </td>
   <td style="text-align:center;"> b3 </td>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 8.995 </td>
   <td style="text-align:center;"> 8.995 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x4 </td>
   <td style="text-align:center;"> b4 </td>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 10.032 </td>
   <td style="text-align:center;"> 10.032 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x5 </td>
   <td style="text-align:center;"> b5 </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 11.001 </td>
   <td style="text-align:center;"> 11.001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> x6 </td>
   <td style="text-align:center;"> b6 </td>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 12.025 </td>
   <td style="text-align:center;"> 12.025 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Matrix Properties

In this section, we list some common matrix properties which are useful in econometrics.

For much more information, see "The Matrix Cookbook" at <http://matrixcookbook.com>

Vectors are lowercase letters while capital letters represent matrices.

The transpose operation reverses the order of an operation so that the transpose of a product of two matrices can be expressed as the product of the transposed matrices in reverse order.

$$(AB)^T = B^T A^T$$

Let $a$ and $c$ be vectors.

The transpose operator reverses the order of the operation.

The transpose of a transposed vector is the vector.\
$$(a^T B c)^T = c^T B^T a$$

$(A + B)$ pre-multiplies $C$ is equal to the distributed addition of $AC$ and $BC$.

$$(A+B)C = AC + BC$$

The transpose of the addition of the vectors $a$ and $b$ pre-multiplies $C$ and the multiplication can be distributed.

$$(a+b)^T C = a^T C + b^T C$$

As noted previously, in almost every case, matrix multiplication is not communicative.

$$AB \ne BA$$

In the code below, we demonstrate each of these properties.

```{webr-r}

rm(list = ls())

a <- c(1, 2, 3, 4)
b <- c(6, 7, 8, 9)

mat_A <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
                nrow = 4,
                ncol = 4)

mat_B <- matrix(c(2,4,6,8,10,12,14,16,18,20,22,24,26,27,30,32),
                nrow = 4,
                ncol = 4)

mat_C <- matrix(c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31),
                nrow = 4,
                ncol = 4)

# (AB)^T = B^T A^T
lhs <- t(mat_A %*% mat_B)
rhs <- t(mat_B) %*% t(mat_A)

all.equal(lhs, rhs)
  
# (a^T B c)^T = c^t B^T a
lhs <- t(t(a) %*% mat_B %*% b)
rhs <- t(b) %*% t(mat_B) %*% a

all.equal(lhs, rhs)

#(A + B)C = AC + BC

lhs <- (mat_A + mat_B) %*% mat_C
rhs <- (mat_A %*% mat_C) + (mat_B %*% mat_C)

all.equal(lhs, rhs)

#(a+b)^T C = a^T C + b^T C

lhs <- t(a + b) %*% mat_C
rhs <- (t(a) %*% mat_C) + (t(b) %*% mat_C)

all.equal(lhs, rhs)

# AB is not equal to BA

lhs <- mat_A %*% mat_B
rhs <- mat_B %*% mat_A

all.equal(lhs, rhs)
```

## Common Derivative Properties

The following rules are often quite useful

First, let's state derivative properties with respect to scalars.

The derivative of $bx$ with respect to $x$ is $b$.

$$\frac{\partial bx}{\partial x} = b $$

The derivative of $x^2$ with respect to $x$ is $2x$.

$$\frac{\partial x^2}{\partial x} = 2x$$

The derivative of $bx^2$ with respect to $x$ is $2bx$.

$$\frac{\partial b x^2}{\partial x} = 2bx$$

We can take these scalar rules and apply them in a similar fashion of vectors and matrices.

$$\frac{\partial x^T B}{\partial x} = B$$

$$\frac{\partial x^T b}{\partial x} = b$$

$$\frac{\partial x^T x}{\partial x} = 2x$$ $$\frac{\partial x^T B x}{\partial x} = (B + B^T)x = 2Bx $$

## Residual Maker Matrix

$M$ is a $(n \times n)$ symmetric, idempotent matrix where

$$M = I - X(X^T X)^{-1}X^T$$

$M$ is the **residual maker** matrix.

When $y$ is pre-multiplied by $M$, it produces the vector of least squares residuals, $e$.

$$
My = (I - X(X^TX)^{-1} X^T = y - X(X^TX)^{-1}X^Ty = y - Xb = e
$$

We also note that $MX = 0$ because the regression of $X$ on $X$ is a 'perfect fit' or:

$$
MX = (I - X(X^TX)^{-1} X^T)X = X - X(X^TX)^{-1}X^TX = X - X = 0
$$


::: {.cell}

```{.r .cell-code}
rm(list = ls())

nrows <- 7000
ncols <- 9
 
x <- matrix(0, nrows, ncols)
      x[,1] <- 1
      for (l in 2: ncols) {
        x[,l] <- rnorm(nrows)
        gc()
      }

beta <- c(1,2,3,4,5,6,7,8,9)

epsilon <- rnorm(nrows, 0, 1)

y <- x %*% beta + epsilon

bhat <- solve(t(x) %*% x) %*% t(x) %*% y

yhat <- x %*% bhat

# Solve for residuals

e <- y - yhat

# Create Identity Matrix

I <- diag(nrows)

# Create Residual Maker

M <- I - x %*% solve(t(x) %*% x) %*% t(x)

# Residuals = My

e1 <- M %*% y

df_e <- data.frame(e, e1)

kable(df_e[1:10,],
      col.names = c("e (y - yhat)",
                    "e (My)")) %>%
kable_styling(font_size = 12,
              position  = "center",
              full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> e (y - yhat) </th>
   <th style="text-align:right;"> e (My) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.7044752 </td>
   <td style="text-align:right;"> 0.7044752 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.5357223 </td>
   <td style="text-align:right;"> 0.5357223 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -1.3074348 </td>
   <td style="text-align:right;"> -1.3074348 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.6846867 </td>
   <td style="text-align:right;"> -0.6846867 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.8841051 </td>
   <td style="text-align:right;"> 0.8841051 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 1.8214250 </td>
   <td style="text-align:right;"> 1.8214250 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 1.2356265 </td>
   <td style="text-align:right;"> 1.2356265 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -1.0181093 </td>
   <td style="text-align:right;"> -1.0181093 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.3338123 </td>
   <td style="text-align:right;"> -0.3338123 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.2873728 </td>
   <td style="text-align:right;"> 0.2873728 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## M is Symmetric

Recall that $M$ is symmetric, that is,

$$M = M^T$$


::: {.cell}

```{.r .cell-code}
rm(list = ls())

nrows <- 3000
ncols <- 5
 
x <- matrix(0, nrows, ncols)
      x[,1] <- 1
      for (l in 2: ncols) {
        x[,l] <- rnorm(nrows)
      }

beta <- c(1,2,3,4,5)

epsilon <- rnorm(nrows, 0, 1)

y <- x %*% beta + epsilon

# Create Identity Matrix

I <- diag(nrows)

# Create Residual Maker

M <- I - x %*% solve(t(x) %*% x) %*% t(x)

# Transpose of M

t_M <- t(M)

kable(data.frame(M)[1:5,1:5]) %>%
kable_styling(font_size = 12,
              position  = "center",
              full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> X1 </th>
   <th style="text-align:right;"> X2 </th>
   <th style="text-align:right;"> X3 </th>
   <th style="text-align:right;"> X4 </th>
   <th style="text-align:right;"> X5 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.9966323 </td>
   <td style="text-align:right;"> 0.0001332 </td>
   <td style="text-align:right;"> 0.0000138 </td>
   <td style="text-align:right;"> 0.0006391 </td>
   <td style="text-align:right;"> -0.0000068 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.0001332 </td>
   <td style="text-align:right;"> 0.9990577 </td>
   <td style="text-align:right;"> 0.0000818 </td>
   <td style="text-align:right;"> -0.0008291 </td>
   <td style="text-align:right;"> -0.0006630 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.0000138 </td>
   <td style="text-align:right;"> 0.0000818 </td>
   <td style="text-align:right;"> 0.9987474 </td>
   <td style="text-align:right;"> 0.0002162 </td>
   <td style="text-align:right;"> 0.0000077 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.0006391 </td>
   <td style="text-align:right;"> -0.0008291 </td>
   <td style="text-align:right;"> 0.0002162 </td>
   <td style="text-align:right;"> 0.9988237 </td>
   <td style="text-align:right;"> -0.0007385 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.0000068 </td>
   <td style="text-align:right;"> -0.0006630 </td>
   <td style="text-align:right;"> 0.0000077 </td>
   <td style="text-align:right;"> -0.0007385 </td>
   <td style="text-align:right;"> 0.9994385 </td>
  </tr>
</tbody>
</table>

`````
:::

```{.r .cell-code}
kable(data.frame(t_M)[1:5,1:5]) %>%
kable_styling(font_size = 12,
              position  = "center",
              full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> X1 </th>
   <th style="text-align:right;"> X2 </th>
   <th style="text-align:right;"> X3 </th>
   <th style="text-align:right;"> X4 </th>
   <th style="text-align:right;"> X5 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.9966323 </td>
   <td style="text-align:right;"> 0.0001332 </td>
   <td style="text-align:right;"> 0.0000138 </td>
   <td style="text-align:right;"> 0.0006391 </td>
   <td style="text-align:right;"> -0.0000068 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.0001332 </td>
   <td style="text-align:right;"> 0.9990577 </td>
   <td style="text-align:right;"> 0.0000818 </td>
   <td style="text-align:right;"> -0.0008291 </td>
   <td style="text-align:right;"> -0.0006630 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.0000138 </td>
   <td style="text-align:right;"> 0.0000818 </td>
   <td style="text-align:right;"> 0.9987474 </td>
   <td style="text-align:right;"> 0.0002162 </td>
   <td style="text-align:right;"> 0.0000077 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.0006391 </td>
   <td style="text-align:right;"> -0.0008291 </td>
   <td style="text-align:right;"> 0.0002162 </td>
   <td style="text-align:right;"> 0.9988237 </td>
   <td style="text-align:right;"> -0.0007385 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.0000068 </td>
   <td style="text-align:right;"> -0.0006630 </td>
   <td style="text-align:right;"> 0.0000077 </td>
   <td style="text-align:right;"> -0.0007385 </td>
   <td style="text-align:right;"> 0.9994385 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## M is Idempotent

$M$ is also idempotent, that is, it is a matrix that, when multiplied by itself, yields itself.

$$
M M = M
$$


::: {.cell}

```{.r .cell-code}
rm(list = ls())

nrows <- 1000
ncols <- 5
 
x <- matrix(0, nrows, ncols)
      x[,1] <- 1
      for (l in 2: ncols) {
        x[,l] <- rnorm(nrows)
      }

beta <- c(1,2,3,4,5)

epsilon <- rnorm(nrows, 0, 1)

y <- x %*% beta + epsilon

# Create Identity Matrix

I <- diag(nrows)

# Create Residual Maker

M <- I - x %*% solve(t(x) %*% x) %*% t(x)


MM <- M %*% M

kable(data.frame(M)[1:5,1:5]) %>%
kable_styling(font_size = 12,
              position  = "center",
              full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> X1 </th>
   <th style="text-align:right;"> X2 </th>
   <th style="text-align:right;"> X3 </th>
   <th style="text-align:right;"> X4 </th>
   <th style="text-align:right;"> X5 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.9937489 </td>
   <td style="text-align:right;"> -0.0013963 </td>
   <td style="text-align:right;"> -0.0018414 </td>
   <td style="text-align:right;"> -0.0023038 </td>
   <td style="text-align:right;"> -0.0021625 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.0013963 </td>
   <td style="text-align:right;"> 0.9958213 </td>
   <td style="text-align:right;"> -0.0004577 </td>
   <td style="text-align:right;"> 0.0008719 </td>
   <td style="text-align:right;"> -0.0026584 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.0018414 </td>
   <td style="text-align:right;"> -0.0004577 </td>
   <td style="text-align:right;"> 0.9982273 </td>
   <td style="text-align:right;"> -0.0014619 </td>
   <td style="text-align:right;"> -0.0007969 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.0023038 </td>
   <td style="text-align:right;"> 0.0008719 </td>
   <td style="text-align:right;"> -0.0014619 </td>
   <td style="text-align:right;"> 0.9973824 </td>
   <td style="text-align:right;"> -0.0003145 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.0021625 </td>
   <td style="text-align:right;"> -0.0026584 </td>
   <td style="text-align:right;"> -0.0007969 </td>
   <td style="text-align:right;"> -0.0003145 </td>
   <td style="text-align:right;"> 0.9979395 </td>
  </tr>
</tbody>
</table>

`````
:::

```{.r .cell-code}
kable(data.frame(MM)[1:5,1:5]) %>%
kable_styling(font_size = 12,
              position  = "center",
              full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> X1 </th>
   <th style="text-align:right;"> X2 </th>
   <th style="text-align:right;"> X3 </th>
   <th style="text-align:right;"> X4 </th>
   <th style="text-align:right;"> X5 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.9937489 </td>
   <td style="text-align:right;"> -0.0013963 </td>
   <td style="text-align:right;"> -0.0018414 </td>
   <td style="text-align:right;"> -0.0023038 </td>
   <td style="text-align:right;"> -0.0021625 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.0013963 </td>
   <td style="text-align:right;"> 0.9958213 </td>
   <td style="text-align:right;"> -0.0004577 </td>
   <td style="text-align:right;"> 0.0008719 </td>
   <td style="text-align:right;"> -0.0026584 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.0018414 </td>
   <td style="text-align:right;"> -0.0004577 </td>
   <td style="text-align:right;"> 0.9982273 </td>
   <td style="text-align:right;"> -0.0014619 </td>
   <td style="text-align:right;"> -0.0007969 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.0023038 </td>
   <td style="text-align:right;"> 0.0008719 </td>
   <td style="text-align:right;"> -0.0014619 </td>
   <td style="text-align:right;"> 0.9973824 </td>
   <td style="text-align:right;"> -0.0003145 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.0021625 </td>
   <td style="text-align:right;"> -0.0026584 </td>
   <td style="text-align:right;"> -0.0007969 </td>
   <td style="text-align:right;"> -0.0003145 </td>
   <td style="text-align:right;"> 0.9979395 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## MX is Zero

$MX$ is equal to zero.


::: {.cell}

```{.r .cell-code}
rm(list = ls())

nrows <- 1000
ncols <- 5
 
x <- matrix(0, nrows, ncols)
      x[,1] <- 1
      for (l in 2: ncols) {
        x[,l] <- rnorm(nrows)
      }

beta <- c(1,2,3,4,5)

epsilon <- rnorm(nrows, 0, 1)

y <- x %*% beta + epsilon

# Create Identity Matrix

I <- diag(nrows)

# Create Residual Maker

M <- I - x %*% solve(t(x) %*% x) %*% t(x)

MX <- M %*% x

kable(data.frame(MX)[1:5,1:5]) %>%
kable_styling(font_size = 12,
              position  = "center",
              full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> X1 </th>
   <th style="text-align:right;"> X2 </th>
   <th style="text-align:right;"> X3 </th>
   <th style="text-align:right;"> X4 </th>
   <th style="text-align:right;"> X5 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Projection Matrix

While $My$ captures the **unexplained variation** in $y$, $Py$ captured the explained variation in $y$ or:

$$M = I - X(X^T X)^{-1}X^T$$

$$P = X(X^T X)^{-1}X^T$$

The projection matrix $P$ is symmetric and idempotent.

$PM$ and $MP$ are equivalent and equal to zero and $P$ and $M$ are orthogonal.

Given the orthogonality, we can **partition** y into two parts:

$$y = Py + My$$ $$y = X(X^T X)^{-1} X^T y + (I-X(X^T X)^{-1}X^T)y$$

Recall that $b = (X^T X)^{-1} X^T y$ and substitute:

$$y = Xb + (I-X(X^T X)^{-1}X^T)y$$

Now multiply $y$ through the second term:

$$y = X b + (y - X(X^T X)^{-1} X^T y)$$

Substitute again for $b$

$$y = X b + (y - X b)$$

Recognize that $Xb = \hat{y}$

$$y = Xb + (y - \hat{y})$$ $$y = Xb + e$$


::: {.cell}

```{.r .cell-code}
rm(list = ls())

nrows <- 1000
ncols <- 5
 
x <- matrix(0, nrows, ncols)
      x[,1] <- 1
      for (l in 2: ncols) {
        x[,l] <- rnorm(nrows)
      }

beta <- c(1,2,3,4,5)

epsilon <- rnorm(nrows, 0, 1)

y <- x %*% beta + epsilon

bhat <- solve(t(x) %*% x) %*% t(x) %*% y

yhat <- x %*% bhat

# Create Identity Matrix

I <- diag(nrows)

M <- I - x %*% solve(t(x) %*% x) %*% t(x)

P <- x %*% solve(t(x) %*% x) %*% t(x)

Py <- P %*% y

df_yhat <- data.frame(yhat, Py)

kable(df_yhat[1:10,]) %>%
kable_styling(font_size = 12,
              position  = "center",
              full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> yhat </th>
   <th style="text-align:right;"> Py </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 3.2036752 </td>
   <td style="text-align:right;"> 3.2036752 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -12.0361510 </td>
   <td style="text-align:right;"> -12.0361510 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 1.5782350 </td>
   <td style="text-align:right;"> 1.5782350 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -4.1662385 </td>
   <td style="text-align:right;"> -4.1662385 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 8.6660633 </td>
   <td style="text-align:right;"> 8.6660633 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -2.1167407 </td>
   <td style="text-align:right;"> -2.1167407 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -0.0707777 </td>
   <td style="text-align:right;"> -0.0707777 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 8.4587151 </td>
   <td style="text-align:right;"> 8.4587151 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 1.8995740 </td>
   <td style="text-align:right;"> 1.8995740 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 5.1218411 </td>
   <td style="text-align:right;"> 5.1218411 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


## Partitioned Y

We can now partition $y$ into explained variation and unexplained variation or:

$$
y = Py + My = Xb + e
$$

Since **P** and **M** are orthogonal, **Py** and **My** are orthogonal.


::: {.cell}

```{.r .cell-code}
rm(list = ls())

nrows <- 1000
ncols <- 5
 
x <- matrix(0, nrows, ncols)
      x[,1] <- 1
      for (l in 2: ncols) {
        x[,l] <- rnorm(nrows)
      }

beta <- c(1,2,3,4,5)

epsilon <- rnorm(nrows, 0, 1)

y <- x %*% beta + epsilon

# Create Identity Matrix

I <- diag(nrows)

M <- I - x %*% solve(t(x) %*% x) %*% t(x)

P <- x %*% solve(t(x) %*% x) %*% t(x)

ypart <- P %*% y + M %*% y

df_y <- data.frame(y, ypart)

kable(df_y[1:10,]) %>%
kable_styling(font_size = 12,
              position  = "center",
              full_width = FALSE)
```

::: {.cell-output-display}
`````{=html}
<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> y </th>
   <th style="text-align:right;"> ypart </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 10.860336 </td>
   <td style="text-align:right;"> 10.860336 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -7.490902 </td>
   <td style="text-align:right;"> -7.490902 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -2.729709 </td>
   <td style="text-align:right;"> -2.729709 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -10.164991 </td>
   <td style="text-align:right;"> -10.164991 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -4.621012 </td>
   <td style="text-align:right;"> -4.621012 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 6.230111 </td>
   <td style="text-align:right;"> 6.230111 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -9.098625 </td>
   <td style="text-align:right;"> -9.098625 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -9.537042 </td>
   <td style="text-align:right;"> -9.537042 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> -11.211872 </td>
   <td style="text-align:right;"> -11.211872 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 14.520627 </td>
   <td style="text-align:right;"> 14.520627 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


