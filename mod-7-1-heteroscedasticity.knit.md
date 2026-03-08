---
title: "Heteroscedasticity"
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
editor: 
  markdown: 
    wrap: 72
---


::: {.cell}

:::


## Introduction

In this section, we develop the variance estimator for the least squares
estimator. We then discuss the problem of heteroscedastic disturbances
and the impact of heteroscedasticity on the standard errors of the
parameter estimates. We discuss the need for a consistent estimator of
the variance of the least squares estimator.

## Finite Sample Variance

Recall, in finite samples, that
$$var(\boldsymbol{b}) = E[\boldsymbol{(b - \beta)(b - \beta)^T}]$$
$$var(\boldsymbol{b|X}) = E[\boldsymbol{[(X^T X)^{-1} X^T \epsilon]\ [(X^T X)^{-1} X^T \epsilon]^T | X}]$$
$$var(\boldsymbol{b|X}) = E[\boldsymbol{(X^T X)^{-1} X^T \epsilon \epsilon^T X (X^T X)^{-1} | X}]$$
$$var(\boldsymbol{b|X}) = \boldsymbol{(X^T X)^{-1} X^T\ E[\epsilon \epsilon^T|X]\ X (X^T X)^{-1} }$$

What is $E[\boldsymbol{\epsilon \epsilon^T |X}]?$
$$E[\boldsymbol{\epsilon \epsilon^T|X}] = 
  \begin{bmatrix}
  E[\epsilon_1^2|X] & E[\epsilon_1 \epsilon_2|X] & \dots  & E[\epsilon_1 \epsilon_n|X] \\
  E[\epsilon_2 \epsilon_1|X] & E[\epsilon_2^2|X] & \dots  & E[\epsilon_2 \epsilon_n|X] \\
  \vdots & \vdots & \ddots & \vdots \\
  E[\epsilon_n \epsilon_1|X] & E[\epsilon_n \epsilon_1|X] & \dots  & E[\epsilon_n^2|X]
  \end{bmatrix}$$

## Homoscedastic Disturbances

If the assumption of nonautocorrelation holds, then:

$$
E[\epsilon_i, \epsilon_j] = 0 \quad \forall i \ne j
$$

$$E[\epsilon \epsilon^T|X] =
      \begin{bmatrix}
      E[\epsilon_1^2|X] & 0 & \dots  & 0 \\
      0 & E[\epsilon_2^2|X] & \dots  & 0 \\
      \vdots & \vdots & \ddots & \vdots \\
      0 & 0 & \dots  & E[\epsilon_n^2|X]
      \end{bmatrix}$$

If the assumption of constant variance of the disturbances holds, then:

$$
E[\epsilon_i, \epsilon_j | X] = E[\epsilon_i, \epsilon_j] = \sigma^2 \quad \forall i=j
$$

$$E[\epsilon \epsilon^T|X] =
    \begin{bmatrix}
    \sigma^2 & 0 & \dots  & 0 \\
    0 & \sigma^2 & \dots  & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \dots  & \sigma^2
    \end{bmatrix} = \sigma^2 \boldsymbol{I}$$

-   Given the assumptions of nonautocorrelation and homoscedasticity:
    $$var(\boldsymbol{b|X}) = \boldsymbol{(X^T X)^{-1} X^T\ E[\epsilon \epsilon^T|X]\ X (X^T X)^{-1} } = \sigma^2 \boldsymbol{(X^T X)^{-1}}$$

## Visual Representation

In the code chunk below, we generate 10,000 observations over 4
independent variables and 1 constant term. We define the disturbances
as:

$$
\epsilon \sim N(0, 100)
$$ In the resulting graphic, we plot the disturbances with respect to
**X2**. The disturbances appear to be 'white noise', that is, randomly
distributed around the mean of 0 with a variance of 100.

If the disturbances were **heteroscedastic** with respect to **X2**, the
variance of the disturbances would change depending on the values of
**X2**. We'll explore that in coming sections.


::: {.cell}

```{.r .cell-code}
rm(list = ls())
set.seed(1234)

nrows <- 10000
ncols <- 5
 
x <- matrix(0, nrows, ncols)
      x[,1] <- 1
      for (l in 2: ncols) {
        x[,l] <- runif(nrows, min = 0, max = 20)
        gc()
      }

beta <- c(1, 3, 5, 7, 12)

epsilon <- rnorm(nrows, sd = 10)

y <- x %*% beta + epsilon

plotdata <- data.frame(y, x, epsilon)

ggplot(data = plotdata,
       aes(x = X2,
           y = epsilon)) + 
geom_point(color = "blue") +
theme_minimal() +
labs(y = "Disturbances",
     x = "X2",
     title = "Homoscedastic Disturbances")
```

::: {.cell-output-display}
![](mod-7-1-heteroscedasticity_files/figure-html/unnamed-chunk-1-1.png){width=672}
:::
:::


## Heteroscedasticity

First, we continue to assume the assumption of nonautocorrelation
continues to hold.

Second, we now assume the assumption of homoscedasticity is violated,
that is, the disturbances are **heteroscedastic** or:

$$
E[\epsilon_i,\epsilon_j | \boldsymbol{X}] = 0 \quad \forall i \ne j
$$\

$$
E[\epsilon_i,\epsilon_j | \boldsymbol{X}] = \sigma_i^2 \quad \forall i = i
$$

In other words, if the assumption of the constant variances of the
disturbances is violated, the the variance of each disturbance is a
function of one or more of the independent variables or:
$$E[\boldsymbol{\epsilon \epsilon^T|X}] = \boldsymbol{\Omega} 
  =
  \begin{bmatrix}
  \sigma_1^2 & 0 & \dots  & 0 \\
  0 & \sigma_2^2 & \dots  & 0 \\
  \vdots & \vdots & \ddots & \vdots \\
  0 & 0 & \dots  & \sigma_n^2
  \end{bmatrix}\quad \sigma_i^2 \ne \sigma^2,\ \forall i$$

Given heteroscedastic disturbances such that
$var(\boldsymbol{\epsilon}) = \boldsymbol{\Omega}$:

$$var(\boldsymbol{b|X}) = \boldsymbol{(X^T X)^{-1} X^T\ E[\epsilon \epsilon^T|X]\ X (X^T X)^{-1}}$$
$$var(\boldsymbol{b|X}) = \boldsymbol{(X^T X)^{-1} \left( X^T\ \Omega\ X \right) (X^T X)^{-1}}$$

As you can see, the $\boldsymbol{\Omega}$ is the variance-covariance
matrix of the disturbances. Since the diagonal terms (the variances of
the disturbances) are no longer equal to a constant, the
variance-covariance matrix of the disturbances no longer 'collapses' to
$\sigma^2 \boldsymbol{I}$.

In other words, we cannot factor out $\sigma^2 \boldsymbol{I}$ which
means that the variance of the least squares estimator changes as the
values of the independent variables change.

**The least squares estimator remains unbiased and consistent, however,
the standard errors of the least squares estimator are biased, and thus
our inferences are unreliable.**

The question is: how can we test for the presence of heteroscedasticity?

## Visual Representation


::: {.cell}

```{.r .cell-code}
rm(list = ls())
set.seed(1234)

nrows <- 10000
ncols <- 5
 
x <- matrix(0, nrows, ncols)
      x[,1] <- 1
      for (l in 2: ncols) {
        x[,l] <- runif(nrows, min = 0, max = 20)
        gc()
      }

beta <- c(1, 3, 5, 7, 12)

epsilon_het <- rnorm(nrows, sd = x[,2])

y <- x %*% beta + epsilon_het

plotdata <- data.frame(y, x, epsilon_het)

ggplot(data = plotdata,
       aes(x = X2,
           y = epsilon_het)) + 
geom_point(color = "red") +
theme_minimal() +
labs(y = "Disturbances",
     x = "X2",
     title = "Heteroscedastic Disturbances")
```

::: {.cell-output-display}
![](mod-7-1-heteroscedasticity_files/figure-html/unnamed-chunk-2-1.png){width=672}
:::
:::


## Breusch-Pagan Test

A more formal, mathematical way of detecting heteroscedasticity is the
**Breusch-Pagan test**. It involves using a variance function and using
a $\chi^2$ test to test the null hypothesis that heteroskedasticity is
not present (i.e. homoscedastic) against the alternative hypothesis that
heteroscedasticity is present.

To examine the null and alternative hypotheses, we need a **variance
function**. A variance function specifies a relationship between the
variance of the errors and a set of explanatory variables.

Given that $y$ is our dependent variable, $\epsilon$ are our population
disturbances, $\sigma^2_i$ is the variance of $\epsilon_i$ and $Z$ is a
set of exogenous regressors, then we can state:

$$var(\epsilon_i) = \sigma^2_i = E[\epsilon_i^2] = h(Z)$$

We note that $Z$ and $X$ may be the same or may contain different
elements, that is, the set of exogenous regressors for the regression

$$y = X\beta + \epsilon$$

may be different then the set of exogenous regressors for

$$var(\epsilon) = Z\delta + \theta$$ Given that

$$var(\epsilon_i) = \sigma^2_i = E[\epsilon_i^2]$$

we can now specify the model for the $\sigma^2_i$ as:

$$var(\epsilon_i) = E[\epsilon_i^2] = \sigma^2_i = h(z_0 + \delta_1 z_{i1} + \delta_2 z_{i2} + ... + \delta_k z_{ik})$$

We can now observe how $\sigma^2_i$ changes for each observation
depending on the values of the $z_i$ variables.

If $\delta_1 = \delta_2 = .. = \delta_k = 0$, then the $z_i$ have no
influence and $var(\epsilon_i) = \sigma^2_i = h(z_0)$, that is, we have
constant variance and homoscedasticity.

If, on the other hand, $\delta_i \ne 0$ for any , then at least one
$z_i$ influences $var(\epsilon_i) = \sigma^2_i$.

For example, if $z_1$ is influential, then

$$var(\epsilon_i) = E[\epsilon_i^2]= \sigma^2_i = h(z_0 + \delta_1 z_1)$$

We no longer have constant variance and the errors are heteroscedastic.

The above discussion leads to the null and alternative hypotheses of the
Breusch-Pagan test.

The null and alternative hypotheses are:
$$H_0: \delta_1 = \delta_2 = \delta_3 = ... + \delta_k = 0$$
$$H_A: \text{At least one} \  \delta_i \ \text{is not zero}$$

We now specify a linear variance function to obtain:
$$\epsilon^2_i = E[\epsilon_i^2] = z_0 + \delta_1 z_{i1} + \delta_2 z_{i2} + ... + \delta_k z_{ik}$$

Assume that $\upsilon_i = \epsilon^2_i - E[\epsilon_i^2]$, that is, we
can denote $\upsilon$ as difference between a squared error term and its
mean. From the above equation, we can use $\upsilon$ and write
$$\epsilon^2_i = E[\epsilon_i^2] + \upsilon_i = z_0 + \delta_1 z_{i1} + \delta_2 z_{i2} + ... + \delta_k z_{ik} + \upsilon_i$$

Recall that the population errors $\epsilon_i$ are unobservable,
however, we can use the least squares estimates of the population
disturbances, that is, the residuals. We can then estimate the following
model:
$$e^2_i = \hat{z_0} + \hat{\delta_1} z_{i1} + \hat{\delta_2} z_{i2} + ... + \hat{\delta_k} z_{ik} + \hat{\upsilon_i}$$

The intuition for the Breusch-Pagan test is straightforward.

**Do the independent variables explain the variation in the least
squares residuals** $e^2_i$**.**

Since $R^2$ from the above regression is a measure of the variation in
$e_i^2$ that is explained by the regressors, then $R^2$ can be used to
form a test statistic.

Under the null (homoscedasticity), the sample size $N$ times $R^2$ is
distributed $\chi^2$ with $K-1$ degrees of freedom or
$\chi^2 = N \times R^2 \sim \chi^2_{K-1}$.

## Testing for Heteroscedasticity

Using the *CollegeDistance* data, we want to explore whether educational
attainment is a function of proximity to an institution of higher
education or:

$$
education_i = f(distance_i)
$$

We have obtained 4,379 observations which cover students from
approximately 1,100 high schools and is based on survey data conducted
by the Department of Education in 1980, with a follow-up in 1986. We
argue that educational attainment in a function of distance, gender,
whether either parent is a college graduate, income, and region. We
specify a linear regression model as:

$$
education_i = \beta_0 + \beta_1 \, distance_i + \beta_2 \, gender_i + \beta_3 \, father_i + \beta_4 \, mother_i + \\ \beta_5 \, income_i + \beta_6 \, region_i + \epsilon_i
$$

Our testable hypothesis of interest is:

$$
H_o: \beta_1 = 0 \\
H_a: \beta_1 \ne 0
$$

We first estimate the model under the assumption of homoscedasticity. We
present the results in Table 1. To present these results, we use the
**modelsummary** package. You can read more about the package here:
<https://vincentarelbundock.github.io/modelsummary/reference/modelsummary.html>


::: {.cell}

```{.r .cell-code}
rm(list = ls())
library(AER)
library(modelsummary)

data("CollegeDistance")

mod_1 <- lm(education ~ distance + gender + fcollege + mcollege + 
              income + region, data = CollegeDistance)

modelsummary(mod_1,
             stars = TRUE,
             title = "Table 1")
```

::: {.cell-output-display}

```{=html}
<!-- preamble start -->

    <script src="https://cdn.jsdelivr.net/gh/vincentarelbundock/tinytable@main/inst/tinytable.js"></script>

    <script>
      // Create table-specific functions using external factory
      const tableFns_s0c0wp0wgzp2mtekdlcd = TinyTable.createTableFunctions("tinytable_s0c0wp0wgzp2mtekdlcd");
      // tinytable span after
      window.addEventListener('load', function () {
          var cellsToStyle = [
            // tinytable style arrays after
          { positions: [ { i: '22', j: 2 } ], css_id: 'tinytable_css_8vorip105zdm1h5cui6u',}, 
          { positions: [ { i: '14', j: 2 } ], css_id: 'tinytable_css_2g28wlkxel9rteazuwh0',}, 
          { positions: [ { i: '1', j: 2 }, { i: '2', j: 2 }, { i: '3', j: 2 }, { i: '4', j: 2 }, { i: '5', j: 2 }, { i: '6', j: 2 }, { i: '7', j: 2 }, { i: '8', j: 2 }, { i: '9', j: 2 }, { i: '10', j: 2 }, { i: '11', j: 2 }, { i: '12', j: 2 }, { i: '13', j: 2 }, { i: '15', j: 2 }, { i: '16', j: 2 }, { i: '17', j: 2 }, { i: '18', j: 2 }, { i: '19', j: 2 }, { i: '20', j: 2 }, { i: '21', j: 2 } ], css_id: 'tinytable_css_8x41nuhm8x4ojvsbk0uj',}, 
          { positions: [ { i: '0', j: 2 } ], css_id: 'tinytable_css_d44smb0golwdagwxoe6i',}, 
          { positions: [ { i: '22', j: 1 } ], css_id: 'tinytable_css_7h4dgmk8v7fk927cbetr',}, 
          { positions: [ { i: '14', j: 1 } ], css_id: 'tinytable_css_9ssbaixmupv0hlmfzep8',}, 
          { positions: [ { i: '1', j: 1 }, { i: '2', j: 1 }, { i: '3', j: 1 }, { i: '4', j: 1 }, { i: '5', j: 1 }, { i: '6', j: 1 }, { i: '7', j: 1 }, { i: '8', j: 1 }, { i: '9', j: 1 }, { i: '10', j: 1 }, { i: '11', j: 1 }, { i: '12', j: 1 }, { i: '13', j: 1 }, { i: '15', j: 1 }, { i: '16', j: 1 }, { i: '17', j: 1 }, { i: '18', j: 1 }, { i: '19', j: 1 }, { i: '20', j: 1 }, { i: '21', j: 1 } ], css_id: 'tinytable_css_cf57j9a2hf68o73ggcux',}, 
          { positions: [ { i: '0', j: 1 } ], css_id: 'tinytable_css_h1bbpylh99ipzl4fbvjq',}, 
          ];

          // Loop over the arrays to style the cells
          cellsToStyle.forEach(function (group) {
              group.positions.forEach(function (cell) {
                  tableFns_s0c0wp0wgzp2mtekdlcd.styleCell(cell.i, cell.j, group.css_id);
              });
          });
      });
    </script>

    <link rel="stylesheet" href="https://cdn.jsdelivr.net/gh/vincentarelbundock/tinytable@main/inst/tinytable.css">
    <style>
    /* tinytable css entries after */
    #tinytable_s0c0wp0wgzp2mtekdlcd td.tinytable_css_8vorip105zdm1h5cui6u, #tinytable_s0c0wp0wgzp2mtekdlcd th.tinytable_css_8vorip105zdm1h5cui6u {  position: relative; --border-bottom: 1; --border-left: 0; --border-right: 0; --border-top: 0; --line-color-bottom: var(--tt-line-color); --line-color-left: var(--tt-line-color); --line-color-right: var(--tt-line-color); --line-color-top: var(--tt-line-color); --line-width-bottom: 0.1em; --line-width-left: 0.1em; --line-width-right: 0.1em; --line-width-top: 0.1em; --trim-bottom-left: 0%; --trim-bottom-right: 0%; --trim-left-bottom: 0%; --trim-left-top: 0%; --trim-right-bottom: 0%; --trim-right-top: 0%; --trim-top-left: 0%; --trim-top-right: 0%; ; text-align: center }
    #tinytable_s0c0wp0wgzp2mtekdlcd td.tinytable_css_2g28wlkxel9rteazuwh0, #tinytable_s0c0wp0wgzp2mtekdlcd th.tinytable_css_2g28wlkxel9rteazuwh0 {  position: relative; --border-bottom: 1; --border-left: 0; --border-right: 0; --border-top: 0; --line-color-bottom: var(--tt-line-color); --line-color-left: var(--tt-line-color); --line-color-right: var(--tt-line-color); --line-color-top: var(--tt-line-color); --line-width-bottom: 0.05em; --line-width-left: 0.1em; --line-width-right: 0.1em; --line-width-top: 0.1em; --trim-bottom-left: 0%; --trim-bottom-right: 0%; --trim-left-bottom: 0%; --trim-left-top: 0%; --trim-right-bottom: 0%; --trim-right-top: 0%; --trim-top-left: 0%; --trim-top-right: 0%; ; text-align: center }
    #tinytable_s0c0wp0wgzp2mtekdlcd td.tinytable_css_8x41nuhm8x4ojvsbk0uj, #tinytable_s0c0wp0wgzp2mtekdlcd th.tinytable_css_8x41nuhm8x4ojvsbk0uj { text-align: center }
    #tinytable_s0c0wp0wgzp2mtekdlcd td.tinytable_css_d44smb0golwdagwxoe6i, #tinytable_s0c0wp0wgzp2mtekdlcd th.tinytable_css_d44smb0golwdagwxoe6i {  position: relative; --border-bottom: 1; --border-left: 0; --border-right: 0; --border-top: 1; --line-color-bottom: var(--tt-line-color); --line-color-left: var(--tt-line-color); --line-color-right: var(--tt-line-color); --line-color-top: var(--tt-line-color); --line-width-bottom: 0.05em; --line-width-left: 0.1em; --line-width-right: 0.1em; --line-width-top: 0.1em; --trim-bottom-left: 0%; --trim-bottom-right: 0%; --trim-left-bottom: 0%; --trim-left-top: 0%; --trim-right-bottom: 0%; --trim-right-top: 0%; --trim-top-left: 0%; --trim-top-right: 0%; ; text-align: center }
    #tinytable_s0c0wp0wgzp2mtekdlcd td.tinytable_css_7h4dgmk8v7fk927cbetr, #tinytable_s0c0wp0wgzp2mtekdlcd th.tinytable_css_7h4dgmk8v7fk927cbetr {  position: relative; --border-bottom: 1; --border-left: 0; --border-right: 0; --border-top: 0; --line-color-bottom: var(--tt-line-color); --line-color-left: var(--tt-line-color); --line-color-right: var(--tt-line-color); --line-color-top: var(--tt-line-color); --line-width-bottom: 0.1em; --line-width-left: 0.1em; --line-width-right: 0.1em; --line-width-top: 0.1em; --trim-bottom-left: 0%; --trim-bottom-right: 0%; --trim-left-bottom: 0%; --trim-left-top: 0%; --trim-right-bottom: 0%; --trim-right-top: 0%; --trim-top-left: 0%; --trim-top-right: 0%; ; text-align: left }
    #tinytable_s0c0wp0wgzp2mtekdlcd td.tinytable_css_9ssbaixmupv0hlmfzep8, #tinytable_s0c0wp0wgzp2mtekdlcd th.tinytable_css_9ssbaixmupv0hlmfzep8 {  position: relative; --border-bottom: 1; --border-left: 0; --border-right: 0; --border-top: 0; --line-color-bottom: var(--tt-line-color); --line-color-left: var(--tt-line-color); --line-color-right: var(--tt-line-color); --line-color-top: var(--tt-line-color); --line-width-bottom: 0.05em; --line-width-left: 0.1em; --line-width-right: 0.1em; --line-width-top: 0.1em; --trim-bottom-left: 0%; --trim-bottom-right: 0%; --trim-left-bottom: 0%; --trim-left-top: 0%; --trim-right-bottom: 0%; --trim-right-top: 0%; --trim-top-left: 0%; --trim-top-right: 0%; ; text-align: left }
    #tinytable_s0c0wp0wgzp2mtekdlcd td.tinytable_css_cf57j9a2hf68o73ggcux, #tinytable_s0c0wp0wgzp2mtekdlcd th.tinytable_css_cf57j9a2hf68o73ggcux { text-align: left }
    #tinytable_s0c0wp0wgzp2mtekdlcd td.tinytable_css_h1bbpylh99ipzl4fbvjq, #tinytable_s0c0wp0wgzp2mtekdlcd th.tinytable_css_h1bbpylh99ipzl4fbvjq {  position: relative; --border-bottom: 1; --border-left: 0; --border-right: 0; --border-top: 1; --line-color-bottom: var(--tt-line-color); --line-color-left: var(--tt-line-color); --line-color-right: var(--tt-line-color); --line-color-top: var(--tt-line-color); --line-width-bottom: 0.05em; --line-width-left: 0.1em; --line-width-right: 0.1em; --line-width-top: 0.1em; --trim-bottom-left: 0%; --trim-bottom-right: 0%; --trim-left-bottom: 0%; --trim-left-top: 0%; --trim-right-bottom: 0%; --trim-right-top: 0%; --trim-top-left: 0%; --trim-top-right: 0%; ; text-align: left }
    </style>
    <div class="container">
      <table class="tinytable" id="tinytable_s0c0wp0wgzp2mtekdlcd" style="width: auto; margin-left: auto; margin-right: auto;" data-quarto-disable-processing='true'>
        <caption>Table 1</caption>
        <thead>
              <tr>
                <th scope="col" data-row="0" data-col="1"> </th>
                <th scope="col" data-row="0" data-col="2">(1)</th>
              </tr>
        </thead>
        <tfoot><tr><td colspan='2'>+ p < 0.1, * p < 0.05, ** p < 0.01, *** p < 0.001</td></tr></tfoot>
        <tbody>
                <tr>
                  <td data-row="1" data-col="1">(Intercept)</td>
                  <td data-row="1" data-col="2">13.499***</td>
                </tr>
                <tr>
                  <td data-row="2" data-col="1"></td>
                  <td data-row="2" data-col="2">(0.048)</td>
                </tr>
                <tr>
                  <td data-row="3" data-col="1">distance</td>
                  <td data-row="3" data-col="2">-0.041***</td>
                </tr>
                <tr>
                  <td data-row="4" data-col="1"></td>
                  <td data-row="4" data-col="2">(0.011)</td>
                </tr>
                <tr>
                  <td data-row="5" data-col="1">genderfemale</td>
                  <td data-row="5" data-col="2">0.028</td>
                </tr>
                <tr>
                  <td data-row="6" data-col="1"></td>
                  <td data-row="6" data-col="2">(0.049)</td>
                </tr>
                <tr>
                  <td data-row="7" data-col="1">fcollegeyes</td>
                  <td data-row="7" data-col="2">0.837***</td>
                </tr>
                <tr>
                  <td data-row="8" data-col="1"></td>
                  <td data-row="8" data-col="2">(0.070)</td>
                </tr>
                <tr>
                  <td data-row="9" data-col="1">mcollegeyes</td>
                  <td data-row="9" data-col="2">0.568***</td>
                </tr>
                <tr>
                  <td data-row="10" data-col="1"></td>
                  <td data-row="10" data-col="2">(0.079)</td>
                </tr>
                <tr>
                  <td data-row="11" data-col="1">incomehigh</td>
                  <td data-row="11" data-col="2">0.480***</td>
                </tr>
                <tr>
                  <td data-row="12" data-col="1"></td>
                  <td data-row="12" data-col="2">(0.058)</td>
                </tr>
                <tr>
                  <td data-row="13" data-col="1">regionwest</td>
                  <td data-row="13" data-col="2">-0.116+</td>
                </tr>
                <tr>
                  <td data-row="14" data-col="1"></td>
                  <td data-row="14" data-col="2">(0.062)</td>
                </tr>
                <tr>
                  <td data-row="15" data-col="1">Num.Obs.</td>
                  <td data-row="15" data-col="2">4739</td>
                </tr>
                <tr>
                  <td data-row="16" data-col="1">R2</td>
                  <td data-row="16" data-col="2">0.111</td>
                </tr>
                <tr>
                  <td data-row="17" data-col="1">R2 Adj.</td>
                  <td data-row="17" data-col="2">0.110</td>
                </tr>
                <tr>
                  <td data-row="18" data-col="1">AIC</td>
                  <td data-row="18" data-col="2">18421.1</td>
                </tr>
                <tr>
                  <td data-row="19" data-col="1">BIC</td>
                  <td data-row="19" data-col="2">18472.8</td>
                </tr>
                <tr>
                  <td data-row="20" data-col="1">Log.Lik.</td>
                  <td data-row="20" data-col="2">-9202.551</td>
                </tr>
                <tr>
                  <td data-row="21" data-col="1">F</td>
                  <td data-row="21" data-col="2">98.196</td>
                </tr>
                <tr>
                  <td data-row="22" data-col="1">RMSE</td>
                  <td data-row="22" data-col="2">1.69</td>
                </tr>
        </tbody>
      </table>
    </div>
<!-- hack to avoid NA insertion in last line -->
```

:::
:::


However, before we can comment on the sign or magnitude of the estimated
coefficients, we must examine our null hypothesis that the coefficient
for distance is equal to zero. However, we can no longer assume that the
homoscedasticity assumption holds. Instead, we must test the null
hypothesis of homoscedasticity or:

$$
H_o: \sigma_i^2 = \sigma^2 \\
H_a: \sigma_i^2 \ne \sigma^2
$$


::: {.cell}

```{.r .cell-code}
bptest <- bptest(mod_1)

bp <- tibble(bptest$statistic,
             bptest$p.value)

kable(bp,
      col.names = c('BP Test', 'P-Value'),
      align     = 'c',
      digits    = 3) %>%
kable_classic()
```

::: {.cell-output-display}
`````{=html}
<table class=" lightable-classic" style='font-family: "Arial Narrow", "Source Sans Pro", sans-serif; margin-left: auto; margin-right: auto;'>
 <thead>
  <tr>
   <th style="text-align:center;"> BP Test </th>
   <th style="text-align:center;"> P-Value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> 18.252 </td>
   <td style="text-align:center;"> 0.006 </td>
  </tr>
</tbody>
</table>

`````
:::
:::


We thus reject the null hypothesis of homoscedasticity at the 1% level
of statistical significance. We must correct the standard errors before
conducting a hypothesis test.

We can use the **modelsummary** package to report the uncorrected and
corrected standard errors. Note that we have a variety of options to
use. The first column contains the uncorrected standard errors, the 2nd,
and 3rd column contain robust to heteroscedasticity standard errors.
While the standard errors may not have changed much, we must use the 2nd
or 3rd columns to conduct inferences.

We will discuss more about correcting standard errors in the next
section.


::: {.cell}

```{.r .cell-code}
modelsummary(mod_1,
             vcov = list("classical", "HC0", "NeweyWest"),
             stars = TRUE,
             title = "Table 2")
```

::: {.cell-output-display}

```{=html}
<!-- preamble start -->

    <script src="https://cdn.jsdelivr.net/gh/vincentarelbundock/tinytable@main/inst/tinytable.js"></script>

    <script>
      // Create table-specific functions using external factory
      const tableFns_36t5i932ydamvl4zq2c2 = TinyTable.createTableFunctions("tinytable_36t5i932ydamvl4zq2c2");
      // tinytable span after
      window.addEventListener('load', function () {
          var cellsToStyle = [
            // tinytable style arrays after
          { positions: [ { i: '23', j: 2 }, { i: '23', j: 3 }, { i: '23', j: 4 } ], css_id: 'tinytable_css_chk2bf5od93arhnkrgbs',}, 
          { positions: [ { i: '14', j: 2 }, { i: '14', j: 3 }, { i: '14', j: 4 } ], css_id: 'tinytable_css_xozmli066a4omgc90uaz',}, 
          { positions: [ { i: '1', j: 2 }, { i: '2', j: 2 }, { i: '3', j: 2 }, { i: '4', j: 2 }, { i: '5', j: 2 }, { i: '6', j: 2 }, { i: '7', j: 2 }, { i: '8', j: 2 }, { i: '9', j: 2 }, { i: '10', j: 2 }, { i: '11', j: 2 }, { i: '12', j: 2 }, { i: '13', j: 2 }, { i: '15', j: 2 }, { i: '16', j: 2 }, { i: '17', j: 2 }, { i: '18', j: 2 }, { i: '19', j: 2 }, { i: '20', j: 2 }, { i: '21', j: 2 }, { i: '22', j: 2 }, { i: '1', j: 3 }, { i: '2', j: 3 }, { i: '3', j: 3 }, { i: '4', j: 3 }, { i: '5', j: 3 }, { i: '6', j: 3 }, { i: '7', j: 3 }, { i: '8', j: 3 }, { i: '9', j: 3 }, { i: '10', j: 3 }, { i: '11', j: 3 }, { i: '12', j: 3 }, { i: '13', j: 3 }, { i: '15', j: 3 }, { i: '16', j: 3 }, { i: '17', j: 3 }, { i: '18', j: 3 }, { i: '19', j: 3 }, { i: '20', j: 3 }, { i: '21', j: 3 }, { i: '22', j: 3 }, { i: '1', j: 4 }, { i: '2', j: 4 }, { i: '3', j: 4 }, { i: '4', j: 4 }, { i: '5', j: 4 }, { i: '6', j: 4 }, { i: '7', j: 4 }, { i: '8', j: 4 }, { i: '9', j: 4 }, { i: '10', j: 4 }, { i: '11', j: 4 }, { i: '12', j: 4 }, { i: '13', j: 4 }, { i: '15', j: 4 }, { i: '16', j: 4 }, { i: '17', j: 4 }, { i: '18', j: 4 }, { i: '19', j: 4 }, { i: '20', j: 4 }, { i: '21', j: 4 }, { i: '22', j: 4 } ], css_id: 'tinytable_css_x9hpzlcx137u765iietf',}, 
          { positions: [ { i: '0', j: 2 }, { i: '0', j: 3 }, { i: '0', j: 4 } ], css_id: 'tinytable_css_bkby0znnwjx5e4s1vzre',}, 
          { positions: [ { i: '23', j: 1 } ], css_id: 'tinytable_css_ztob0nuw2zjswsr545tn',}, 
          { positions: [ { i: '14', j: 1 } ], css_id: 'tinytable_css_n5xd3q11s98h38sjfhej',}, 
          { positions: [ { i: '1', j: 1 }, { i: '2', j: 1 }, { i: '3', j: 1 }, { i: '4', j: 1 }, { i: '5', j: 1 }, { i: '6', j: 1 }, { i: '7', j: 1 }, { i: '8', j: 1 }, { i: '9', j: 1 }, { i: '10', j: 1 }, { i: '11', j: 1 }, { i: '12', j: 1 }, { i: '13', j: 1 }, { i: '15', j: 1 }, { i: '16', j: 1 }, { i: '17', j: 1 }, { i: '18', j: 1 }, { i: '19', j: 1 }, { i: '20', j: 1 }, { i: '21', j: 1 }, { i: '22', j: 1 } ], css_id: 'tinytable_css_5ypjbsknb8cu51e73d7h',}, 
          { positions: [ { i: '0', j: 1 } ], css_id: 'tinytable_css_iwismptixnys3tql1um4',}, 
          ];

          // Loop over the arrays to style the cells
          cellsToStyle.forEach(function (group) {
              group.positions.forEach(function (cell) {
                  tableFns_36t5i932ydamvl4zq2c2.styleCell(cell.i, cell.j, group.css_id);
              });
          });
      });
    </script>

    <link rel="stylesheet" href="https://cdn.jsdelivr.net/gh/vincentarelbundock/tinytable@main/inst/tinytable.css">
    <style>
    /* tinytable css entries after */
    #tinytable_36t5i932ydamvl4zq2c2 td.tinytable_css_chk2bf5od93arhnkrgbs, #tinytable_36t5i932ydamvl4zq2c2 th.tinytable_css_chk2bf5od93arhnkrgbs {  position: relative; --border-bottom: 1; --border-left: 0; --border-right: 0; --border-top: 0; --line-color-bottom: var(--tt-line-color); --line-color-left: var(--tt-line-color); --line-color-right: var(--tt-line-color); --line-color-top: var(--tt-line-color); --line-width-bottom: 0.1em; --line-width-left: 0.1em; --line-width-right: 0.1em; --line-width-top: 0.1em; --trim-bottom-left: 0%; --trim-bottom-right: 0%; --trim-left-bottom: 0%; --trim-left-top: 0%; --trim-right-bottom: 0%; --trim-right-top: 0%; --trim-top-left: 0%; --trim-top-right: 0%; ; text-align: center }
    #tinytable_36t5i932ydamvl4zq2c2 td.tinytable_css_xozmli066a4omgc90uaz, #tinytable_36t5i932ydamvl4zq2c2 th.tinytable_css_xozmli066a4omgc90uaz {  position: relative; --border-bottom: 1; --border-left: 0; --border-right: 0; --border-top: 0; --line-color-bottom: var(--tt-line-color); --line-color-left: var(--tt-line-color); --line-color-right: var(--tt-line-color); --line-color-top: var(--tt-line-color); --line-width-bottom: 0.05em; --line-width-left: 0.1em; --line-width-right: 0.1em; --line-width-top: 0.1em; --trim-bottom-left: 0%; --trim-bottom-right: 0%; --trim-left-bottom: 0%; --trim-left-top: 0%; --trim-right-bottom: 0%; --trim-right-top: 0%; --trim-top-left: 0%; --trim-top-right: 0%; ; text-align: center }
    #tinytable_36t5i932ydamvl4zq2c2 td.tinytable_css_x9hpzlcx137u765iietf, #tinytable_36t5i932ydamvl4zq2c2 th.tinytable_css_x9hpzlcx137u765iietf { text-align: center }
    #tinytable_36t5i932ydamvl4zq2c2 td.tinytable_css_bkby0znnwjx5e4s1vzre, #tinytable_36t5i932ydamvl4zq2c2 th.tinytable_css_bkby0znnwjx5e4s1vzre {  position: relative; --border-bottom: 1; --border-left: 0; --border-right: 0; --border-top: 1; --line-color-bottom: var(--tt-line-color); --line-color-left: var(--tt-line-color); --line-color-right: var(--tt-line-color); --line-color-top: var(--tt-line-color); --line-width-bottom: 0.05em; --line-width-left: 0.1em; --line-width-right: 0.1em; --line-width-top: 0.1em; --trim-bottom-left: 0%; --trim-bottom-right: 0%; --trim-left-bottom: 0%; --trim-left-top: 0%; --trim-right-bottom: 0%; --trim-right-top: 0%; --trim-top-left: 0%; --trim-top-right: 0%; ; text-align: center }
    #tinytable_36t5i932ydamvl4zq2c2 td.tinytable_css_ztob0nuw2zjswsr545tn, #tinytable_36t5i932ydamvl4zq2c2 th.tinytable_css_ztob0nuw2zjswsr545tn {  position: relative; --border-bottom: 1; --border-left: 0; --border-right: 0; --border-top: 0; --line-color-bottom: var(--tt-line-color); --line-color-left: var(--tt-line-color); --line-color-right: var(--tt-line-color); --line-color-top: var(--tt-line-color); --line-width-bottom: 0.1em; --line-width-left: 0.1em; --line-width-right: 0.1em; --line-width-top: 0.1em; --trim-bottom-left: 0%; --trim-bottom-right: 0%; --trim-left-bottom: 0%; --trim-left-top: 0%; --trim-right-bottom: 0%; --trim-right-top: 0%; --trim-top-left: 0%; --trim-top-right: 0%; ; text-align: left }
    #tinytable_36t5i932ydamvl4zq2c2 td.tinytable_css_n5xd3q11s98h38sjfhej, #tinytable_36t5i932ydamvl4zq2c2 th.tinytable_css_n5xd3q11s98h38sjfhej {  position: relative; --border-bottom: 1; --border-left: 0; --border-right: 0; --border-top: 0; --line-color-bottom: var(--tt-line-color); --line-color-left: var(--tt-line-color); --line-color-right: var(--tt-line-color); --line-color-top: var(--tt-line-color); --line-width-bottom: 0.05em; --line-width-left: 0.1em; --line-width-right: 0.1em; --line-width-top: 0.1em; --trim-bottom-left: 0%; --trim-bottom-right: 0%; --trim-left-bottom: 0%; --trim-left-top: 0%; --trim-right-bottom: 0%; --trim-right-top: 0%; --trim-top-left: 0%; --trim-top-right: 0%; ; text-align: left }
    #tinytable_36t5i932ydamvl4zq2c2 td.tinytable_css_5ypjbsknb8cu51e73d7h, #tinytable_36t5i932ydamvl4zq2c2 th.tinytable_css_5ypjbsknb8cu51e73d7h { text-align: left }
    #tinytable_36t5i932ydamvl4zq2c2 td.tinytable_css_iwismptixnys3tql1um4, #tinytable_36t5i932ydamvl4zq2c2 th.tinytable_css_iwismptixnys3tql1um4 {  position: relative; --border-bottom: 1; --border-left: 0; --border-right: 0; --border-top: 1; --line-color-bottom: var(--tt-line-color); --line-color-left: var(--tt-line-color); --line-color-right: var(--tt-line-color); --line-color-top: var(--tt-line-color); --line-width-bottom: 0.05em; --line-width-left: 0.1em; --line-width-right: 0.1em; --line-width-top: 0.1em; --trim-bottom-left: 0%; --trim-bottom-right: 0%; --trim-left-bottom: 0%; --trim-left-top: 0%; --trim-right-bottom: 0%; --trim-right-top: 0%; --trim-top-left: 0%; --trim-top-right: 0%; ; text-align: left }
    </style>
    <div class="container">
      <table class="tinytable" id="tinytable_36t5i932ydamvl4zq2c2" style="width: auto; margin-left: auto; margin-right: auto;" data-quarto-disable-processing='true'>
        <caption>Table 2</caption>
        <thead>
              <tr>
                <th scope="col" data-row="0" data-col="1"> </th>
                <th scope="col" data-row="0" data-col="2">(1)</th>
                <th scope="col" data-row="0" data-col="3">(2)</th>
                <th scope="col" data-row="0" data-col="4">(3)</th>
              </tr>
        </thead>
        <tfoot><tr><td colspan='4'>+ p < 0.1, * p < 0.05, ** p < 0.01, *** p < 0.001</td></tr></tfoot>
        <tbody>
                <tr>
                  <td data-row="1" data-col="1">(Intercept)</td>
                  <td data-row="1" data-col="2">13.499***</td>
                  <td data-row="1" data-col="3">13.499***</td>
                  <td data-row="1" data-col="4">13.499***</td>
                </tr>
                <tr>
                  <td data-row="2" data-col="1"></td>
                  <td data-row="2" data-col="2">(0.048)</td>
                  <td data-row="2" data-col="3">(0.048)</td>
                  <td data-row="2" data-col="4">(0.062)</td>
                </tr>
                <tr>
                  <td data-row="3" data-col="1">distance</td>
                  <td data-row="3" data-col="2">-0.041***</td>
                  <td data-row="3" data-col="3">-0.041***</td>
                  <td data-row="3" data-col="4">-0.041**</td>
                </tr>
                <tr>
                  <td data-row="4" data-col="1"></td>
                  <td data-row="4" data-col="2">(0.011)</td>
                  <td data-row="4" data-col="3">(0.010)</td>
                  <td data-row="4" data-col="4">(0.014)</td>
                </tr>
                <tr>
                  <td data-row="5" data-col="1">genderfemale</td>
                  <td data-row="5" data-col="2">0.028</td>
                  <td data-row="5" data-col="3">0.028</td>
                  <td data-row="5" data-col="4">0.028</td>
                </tr>
                <tr>
                  <td data-row="6" data-col="1"></td>
                  <td data-row="6" data-col="2">(0.049)</td>
                  <td data-row="6" data-col="3">(0.049)</td>
                  <td data-row="6" data-col="4">(0.053)</td>
                </tr>
                <tr>
                  <td data-row="7" data-col="1">fcollegeyes</td>
                  <td data-row="7" data-col="2">0.837***</td>
                  <td data-row="7" data-col="3">0.837***</td>
                  <td data-row="7" data-col="4">0.837***</td>
                </tr>
                <tr>
                  <td data-row="8" data-col="1"></td>
                  <td data-row="8" data-col="2">(0.070)</td>
                  <td data-row="8" data-col="3">(0.071)</td>
                  <td data-row="8" data-col="4">(0.076)</td>
                </tr>
                <tr>
                  <td data-row="9" data-col="1">mcollegeyes</td>
                  <td data-row="9" data-col="2">0.568***</td>
                  <td data-row="9" data-col="3">0.568***</td>
                  <td data-row="9" data-col="4">0.568***</td>
                </tr>
                <tr>
                  <td data-row="10" data-col="1"></td>
                  <td data-row="10" data-col="2">(0.079)</td>
                  <td data-row="10" data-col="3">(0.080)</td>
                  <td data-row="10" data-col="4">(0.081)</td>
                </tr>
                <tr>
                  <td data-row="11" data-col="1">incomehigh</td>
                  <td data-row="11" data-col="2">0.480***</td>
                  <td data-row="11" data-col="3">0.480***</td>
                  <td data-row="11" data-col="4">0.480***</td>
                </tr>
                <tr>
                  <td data-row="12" data-col="1"></td>
                  <td data-row="12" data-col="2">(0.058)</td>
                  <td data-row="12" data-col="3">(0.060)</td>
                  <td data-row="12" data-col="4">(0.062)</td>
                </tr>
                <tr>
                  <td data-row="13" data-col="1">regionwest</td>
                  <td data-row="13" data-col="2">-0.116+</td>
                  <td data-row="13" data-col="3">-0.116*</td>
                  <td data-row="13" data-col="4">-0.116</td>
                </tr>
                <tr>
                  <td data-row="14" data-col="1"></td>
                  <td data-row="14" data-col="2">(0.062)</td>
                  <td data-row="14" data-col="3">(0.059)</td>
                  <td data-row="14" data-col="4">(0.072)</td>
                </tr>
                <tr>
                  <td data-row="15" data-col="1">Num.Obs.</td>
                  <td data-row="15" data-col="2">4739</td>
                  <td data-row="15" data-col="3">4739</td>
                  <td data-row="15" data-col="4">4739</td>
                </tr>
                <tr>
                  <td data-row="16" data-col="1">R2</td>
                  <td data-row="16" data-col="2">0.111</td>
                  <td data-row="16" data-col="3">0.111</td>
                  <td data-row="16" data-col="4">0.111</td>
                </tr>
                <tr>
                  <td data-row="17" data-col="1">R2 Adj.</td>
                  <td data-row="17" data-col="2">0.110</td>
                  <td data-row="17" data-col="3">0.110</td>
                  <td data-row="17" data-col="4">0.110</td>
                </tr>
                <tr>
                  <td data-row="18" data-col="1">AIC</td>
                  <td data-row="18" data-col="2">18421.1</td>
                  <td data-row="18" data-col="3">18421.1</td>
                  <td data-row="18" data-col="4">18421.1</td>
                </tr>
                <tr>
                  <td data-row="19" data-col="1">BIC</td>
                  <td data-row="19" data-col="2">18472.8</td>
                  <td data-row="19" data-col="3">18472.8</td>
                  <td data-row="19" data-col="4">18472.8</td>
                </tr>
                <tr>
                  <td data-row="20" data-col="1">Log.Lik.</td>
                  <td data-row="20" data-col="2">-9202.551</td>
                  <td data-row="20" data-col="3">-9202.551</td>
                  <td data-row="20" data-col="4">-9202.551</td>
                </tr>
                <tr>
                  <td data-row="21" data-col="1">F</td>
                  <td data-row="21" data-col="2">98.196</td>
                  <td data-row="21" data-col="3">106.909</td>
                  <td data-row="21" data-col="4"></td>
                </tr>
                <tr>
                  <td data-row="22" data-col="1">RMSE</td>
                  <td data-row="22" data-col="2">1.69</td>
                  <td data-row="22" data-col="3">1.69</td>
                  <td data-row="22" data-col="4">1.69</td>
                </tr>
                <tr>
                  <td data-row="23" data-col="1">Std.Errors</td>
                  <td data-row="23" data-col="2">IID</td>
                  <td data-row="23" data-col="3">HC0</td>
                  <td data-row="23" data-col="4">NeweyWest</td>
                </tr>
        </tbody>
      </table>
    </div>
<!-- hack to avoid NA insertion in last line -->
```

:::
:::


We note that the coefficient for *distance* is statistically significant
at the 1% level and we thus reject the null hypothesis of zero. We must
take care to note that distance is measured in 10 mile increments, that
is, 1 = 10 miles, 10 = 100 miles. We find that a 10 mile increase in
distance reduces the predicted years of education by 0.041 years.

It appears that having a parent who was a college graduate had a
positive influence on educational attainment. We reject the null
hypothesis of zero for the coefficient for whether the father graduated
college at the 1% level of statistical significance. We also reject the
null hypothesis of zero for the coefficient for whether the mother
graduated college at the 1% level of statistical significance. If the
father graduated college, estimated educational attainment increased by
approximately 0.84 years, while if the mother graduated college,
estimated educational attainment increased by approximately 0.57 years.

We reject the null hypothesis that the coefficient for income above
25,000 is equal to zero at the 1% level of statistical significance. For
households with income above 25,000, predicted educational attainment
increased by approximately 0.48 years.

While the coefficient for the West region is statistically significant
in the uncorrected model, it is only weakly significant at the 10% level
in the second model, and insignificant in the third model. We thus
conclude that we do not have sufficient evidence to reject the null
hypothesis of zero for the coefficient for the West region as it appears
to be fragile to the choice of correction for heteroscedasticity.

We note that gender does not appear to statistically significantly
influence educational attainment as we cannot reject the null hypothesis
the coefficient is equal to zero.

