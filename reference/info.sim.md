# Simulated Information IV Dataset

A simulated dataset for demonstrating information IV methods. Contains
experimental data with randomized information treatments and Bayesian
updating.

## Usage

``` r
info.sim
```

## Format

A data.table with 500 observations and 9 variables:

- tau:

  Individual-specific treatment effect parameters (normalized to mean 1)

- alpha:

  Learning rates, negatively correlated with tau

- Z:

  Binary treatment assignment (0/1)

- signal:

  Information signal (1 for high signal, -1 for low signal)

- prior:

  Prior beliefs (affected by unobserved confounder V)

- posterior:

  Posterior beliefs after updating

- Y:

  Outcome variable

- Y0:

  Counterfactual outcome (pre-treatment)

- dX:

  Change in beliefs (posterior - prior)

- dY:

  Change in outcomes (Y - Y0)

## Source

Simulated data

## Details

This dataset is generated with endogeneity where unobserved factors (V,
U) affect both beliefs and outcomes. The true average partial effect
(APE) is 1. People with higher belief effects (tau) have lower learning
rates (alpha).

## Examples

``` r
data(info.sim)
head(info.sim)
#>          tau     alpha      Z signal      prior  posterior          Y
#>        <num>     <num> <lgcl>  <num>      <num>      <num>      <num>
#> 1: 1.1871082 0.4150711   TRUE      1 -1.3948654 -0.4008260 -1.5443077
#> 2: 1.4630353 0.3184235  FALSE     -1 -1.5366966 -1.3657998 -4.2523432
#> 3: 1.8690124 0.2225582  FALSE     -1 -0.5080131 -0.6175088 -1.0815787
#> 4: 0.4509035 0.8310383  FALSE     -1  1.6910304 -0.5453189  0.2467968
#> 5: 1.5743751 0.2874673   TRUE      1  1.1247404  1.0888816  3.0765761
#> 6: 0.4813998 0.8118558  FALSE     -1  0.3403998 -0.7478116 -2.1605352
#>            Y0          dX          dY
#>         <num>       <num>       <num>
#> 1: -2.7243400  0.99403941  1.18003232
#> 2: -4.5023712  0.17089677  0.25002800
#> 3: -0.8769299 -0.10949571 -0.20464883
#> 4:  1.2551745 -2.23634927 -1.00837778
#> 5:  3.1330313 -0.03585879 -0.05645518
#> 6: -1.6366705 -1.08821137 -0.52386476

# Check the structure
str(info.sim)
#> Classes ‘data.table’ and 'data.frame':   500 obs. of  10 variables:
#>  $ tau      : num  1.187 1.463 1.869 0.451 1.574 ...
#>  $ alpha    : num  0.415 0.318 0.223 0.831 0.287 ...
#>  $ Z        : logi  TRUE FALSE FALSE FALSE TRUE FALSE ...
#>  $ signal   : num  1 -1 -1 -1 1 -1 1 -1 -1 1 ...
#>  $ prior    : num  -1.395 -1.537 -0.508 1.691 1.125 ...
#>  $ posterior: num  -0.401 -1.366 -0.618 -0.545 1.089 ...
#>  $ Y        : num  -1.544 -4.252 -1.082 0.247 3.077 ...
#>  $ Y0       : num  -2.724 -4.502 -0.877 1.255 3.133 ...
#>  $ dX       : num  0.994 0.1709 -0.1095 -2.2363 -0.0359 ...
#>  $ dY       : num  1.18 0.25 -0.2046 -1.0084 -0.0565 ...
#>  - attr(*, ".internal.selfref")=<pointer: (nil)> 

# Summary statistics
summary(info.sim)
#>       tau               alpha            Z               signal  
#>  Min.   :0.001161   Min.   :0.2131   Mode :logical   Min.   :-1  
#>  1st Qu.:0.527604   1st Qu.:0.3144   FALSE:250       1st Qu.:-1  
#>  Median :0.999693   Median :0.5002   TRUE :250       Median : 0  
#>  Mean   :1.000000   Mean   :0.5508                   Mean   : 0  
#>  3rd Qu.:1.476579   3rd Qu.:0.7822                   3rd Qu.: 1  
#>  Max.   :1.921508   Max.   :1.0000                   Max.   : 1  
#>      prior             posterior              Y                 Y0          
#>  Min.   :-3.079675   Min.   :-2.25665   Min.   :-5.6734   Min.   :-7.09037  
#>  1st Qu.:-0.704907   1st Qu.:-0.75815   1st Qu.:-1.4893   1st Qu.:-1.73570  
#>  Median : 0.002541   Median :-0.08406   Median :-0.1045   Median : 0.04646  
#>  Mean   :-0.005173   Mean   :-0.04358   Mean   :-0.1032   Mean   :-0.11653  
#>  3rd Qu.: 0.714311   3rd Qu.: 0.71187   3rd Qu.: 1.1971   3rd Qu.: 1.34209  
#>  Max.   : 3.092885   Max.   : 2.22197   Max.   : 6.8904   Max.   : 7.57786  
#>        dX                  dY           
#>  Min.   :-3.258083   Min.   :-1.645117  
#>  1st Qu.:-0.505976   1st Qu.:-0.390024  
#>  Median : 0.008641   Median : 0.006762  
#>  Mean   :-0.038409   Mean   : 0.013305  
#>  3rd Qu.: 0.490082   3rd Qu.: 0.422442  
#>  Max.   : 2.481694   Max.   : 2.014024  
```
