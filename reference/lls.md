# Local Least Squares (LLS) Estimation

Implements the Local Least Squares (LLS) estimator for identifying
causal effects in information provision experiments.

## Usage

``` r
lls(...)

iv.lls(
  dat,
  y,
  x,
  r,
  weights = NULL,
  FE.fml = NULL,
  control.fml = NULL,
  bandwidth = NULL,
  normalize.r = TRUE,
  r.support.points = nrow(dat),
  bootstrap = FALSE,
  bootstrap.bayesian = TRUE,
  bootstrap.n = 200,
  ci.level = 0.95,
  cluster = NULL,
  n.cores = availableCores(omit = 2),
  kernel = epanechnikov,
  auto.trim = TRUE,
  sd.trim = 0.02,
  pointmass.zero = FALSE,
  trim.zero = NULL
)

panel.lls(
  dat,
  dy,
  dx,
  weights = NULL,
  FE.fml = NULL,
  control.fml = NULL,
  bandwidth = NULL,
  normalize.x = FALSE,
  normalize.r = TRUE,
  r.support.points = nrow(dat),
  bootstrap = FALSE,
  bootstrap.bayesian = TRUE,
  bootstrap.n = 200,
  ci.level = 0.95,
  cluster = NULL,
  n.cores = availableCores(omit = 2),
  kernel = epanechnikov,
  auto.trim = TRUE,
  sd.trim = 0.02,
  pointmass.zero = FALSE,
  trim.zero = NULL
)
```

## Arguments

- ...:

  Arguments passed on to
  [`lls.internal`](https://dballaelliott.github.io/lls/reference/lls.internal.md)

  `mode`

  :   Character string specifying estimation mode: "iv" (default) or
      "panel"

- dat:

  A data frame or data.table containing the analysis data

- y:

  Character string specifying the outcome variable name in \`dat\` (for
  IV mode)

- x:

  Character string specifying the treatment variable name in \`dat\`
  (for IV mode)

- r:

  Optional; name of variable in \`dat\` giving ranks (IV mode) or vector
  of ranks directly if pre-computed (default: NULL)

- weights:

  Optional; character string specifying the weights variable name in
  \`dat\` (default: NULL)

- FE.fml:

  Optional; formula string for fixed effects specification (default:
  NULL)

- control.fml:

  Optional; formula string for control variables (default: NULL)

- bandwidth:

  Numeric; bandwidth value, if NULL will use rule-of-thumb (default:
  NULL)

- normalize.r:

  Logical; whether to normalize the running variable r to be equally
  spaced between 0 and 1 within sign groups (default: TRUE)

- r.support.points:

  Integer; number of support points to use in estimation (default:
  nrow(dat))

- bootstrap:

  Logical; whether to perform bootstrap inference (default: FALSE)

- bootstrap.bayesian:

  Logical; whether to use Bayesian bootstrap (TRUE) or standard
  bootstrap (default: TRUE)

- bootstrap.n:

  Integer; number of bootstrap iterations (default: 200)

- ci.level:

  Numeric; confidence level for intervals between 0 and 1 (default:
  0.95)

- cluster:

  Optional; character string specifying cluster variable for clustered
  standard errors (default: NULL)

- n.cores:

  Integer; number of cores to use for parallel processing (default:
  availableCores(omit=2))

- kernel:

  Function; kernel function to use (default: epanechnikov)

- auto.trim:

  Logical; if TRUE will trim the support points away from zero to avoid
  truncation. Sets trim.zero. (default: TRUE)

- sd.trim:

  Numeric; trimming proportion for SD calculation in bootstrap (default:
  0.02, 0 or FALSE disables trimming)

- pointmass.zero:

  Logical; if TRUE, will use near(x,0) instead of kernel(x,0) (default:
  FALSE)

- trim.zero:

  Numeric; threshold for trimming support points away from near zero
  values (default: NULL, uses bandwidth/2)

- dy:

  Character string specifying the outcome variable name in \`dat\` (for
  panel mode)

- dx:

  Character string specifying the treatment variable name in \`dat\`
  (for panel mode)

- normalize.x:

  Logical; whether to use ranks of x instead of raw values (default:
  FALSE; panel mode only)

## Value

A list with elements: coef, se, micro.dt, ci.percentile, ci.normal,
bs.ests, bs.micro, bandwidth

## Aliases

- `iv.lls`:

  Alias for `lls(..., mode = "iv")`

- `panel.lls`:

  Alias for `lls(..., mode = "panel")`
