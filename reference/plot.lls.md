# Plot LLS Results

Takes the output of the lls function and creates a plot with binned
estimates and confidence intervals, similar to the plots in
Balla-Elliott (2025).

## Usage

``` r
# S3 method for class 'lls'
plot(
  x,
  nbins = 4,
  colors = NULL,
  text_size = 8,
  size_scale = 3,
  add_se = TRUE,
  add_hline = TRUE,
  noplot = FALSE,
  abs.x = FALSE
)
```

## Arguments

- x:

  An lls object or lls result list

- nbins:

  Number of bins for binned estimates (default: 4)

- colors:

  Color palette for the plot (default: uses met.brewer "Lakota")

- text_size:

  Size of text elements (default: 8)

- size_scale:

  Scale factor for plot elements (default: 3)

- add_se:

  Whether to add standard error bars (default: TRUE)

- add_hline:

  Whether to add horizontal line at APE (default: TRUE)

- noplot:

  Whether to return the binned data instead of a plot (default: FALSE)

- abs.x:

  Whether to use absolute value of x variable for binning and plotting
  (default: FALSE)

## Value

A ggplot2 object
