#' Plot LLS Results
#'
#' @description
#' Takes the output of the lls function and creates a plot with binned estimates
#' and confidence intervals, similar to the plots in Balla-Elliott (2025). 
#'
#' @param x An lls object or lls result list
#' @param nbins Number of bins for binned estimates (default: 4)
#' @param colors Color palette for the plot (default: uses met.brewer "Lakota")
#' @param text_size Size of text elements (default: 8)
#' @param size_scale Scale factor for plot elements (default: 3)
#' @param add_se Whether to add standard error bars (default: TRUE)
#' @param add_hline Whether to add horizontal line at APE (default: TRUE)
#' @param noplot Whether to return the binned data instead of a plot (default: FALSE)
#' @param abs.x Whether to use absolute value of x variable for binning and plotting (default: FALSE)
#'
#' @return A ggplot2 object
#' @export
#' @method plot lls
#' @import ggplot2
#' @import data.table
#' @import MetBrewer
#' @import collapse
#' @importFrom stats quantile qnorm
plot.lls <- function(x, 
                     nbins = 4,
                     colors = NULL,
                     text_size = 8,
                     size_scale = 3,
                     add_se = TRUE,
                     add_hline = TRUE,
                     noplot = FALSE,
                     abs.x = FALSE) {
  
  # Extract micro data from lls result
  micro_dt <- data.table::copy(x$micro.dt)
  
  # Check if we have bootstrap micro data for confidence intervals
  if (is.null(x$bs.micro) && add_se) {
    warning("No bootstrap micro data found. Setting add_se = FALSE \n Rerun lls with bootstrap = TRUE to obtain standard errors.")

    add_se <- FALSE
  }
  
  # Create bins for the estimates
  if (nrow(micro_dt) > 0) {
    # Create bins based on quantiles of r or x
    x_var <- if ("r" %in% names(micro_dt)) "r" else "x"
    
    # Use absolute value if requested
    if (abs.x) {
      micro_dt[, (x_var) := abs(get(x_var))]
    }
    
    # Create bins
    probs <- (0:nbins) / nbins
    bins <- micro_dt[, stats::quantile(get(x_var), probs = probs, na.rm = TRUE)]
    bins[1] <- min(bins[1], 0)
    micro_dt[, xbin := cut(get(x_var), bins, include.lowest = TRUE, right = TRUE)]
    
    # Create binned estimates
    ape.est <- x$coef
    bin_dt <- micro_dt[, .(
      tau_r = collapse::fmean(tau_r),
      x = collapse::fmean(get(x_var)),
      ape = ape.est
    ), by = xbin]
    
    # If we have bootstrap data, add confidence intervals
    if (!is.null(x$bs.micro) && add_se) {
      bs_micro <- data.table::copy(x$bs.micro)
      
      if (abs.x) {
        bs_micro[, (x_var) := abs(get(x_var))]
      }
      
      bs_micro[, xbin := cut(get(x_var), bins, include.lowest = TRUE, right = TRUE)]
      bs_bin_dt <- bs_micro[, .(tau_r = collapse::fmean(tau_r)), by = .(xbin, bs_id)]
      se_dt <- bs_bin_dt[, .(se = trunc_sd(tau_r)), by = xbin]
      bin_dt <- merge(bin_dt, se_dt, by = "xbin", all.x = TRUE)
      bin_dt[, ub := tau_r + 2 * se]
      bin_dt[, lb := tau_r - 2 * se]
    }
    
    if (is.null(colors)) {
      colors <- MetBrewer::met.brewer("Lakota", 5)[c(1, 3:5)]
    }

    if (noplot) return(bin_dt)

    p <- ggplot2::ggplot(bin_dt, ggplot2::aes(x = x, y = tau_r)) +
      ggplot2::geom_line(linewidth = size_scale / 2, alpha = 0.5, color = colors[1]) +
      ggplot2::geom_point(alpha = 1, size = 2 * size_scale, color = colors[1]) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = "solid") 

    if (add_hline && !is.null(x$coef)) {
      p <- p + ggplot2::geom_hline(
        yintercept = x$coef,
        linetype = "dashed",
        color = colors[1],
        linewidth = 3/4 * size_scale
      )
    }

    if (add_se && "se" %in% names(bin_dt)) {
      p <- p + ggplot2::geom_pointrange(
        ggplot2::aes(ymin = lb, ymax = ub),
        linewidth = 3/4 * size_scale,
        alpha = 0.8,
        color = colors[1]
      )
    }

    return(p)
  } else {
    stop("No micro data available for plotting")
  }
}

# Helper function for bin SEs
trunc_sd <- function(x, trim = 0.02) {
  x <- x[!is.na(x)]
  x <- x[abs(x) < stats::quantile(abs(x), probs = 1 - trim)]
  return(stats::sd(x))
}