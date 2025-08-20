#' Internal: Epanechnikov kernel function
#' @keywords internal
epanechnikov <- function(x) {
  3 / 4 * (1 - x^2) * (abs(x) <= 1)
}

#' Internal: trimmed standard deviation helper
#' @keywords internal
trunc.sd <- function(x, trim = 0.02) {
  if (trim == 0) return(sd(x, na.rm = TRUE))
  x <- x[!is.na(x)]
  x <- x[abs(x) < quantile(abs(x), probs = 1 - trim)]
  return(sd(x))
}

#' Local Least Squares (LLS) Estimation
#' 
#' If you're reading this, you've found the pre-release version of this package. Please check this repo frequently for updates. Please open an issue if you find a bug!
#' 
#' @name lls
#' @section Aliases:
#' \describe{
#'   \item{\code{iv.lls}}{Alias for \code{lls(..., mode = "iv")}}
#'   \item{\code{panel.lls}}{Alias for \code{lls(..., mode = "panel")}}
#' }
#' @return A list with elements: coef, se, micro.dt, ci.percentile, ci.normal, bs.ests, bs.micro, bandwidth

NULL
#> NULL



#' @keywords internal
#' @rdname lls
#' @inheritDotParams lls.internal
#' @import data.table
#' @import fixest
#' @importFrom pbapply pblapply
#' @importFrom parallelly availableCores
#' @importFrom Matrix Diagonal
#' @importFrom collapse fmean fsum fquantile
#' @importFrom DirichletReg rdirichlet
#' @importFrom dplyr near
#' @importFrom stats approx as.formula coef lm optimize quantile sd setNames qnorm pnorm rnorm runif .lm.fit
lls <- function(...) {
    do.call(lls.internal, list(...))
}

#' IV Mode
#'
#' @rdname lls
#' @inheritParams lls.internal
#' @export
iv.lls <- function(
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
  ci.level = .95,
  cluster = NULL,
  n.cores = availableCores(omit = 2),
  kernel = epanechnikov,
  auto.trim = TRUE,
  sd.trim = 0.02,
  pointmass.zero = FALSE,
  trim.zero = NULL
) {
  lls.internal(
    dat = dat, y = y, x = x, r = r,
    weights = weights,
    FE.fml = FE.fml, control.fml = control.fml,
    bandwidth = bandwidth,
    normalize.x = FALSE,
    normalize.r = normalize.r,
    r.support.points = r.support.points,
    bootstrap = bootstrap,
    bootstrap.bayesian = bootstrap.bayesian,
    bootstrap.n = bootstrap.n,
    ci.level = ci.level,
    cluster = cluster,
    n.cores = n.cores,
    kernel = kernel,
    mode = "iv",
    auto.trim = auto.trim,
    sd.trim = sd.trim,
    pointmass.zero = pointmass.zero,
    trim.zero = trim.zero
  )
}

#' @rdname lls
#' @inheritParams lls.internal
#' @export
panel.lls <- function(
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
  ci.level = .95,
  cluster = NULL,
  n.cores = availableCores(omit = 2),
  kernel = epanechnikov,
  auto.trim = TRUE,
  sd.trim = 0.02,
  pointmass.zero = FALSE,
  trim.zero = NULL
) {
  lls.internal(
    dat = dat, dy = dy, dx = dx, r = NULL,
    weights = weights,
    FE.fml = FE.fml, control.fml = control.fml,
    bandwidth = bandwidth,
    normalize.x = normalize.x,
    normalize.r = normalize.r,
    r.support.points = r.support.points,
    bootstrap = bootstrap,
    bootstrap.bayesian = bootstrap.bayesian,
    bootstrap.n = bootstrap.n,
    ci.level = ci.level,
    cluster = cluster,
    n.cores = n.cores,
    kernel = kernel,
    mode = "panel",
    auto.trim = auto.trim,
    sd.trim = sd.trim,
    pointmass.zero = pointmass.zero,
    trim.zero = trim.zero
  )
}

#' Internal LLS function that handles estimation.
#'
#' Don't call this -- use the exported `lls(), iv.lls(), panel.lls()` functions instead.
#' @keywords internal
#'
#' @param dat A data frame or data.table containing the analysis data
#' @param y Character string specifying the outcome variable name in `dat` (for IV mode)
#' @param x Character string specifying the treatment variable name in `dat` (for IV mode)
#' @param dy Character string specifying the outcome variable name in `dat` (for panel mode)
#' @param dx Character string specifying the treatment variable name in `dat` (for panel mode)
#' @param r Optional; name of variable in `dat` giving ranks (IV mode) or vector of ranks directly if pre-computed (default: NULL)
#' @param weights Optional; character string specifying the weights variable name in `dat` (default: NULL)
#' @param FE.fml Optional; formula string for fixed effects specification (default: NULL)
#' @param control.fml Optional; formula string for control variables (default: NULL)
#' @param bandwidth Numeric; bandwidth value, if NULL will use rule-of-thumb (default: NULL)
#' @param normalize.x Logical; whether to use ranks of x instead of raw values (default: FALSE; panel mode only)
#' @param normalize.r Logical; whether to normalize the running variable r to be equally spaced between 0 and 1 within sign groups (default: TRUE)
#' @param r.support.points Integer; number of support points to use in estimation (default: nrow(dat))
#' @param bootstrap Logical; whether to perform bootstrap inference (default: FALSE)
#' @param bootstrap.bayesian Logical; whether to use Bayesian bootstrap (TRUE) or standard bootstrap (default: TRUE)
#' @param bootstrap.n Integer; number of bootstrap iterations (default: 200)
#' @param ci.level Numeric; confidence level for intervals between 0 and 1 (default: 0.95)
#' @param cluster Optional; character string specifying cluster variable for clustered standard errors (default: NULL)
#' @param n.cores Integer; number of cores to use for parallel processing (default: availableCores(omit=2))
#' @param kernel Function; kernel function to use (default: epanechnikov)
#' @param mode Character string specifying estimation mode: "iv" (default) or "panel"
#' @param sd.trim Numeric; trimming proportion for SD calculation in bootstrap (default: 0.02, 0 or FALSE disables trimming)
#' @param auto.trim Logical; if TRUE will trim the support points away from zero to avoid truncation. Sets trim.zero. (default: TRUE)
#' @param pointmass.zero Logical; if TRUE, will use near(x,0) instead of kernel(x,0) (default: FALSE)
#' @param trim.zero Numeric; threshold for trimming support points away from near zero values (default: NULL, uses bandwidth/2)
#'
lls.internal <- function(dat, y = NULL, x = NULL, dy = NULL, dx = NULL,
    r = NULL, # ranks directly can be passed in
    weights = NULL,
    FE.fml = NULL,
    control.fml = NULL,
    # rule of thumb bandwidth
    pointmass.zero = FALSE,
    trim.zero = NULL,
    bandwidth = NULL,
    # how should we grid up the x variable
    normalize.x = FALSE, # only for panel mode
    # if TRUE, will use the ranks of x instead of the raw values
    normalize.r = TRUE, # if TRUE, will normalize r to be equally spaced between 0 and 1
    r.support.points = nrow(dat),
    # CI options,
    bootstrap = FALSE,
    # by default, it should do the bayesian bootstrap
    bootstrap.bayesian = TRUE,
    bootstrap.n = 200,
    ci.level = .95,
    cluster = NULL,
    n.cores = availableCores(omit=2),
    kernel  = epanechnikov,
    mode = "iv", # "iv" or "panel"
    auto.trim = TRUE, # if TRUE, will trim the support points to be above the bandwidth
    sd.trim = 0.02, # trim for SD calculation in bootstrap
    user.call.level = TRUE # undocumented, used for recursion in the function
    ) {

    # throw an error if the bandwidth is negative
    if (!is.null(bandwidth) && bandwidth < 0) {
        stop("Bandwidth cannot be negative")
    }

    # Handle mode-specific variable assignment
    panel.mode <- FALSE
    if (mode == "panel") {
        panel.mode <- TRUE
        if (is.null(dy) || is.null(dx)) {
            stop("In panel mode, both 'dy' and 'dx' arguments must be provided")
        }
        if (!is.null(y) || !is.null(x)) {
            warning("In panel mode, 'y' and 'x' arguments are ignored. Use 'dy' and 'dx' instead.")
        }
        if (!is.null(r)) {
            stop("In panel mode, 'r' argument is not used. Double check that the mode is set correctly.")
        }
        # Set y and x to panel mode variables for internal use
        y <- dy
        x <- dx
    } else if (mode == "iv") {
        if (is.null(y) || is.null(x) || is.null(r)) {
            stop("In IV mode, 'y', 'x', and 'r' arguments must be provided")
        }
        if (!is.null(dy) || !is.null(dx)) {
            warning("In IV mode, 'dy' and 'dx' arguments are ignored. Use 'y' and 'x' instead.")
        }
        if (normalize.x){
            warning("In IV mode, 'normalize.x' only has a meaning in panel mode.\nLocal regressions in iv mode condition on R.")
            warning("If you want to use the ranks of r, use `normalize.r = TRUE`.")
        }
    } else {
        stop("Mode must be either 'iv' or 'panel'")
    }


    if(bootstrap & n.cores > 1){
        initial.nthreads <- getFixest_nthreads()
        setFixest_nthreads(1)
    }

    mc_for_recursion <- match.call()
    mc_for_recursion$user.call.level <- FALSE

    if (!bootstrap & user.call.level) {
        message("Set 'bootstrap = TRUE' for standard errors and confidence intervals.")
    }
    # only bootstrap in the outermost call
    mc_for_recursion$bootstrap <- FALSE

    # If bandwidth is not provided, use rule of thumb
    if (is.null(bandwidth)) {
        if (user.call.level) {
            message("No bandwidth provided, using defaults. You should experiment with the bandwidth.")
        }

        if (!normalize.r & !normalize.x) {
            bandwidth <- .9 / sqrt(12) / nrow(dat)^(1 / 5)
            if(user.call.level) {
                message("Using rule-of-thumb bandwidth: ", round(bandwidth, 4))
            }
        }
        else if (normalize.r | normalize.x) {
            bandwidth <- 0.05 # default bandwidth for normalized variables
            if(user.call.level) {
                message("Using default bandwidth of 0.05 for normalized variables.\n",
                        "This is roughly 5% of the data in each local regression.")
            }
        }
    }

    ## end of cross validation
    if (is.null(trim.zero) & pointmass.zero)  trim.zero <- .Machine$double.eps^0.5
    if (is.null(trim.zero) & !pointmass.zero) trim.zero <- bandwidth / 2


    # --------------------------- BEGIN MAIN ESTIMATES --------------------------- #
    dt <- data.table::as.data.table(dat)

    ## set up the data/ranks that we're going to use
    # x.trt.min <-  .Machine$double.eps^0.5 # default tolerance for e.g. near()
    # if (!pointmass.zero) x.trt.min <- bandwidth / 2 # trimming away from zero
    if (is.null(r) & mode == "panel") {
        dt$LLS.INT.R <- dt[[x]]
    } else {
        dt$LLS.INT.R <- dt[[r]]
    }


    if (normalize.r | normalize.x) {
        dt[, r := sign(LLS.INT.R) * frank(abs(LLS.INT.R), ties.method = "dense") / .N, by = sign(LLS.INT.R)]
    }
    else {
        dt$r <- dt$LLS.INT.R
    }


    suppR <- sort(unique(dt$r))
    if (r.support.points < length(suppR)) {
        probs <- seq(0, 1, length.out = r.support.points + 2L)
        probs <- probs[probs > 0 & probs < 1]
        if (normalize.x) {
            vec <- dt[abs(r) >= .Machine$double.eps^0.5, .SD, .SDcols = c("r", weights)]
        } else {
            vec <- dt[abs(r) >= (trim.zero + bandwidth/2), .SD, .SDcols = c("r", weights)]
        }
        if (!is.null(r)){
            vec <- dt[, .SD, .SDcols = c(r, weights)]
        setnames(vec,r, "r")
        }
        suppR <- fquantile(vec$r, probs = probs, w = if (!is.null(weights)) vec[[weights]] else NULL)
    }

    if (auto.trim){
        if (!normalize.x) suppR <- suppR[abs(suppR) > trim.zero]
            # also need to trim away from zero
        if(pointmass.zero | !(panel.mode)) suppR <- suppR[abs(suppR) > bandwidth/2]
        else suppR <- suppR[abs(suppR) > bandwidth]
    }



    ## now we have the ranks, we can estimate the tau
    # message("Estimating tau at support points", print(suppR))
    ## estimate the tau at the support points
    b.eff <- se.eff <- TAU <- tau_dt <- NULL
    sd <- ci.percentile <- ci.normal <- bs_ests <- bs.micro <- NULL
    dt[, rowid := 1:.N]
        tau_dt <- lapply(suppR, function(r) {
            # subset around the point r and the point 0
            wt1 <- kernel((r - dt$r) / (bandwidth / 2))
            wt0 <- NULL
            if (pointmass.zero) { # implies panel mode is true
                wt0 <- 1 * dplyr::near(dt$r, 0) #
            } else if (panel.mode) {
                wt0 <- kernel(dt$r / (bandwidth / 2)) #
            } else {
                # panel mode is false
                wt0 <- rep(0, length(dt$r)) # no control group
            }
            wt0_denom <- fifelse(pointmass.zero | panel.mode, fsum(wt0), 1)
            wt <- wt1 / fsum(wt1) + wt0 / wt0_denom
            # rescale by the size of the treated group
            wt <- wt * fsum(wt1) / fsum(wt)


            if (!is.null(weights)) wt <- wt * dt[[weights]]

            touse <- dt[wt > 0]
            wt.touse <- wt[wt > 0]

            if (nrow(touse) == 0) {
                return(data.table(tau_r = NA_real_, tau_wt = NA_real_, x = NA_real_, r = r))
            }

            # if we have FE, let's just use feols
            tau_r <- NA_real_
            if (!is.null(FE.fml)) {
                # message("FE.fml is ", FE.fml)
                # make the fml for feols
                if (!is.null(control.fml)) control.fml <- paste("+", control.fml)
                if (!is.null(FE.fml)) FE.fml <- paste("|", FE.fml)

                fml <- paste(y, "~", x, control.fml, FE.fml) |> as.formula()


                # try catch in case we can't estimate this with FE
                tryCatch(
                    {
                        tau_r <- coef(feols(fml, data = touse, weights = wt.touse))[[x]]
                    },
                    error = function(e) {
                        message(e, "\n")
                    }
                )
            } else if (!is.null(control.fml)) { # implies no FE, but control

                    fml <- paste(y, "~", x, "+", control.fml) |> as.formula()
                    tau_r <- tryCatch({
                        stats::coef(stats::lm(fml, data = touse, weights = wt.touse))[[x]]
                    }, error = function(e) NA_real_)
            } else { # no FE, or controls, use .lm.fit
                xM <- as.matrix(cbind(1, dt[[x]])) * sqrt(wt)
                tau_r <- .lm.fit(
                    xM, # x
                    sqrt(wt) * dt[[y]] # Y
                )$coefficients[2] # way faster
            }


            x.point <- fmean(dt[[x]], w = wt1)

            data.table(tau_r = tau_r, tau_wt = fsum(wt1), x = x.point, r = r)
        }) |> rbindlist(use.names = TRUE, fill = TRUE)

        # merge the tau_dt with the original data
        # unless we subsampled, in which case the weights have already entered
        # in picking support points
        if (r.support.points == nrow(dat)) {
            tau_dt <- tau_dt[dt[, .(r = r)], on = "r"]
            # return the average of the taus
            if (nrow(tau_dt) > 0) {
                TAU <- fmean(tau_dt$tau_r, w = if (!is.null(weights)) dt[[weights]] else NULL)
            } else {
                TAU <- NA
            }
        }
        else {
            # if we subsampled, just return the average of the taus
            if (nrow(tau_dt) > 0) {
                TAU <- fmean(tau_dt$tau_r)
            } else {
                TAU <- NA
            }
        }

    # tau_dt <- tau_dt[dt[, .(r = r)], on = "r"]



    # ---------------------------- END MAIN ESTIMATES ---------------------------- #


    if (bootstrap) {
        message("Bootstrapping with ", bootstrap.n, " iterations")

        bs_est_list <- pblapply(1:bootstrap.n, \(bs_id){
            dt <- copy(dat)
            
            if (!is.null(cluster)){
                dt$cluster <- dt[[cluster]]
            }

            # # Bayesian bootstrap: draw weights from a dirichlet distribution
            # bayesian bootstrap
            bs_weightname <- weights
            if (bootstrap.bayesian) {
                if (!is.null(cluster)) { # clustering
                    # cluster bootstrap
                    if (bs_id == 1) message("Clustered by ", cluster)
                    dt[, cluster_id := .GRP, by = cluster]
                    # dt[, cluster_id := .GRP, by = get(cluster)]
                    cluster.dt <- dt[, list(cluster_id)] |> unique()
                    cluster.dt$bsweights <- rdirichlet(n = 1, rep(1, nrow(cluster.dt))) |> as.vector()
                    dt <- dt[cluster.dt, on = "cluster_id"]
                }
                else { # no clustering
                    dt$bsweights <- rdirichlet(n = 1, rep(1, nrow(dt))) |> as.vector()
                }

                if (!is.null(weights)) dt$bsweights <- dt$bsweights * dt[[weights]]

                bs_weightname <- "bsweights"
                mc_for_recursion[['weights']] <- bs_weightname

            } else { # standard bootstrap
                if (!is.null(cluster)) { # clustering
                    if (bs_id == 1) message("Clustered by ", cluster)
                    dt[, cluster_id := .GRP, by = cluster]

                    if (!is.null(weights)) {
                        cluster.dt <- dt[, .(wt = fsum(get(weights))), keyby = cluster_id]
                        cluster.dt <- cluster.dt[sample(.N, replace = T, prob = cluster.dt$wt)]
                    }
                    else {
                        cluster.dt <- dt[, list(cluster_id)] |> unique()
                        cluster.dt <- cluster.dt[sample(.N, replace = T)]
                    }
                    dt <- dt[cluster.dt, on = "cluster_id", allow.cartesian=TRUE]
                }
                else { # no clustering
                    if (!is.null(weights)) dt <- dt[sample(.N, replace =T, prob =dt[[weights]])]
                    else dt <- dt[sample(.N,replace = T)]
                }

            }


            mc_for_recursion[['dat']] <- dt

            # need to make sure we don't display any messages
            # sink(file = file(tempfile(pattern = bs_id), open = "wt"), type = "message")
            bs.out <- eval(mc_for_recursion)
            bs_est <- bs.out$coef |> unname()

            bs_micro <- bs.out$micro.dt
            bs_micro[, bs_id := bs_id]

            return(
                list(
                    bs_est = bs_est,
                    bs_micro = bs_micro
                )
            )
        }, cl = n.cores)
        # })

        bs.micro <- lapply(bs_est_list, \(x) x$bs_micro) |> rbindlist(use.names = TRUE, fill = TRUE)
        bs_ests <- lapply(bs_est_list, \(x) x$bs_est) |> unlist(use.names = FALSE)

        valid.bs <- bs_ests[is.numeric(bs_ests)]
        if (length(valid.bs) != bootstrap.n) {
            message("Warning: ", bootstrap.n - length(valid.bs), " bootstrap estimates were not numeric")
        }
        bs_ests <- valid.bs
        # print("can we get the SD?")
        # add cite: bruce hansen suggests trimming the top and bottom 1% of the bootstrap estimates
        sd <- trunc.sd(bs_ests, trim = sd.trim)

        # print(sd)
        # print("can we get the CI?")
        # print(bs_ests)
    ci.percentile <- stats::quantile(bs_ests, probs = c((1 - ci.level) / 2, 1 - (1 - ci.level) / 2), names = FALSE, type = 7)
        ci.normal <- c(TAU + qnorm((1 - ci.level) / 2) * sd, TAU + qnorm(1-(1 - ci.level) / 2) * sd)

    }

    if (bootstrap & n.cores > 1) {
        setFixest_nthreads(initial.nthreads)
    }

    result <- list(
        coef = TAU, se = sd,
        micro.dt = data.table(tau_dt),
        ci.percentile = ci.percentile,
        ci.normal = ci.normal,
        bs.ests = bs_ests,
        bs.micro = bs.micro,
        bandwidth = bandwidth
    )

    # Set class for S3 methods
    class(result) <- "lls"
    return(result)

}


#' Print method for lls objects
#'
#' @param lls.return An lls object
#' @param ... Additional arguments (ignored)
#' @export
#' @method print lls
print.lls <- function(lls.return, ...) {
    x <- lls.return
  cat("Local Least Squares (LLS) Estimation\n")
  cat("====================================\n\n")

  # Main coefficient
  cat("Average Partial Effect (APE):\n")
  if (!is.null(x$coef) && !is.na(x$coef)) {
    cat(sprintf("  Estimate: %8.4f\n", x$coef))

    if (!is.null(x$se) && !is.na(x$se)) {
        cat(sprintf("  Std. Err: %8.4f\n", x$se))
        cat(sprintf("  t-value:  %8.4f\n", x$coef / x$se))

        p_val <- 2 * (1 - pnorm(abs(x$coef / x$se)))
        cat(sprintf("  p-value:   %s\n", ifelse(p_val < 0.001, "<0.001", sprintf("%8.4f", p_val))))
    }
  } else {
    cat("  Estimate: Not available\n")
  }

  # Confidence intervals
  if (!is.null(x$ci.normal) && length(x$ci.normal) == 2) {
    cat(sprintf("\nNormal CI (95%%): [%7.4f, %7.4f]\n", x$ci.normal[1], x$ci.normal[2]))
  }

  if (!is.null(x$ci.percentile) && length(x$ci.percentile) == 2) {
    cat(sprintf("Percentile CI (95%%): [%7.4f, %7.4f]\n", x$ci.percentile[1], x$ci.percentile[2]))
  }

  # Estimation details
  cat("\nEstimation Details:\n")
  if (!is.null(x$bandwidth)) {
    cat(sprintf("  Bandwidth: %8.4f\n", x$bandwidth))
  }

  if (!is.null(x$bs.ests)) {
    cat(sprintf("  Bootstrap reps: %d\n", length(x$bs.ests)))
  }

  if (!is.null(x$micro.dt)) {
    cat(sprintf("  Observations: %d\n", nrow(x$micro.dt)))
    cat(sprintf("  Support points: %d\n", unique(x$micro.dt$r) |> length()))
  }

  invisible(x)
}


#' Summary method for lls objects
#'    just an alias for print
#' @rdname print.lls
#' @param object An lls object
#' @param ... Additional arguments (ignored)
#' @export
#' @method summary lls
summary.lls <- function(lls.return, ...) {

  print(lls.return)
  invisible(lls.return)
}


