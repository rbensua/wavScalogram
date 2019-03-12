#' @title Windowed scale index
#'
#' @description This function computes the windowed scale indices of a signal in the scale
#' interval \eqn{[s_0,s_1]}, for a given set of scale parameters \eqn{s_1} and taking
#' \eqn{s_0} as the minimum scale (see Benítez et al. 2010).
#'
#' The windowed scale index of a signal in the scale interval \eqn{[s_0,s_1]} centered at
#' time \eqn{tc} and with time windows radius \code{windowrad} is given by the quotient
#' \deqn{\frac{WS_{windowrad}(tc,s_{min})}{WS_{windowrad}(tc,s_{max})},}{WS_{windowrad}
#' (tc,s_{min})/WS_{windowrad}(tc,s_{max}),}
#' where \eqn{WS_{windowrad}} is the corresponding windowed scalogram with time windows
#' radius \code{windowrad}, \eqn{s_{max} \in [s_0,s_1]} is the smallest scale such that
#' \eqn{WS_{windowrad}(tc,s)\le WS_{windowrad}(tc,s_{max})} for all \eqn{s \in [s_0,s_1]},
#' and \eqn{s_{min} \in [s_{max},2s_1]} is the smallest scale such that
#' \eqn{WS_{windowrad}(tc,s_{min})\le WS_{windowrad}(tc,s)} for all
#' \eqn{s \in [s_{max},2s_1]}.
#'
#' @usage windowed_scale_index(signal,
#'                             dt = 1,
#'                             scales = NULL,
#'                             powerscales = TRUE,
#'                             s1 = NULL,
#'                             windowrad = NULL,
#'                             delta_t = NULL,
#'                             wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
#'                             wparam = NULL,
#'                             waverad = NULL,
#'                             border_effects = c("BE", "INNER", "PER", "SYM"),
#'                             makefigure = TRUE,
#'                             time_values = NULL,
#'                             figureperiod = TRUE,
#'                             xlab = "Time",
#'                             ylab = NULL,
#'                             main = "Windowed Scale Index")
#'
#' @param signal A vector containing the signal whose windowed scale indices are wanted.
#' @param dt Numeric. The time step of the signal.
#' @param scales A vector containing the wavelet scales at wich the windowed scalograms
#' are computed. This can be either a vector with all the scales, or (if
#' \code{powerscales} is TRUE) following Torrence and Compo 1998, a vector of three
#' elements with the minimum scale, the maximum scale and the number of suboctaves per
#' octave. If NULL, they are automatically computed.
#' @param powerscales Logical. If TRUE (default), construct power 2 scales from
#' \code{scales}. If \code{scales} is NULL, they are automatically computed.
#' @param s1 A vector containing the scales \eqn{s_1}. The windowed scale indices are
#' computed in the intervals \eqn{[s_0,s_1]}, where \eqn{s_0} is the minimum scale in
#' \code{scales}.
#' @param windowrad Integer. Time radius for the windows, measured in dt. By default,
#' it is set to \eqn{ceiling(length(signal) / 20)}.
#' @param delta_t Integer. Increment of time for the construction of windows central
#' times, measured in \code{dt}. By default, it is set to
#' \eqn{ceiling(length(signal) / 256)}.
#' @param wname A string, equal to "MORLET", "DOG", "PAUL", "HAAR" or "HAAR2". The
#' difference between "HAAR" and "HAAR2" is that "HAAR2" is more accurate but slower.
#' @param wparam The corresponding nondimensional parameter for the wavelet function
#' (Morlet, DoG or Paul).
#' @param waverad Numeric. The radius of the wavelet used in the computations for the cone
#' of influence. If it is not specified, it is asumed to be \eqn{\sqrt{2}} for Morlet and DoG,
#' \eqn{1/\sqrt{2}} for Paul and 0.5 for Haar.
#' @param border_effects A string, equal to "BE", "INNER", "PER" or "SYM", which indicates
#' how to manage the border effects which arise usually when a convolution is performed on
#' finite-lenght signals.
#' \itemize{
#' \item "BE": With border effects, padding time series with zeroes.
#' \item "INNER": Normalized inner scalogram with security margin adapted for each
#' different scale.
#' \item "PER": With border effects, using boundary wavelets (periodization of the
#' original time series).
#' \item "SYM": With border effects, using a symmetric catenation of the original time
#' series.
#' }
#' @param makefigure Logical. If TRUE (default), a figure with the windowed scale indices
#' is plotted.
#' @param time_values A numerical vector of length \code{length(signal)} containing custom
#' time values for the figure. If NULL (default), it will be computed starting at 0.
#' @param figureperiod Logical. If TRUE (default), periods are used in the figure instead
#' of scales.
#' @param xlab A string giving a custom X axis label.
#' @param ylab A string giving a custom Y axis label. If NULL (default) the Y label is
#' either "s1" or "Period of s1" depending on the value of \code{figureperiod} if
#' \code{length(s1) > 1}, or "Windowed Scale Index" if \code{length(s1) == 1}.
#' @param main A string giving a custom main title for the figure.
#'
#' @return A list with the following fields:
#' \itemize{
#' \item \code{wsi}: A matrix of size \code{length(tcentral)} x \code{length(s1)}
#' containing the values of the corresponding windowed scale indices.
#' \item \code{s1}: The vector of scales \eqn{s_1}.
#' \item \code{smax}: A matrix of size \code{length(tcentral)} x \code{length(s1)}
#' containing the scales \eqn{s_{max}}.
#' \item \code{smin}: A matrix of size \code{length(tcentral)} x \code{length(s1)}
#' containing the scales \eqn{s_{min}}.
#' \item \code{scalog_smax}: A matrix of size \code{length(tcentral)} x \code{length(s1)}
#' containing the values of the corresponding scalograms at scales \eqn{s_{max}}.
#' \item \code{scalog_smin}: A matrix of size \code{length(tcentral)} x \code{length(s1)}
#' containing the values of the corresponding scalograms at scales \eqn{s_{min}}.
#' \item \code{tcentral}: The vector of central times used in the computation of
#' \code{wsi}.
#' \item \code{fourierfactor}: A factor for converting scales into periods.
#' \item \code{coi_maxscale}: A vector of length \code{length(tcentral)} containing the
#' values of the maximum scale from which there are border effects.
#' }
#'
#' @importFrom fields image.plot
#' @importFrom colorRamps matlab.like
#' @importFrom graphics abline
#'
#' @examples
#' dt <- 0.1
#' time <- seq(0, 50, dt)
#' signal <- c(sin(pi * time), sin(pi * time / 2))
#' # First, we try with default s1 scales (a vector with a wide range of values for s1).
#' wsi_full <- windowed_scale_index(signal = signal, dt = dt, figureperiod = FALSE)
#' # Next, we choose a meaningful s1 value, greater than all relevant scales.
#' wsi <- windowed_scale_index(signal = signal, dt = dt, s1 = 4, figureperiod = FALSE)
#'
#' @section References:
#'
#' R. Benítez, V. J. Bolós, M. E. Ramírez. A wavelet-based tool for studying
#' non-periodicity. Comput. Math. Appl. 60 (2010), no. 3, 634-641.
#'
#' @export
#'

windowed_scale_index <-
  function(signal,
           dt = 1,
           scales = NULL,
           powerscales = TRUE,
           s1 = NULL,
           windowrad = NULL,
           delta_t = NULL,
           wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
           wparam = NULL,
           waverad = NULL,
           border_effects = c("BE", "INNER", "PER", "SYM"),
           makefigure = TRUE,
           time_values = NULL,
           figureperiod = TRUE,
           xlab = "Time",
           ylab = NULL,
           main = "Windowed Scale Index") {

  wname <- toupper(wname)
  wname <- match.arg(wname)

  if (is.null(waverad)) {
    if ((wname == "MORLET") || (wname == "DOG")) {
      waverad <- sqrt(2)
    } else if (wname == "PAUL") {
      waverad <- 1 / sqrt(2)
    } else { # HAAR
      waverad <- 0.5
    }
  }

  border_effects <- toupper(border_effects)
  border_effects <- match.arg(border_effects)

  nt <- length(signal)

  if (is.null(delta_t)){
    delta_t <- ceiling(nt / 256)
  }

  if (is.null(windowrad)) {
    windowrad <- ceiling(nt / 20)
  } else {
    windowrad <- min(windowrad,floor((nt - 1) / 2))
  }

  if (is.unsorted(s1)) {
    s1 <- sort(s1)
  }

  fourierfactor <- fourier_factor(wname = wname, wparam = wparam)

  if (is.null(scales)) {
    scmin <- 2 * dt / fourierfactor
    if (is.null(s1)) {
      scmax <- floor((nt - 2 * windowrad) / (2 * waverad)) * dt
    } else {
      scmax <- 2 * s1[length(s1)]
    }
    if (powerscales) {
      scales <- pow2scales(c(scmin, scmax, ceiling(256 / log2(scmax / scmin))))
    } else {
      scales <- seq(scmin, scmax, by = (scmax - scmin) / 256)
    }
  } else {
    ns <- length(scales)
    if (powerscales && ns == 3) {
      scales <- pow2scales(scales)
    } else {
      if (powerscales && ns != 3) {
        warning("The length of scales is not 3. Powerscales set to FALSE.")
        powerscales <- FALSE
      }
      if (is.unsorted(scales)) {
        warning("Scales were not sorted.")
        scales <- sort(scales)
      }
    }
  }

  ns <- length(scales)

  if (is.null(s1)) {
    s1null <- TRUE
    if (scales[ns] < 2 * scales[1]) {
      stop("We need larger scales.")
    }
    index_s1 <- which(scales <= scales[ns] / 2)
    s1 <- scales[index_s1]
    ns1 <- length(s1)
    index_2s1 <- numeric(ns1)
    for (j in 1:ns1) {
      index_2s1[j] <- max(which(scales <= 2 * s1[j]))
    }
  } else {
    s1null <- FALSE
    ns1 <- length(s1)
    if (s1[1] < scales[1]){
      stop("s1 must be >= the minimum scale.")
    }
    if (2 * s1[ns1] > scales[ns]) {
      stop("2*s1 must be <= the maximum scale.")
    }
    index_s1 <- numeric(ns1)
    index_2s1 <- numeric(ns1)
    for (j in 1:ns1) {
      index_s1[j] <- max(which(scales <= s1[j]))
      index_2s1[j] <- max(which(scales <= 2 * s1[j]))
    }
  }

  scales <- scales[1:index_2s1[ns1]]
  ns <- length(scales)
  if (ns < 4) {
    stop("We need more scales.")
  }

  wsc <-
    windowed_scalogram(signal = signal, dt = dt,
                       scales = scales, powerscales = FALSE,
                       windowrad = windowrad,
                       delta_t = delta_t,
                       wname = wname, wparam = wparam, waverad = waverad,
                       border_effects = border_effects,
                       energy_density = FALSE,
                       makefigure = FALSE)

  tcentraldt <- wsc$tcentral
  nwsi <- length(tcentraldt)

  wsi <- matrix(NA, nrow = nwsi, ncol = ns1)
  scalog_smax <- wsi
  scalog_smin <- wsi
  smax <- wsi
  smin <- wsi
  epsilon <- max(wsc$wsc, na.rm = TRUE) * 1e-6 # This is considered the "numerical zero"
  for (i in 1:nwsi) {
    # Si todo es NA hay que saltarlo:
    if (any(!is.na(wsc$wsc[i, index_2s1]))) {
      ni <- max(which(!is.na(wsc$wsc[i, index_2s1]))) # If border_effects = "INNER" there are NA in wsc.
      for (j in 1:ni) {
        scalog_smax[i, j] <- max(wsc$wsc[i, 1:index_s1[j]])
        if (scalog_smax[i, j] > epsilon) {
          index_smax <- which.max(wsc$wsc[i, 1:index_s1[j]])
          smax[i, j] <- scales[index_smax]

          scalog_smin[i, j] <- min(wsc$wsc[i, index_smax:index_2s1[j]])
          index_smin <- which.min(wsc$wsc[i, index_smax:index_2s1[j]])
          smin[i, j] <- scales[index_smax + index_smin - 1]
          wsi[i, j] <- scalog_smin[i, j] / scalog_smax[i, j]
        }
      }
    }
  }

  # COI
  coi_maxscale <- wsc$coi_maxscale / 2

  if (makefigure) {

    if (figureperiod) {
      Y <- fourierfactor * s1
      coi <- fourierfactor * coi_maxscale
      if (is.null(ylab)) ylab <- expression('Period of s'[1])
    } else {
      Y <- s1
      coi <- coi_maxscale
      if (is.null(ylab)) ylab  <- expression('s'[1])
    }

    if (is.null(time_values)) {
      X <- tcentraldt
    } else {
      if (length(time_values) != nt) {
        warning("Invalid length of time_values vector. Changing to default.")
        X <- tcentraldt
      } else {
        X <- time_values[floor(tcentraldt / dt)]
      }
    }

    if (length(Y) > 1) {
      wavPlot(
        Z = wsi,
        X = X,
        Y = Y,
        Ylog = powerscales,
        coi = coi,
        Xname = xlab,
        Yname = ylab,
        Zname = main
      )
    } else {
      if (is.null(ylab)) ylab <- "Windowed Scale Index"
      plot(X, wsi, type = "l", xlab = xlab, ylab = ylab, main = main, xaxt = "n")
      axis(side = 1, at = X[1 + floor((0:8) * (nwsi - 1) / 8)])
      abline(v = range(X[(coi > Y)]), lty = 2)
    }

  }

  return(list(
    wsi = wsi,
    s1 = s1,
    smax = smax,
    smin = smin,
    scalog_smax = scalog_smax,
    scalog_smin = scalog_smin,
    tcentral = tcentraldt,
    fourierfactor = fourierfactor,
    coi_maxscale = coi_maxscale
  ))
}
