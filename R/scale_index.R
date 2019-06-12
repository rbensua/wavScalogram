#' @title Scale index of a signal
#'
#' @description This function computes the scale index of a signal in the scale interval
#' \eqn{[s_0,s_1]}, for a given set of scale parameters \eqn{s_1} and taking \eqn{s_0} as
#' the minimum scale (see Benítez et al. 2010).
#'
#' The scale index of a signal in the scale interval \eqn{[s_0,s_1]} is given by the
#' quotient \deqn{\frac{S(s_{min})}{S(s_{max})},}{S(s_{min})/S(s_{max})} where \eqn{S} is
#' the scalogram, \eqn{s_{max} \in [s_0,s_1]} is the smallest scale such that
#' \eqn{S(s)\le S(s_{max})} for all \eqn{s \in [s_0,s_1]}, and
#' \eqn{s_{min} \in [s_{max},2s_1]} is the smallest scale such that
#' \eqn{S(s_{min})\le S(s)} for all \eqn{s \in [s_{max},2s_1]}.
#'
#' @usage scale_index(signal = NULL,
#'                    scalog = NULL,
#'                    dt = 1,
#'                    scales = NULL,
#'                    powerscales = TRUE,
#'                    s1 = NULL,
#'                    wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
#'                    wparam = NULL,
#'                    waverad = NULL,
#'                    border_effects = c("BE", "INNER", "PER", "SYM"),
#'                    makefigure = TRUE,
#'                    figureperiod = TRUE,
#'                    plot_scalog = FALSE,
#'                    xlab = NULL,
#'                    ylab = "Scale index",
#'                    main = "Scale Index")
#'
#' @param signal A vector containing the signal whose scale indices are wanted.
#' @param scalog A vector containing the scalogram from which the scale indices are going
#' to be computed. If \code{scalog} is not \code{NULL}, then \code{signal}, \code{waverad}
#' and \code{border_effects} are not necessary and they are ignored.
#' @param dt Numeric. The time step of the signals.
#' @param scales A vector containing the wavelet scales at wich the scalogram
#' is computed. This can be either a vector with all the scales or, following Torrence
#' and Compo 1998, a vector of 3 elements with the minimum scale, the maximum scale and
#' the number of suboctaves per octave (in this case, \code{powerscales} must be TRUE in
#' order to construct power 2 scales using a base 2 logarithmic scale). If \code{scales}
#' is NULL, they are automatically constructed.
#' @param powerscales Logical. It must be TRUE (default) in these cases:
#' \itemize{
#' \item If \code{scales} are power 2 scales, i.e. they use a base 2 logarithmic scale.
#' \item If we want to construct power 2 scales automatically. In this case, \code{scales}
#' must be \code{NULL}.
#' \item If we want to construct power 2 scales from \code{scales}. In this case,
#' \code{length(scales)} must be 3.
#' }
#' @param s1 A vector containing the scales \eqn{s_1}. The scale indices are computed in
#' the intervals \eqn{[s_0,s_1]}, where \eqn{s_0} is the minimum scale in \code{scales}.
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
#' @param makefigure Logical. If TRUE (default), a figure with the scale indices is
#' plotted.
#' @param figureperiod Logical. If TRUE (default), periods are used in the figure instead
#' of scales.
#' @param plot_scalog Logical. If TRUE, it plots the scalogram from which the scale indices
#' are computed.
#' @param xlab A string giving a custom X axis label. If NULL (default) the X label is
#' either "s1" or "Period of s1" depending on the value of \code{figureperiod}.
#' @param ylab A string giving a custom Y axis label.
#' @param main A string giving a custom main title for the figure.
#'
#' @return A list with the following fields:
#' \itemize{
#' \item \code{si}: A vector with the scale indices.
#' \item \code{s0}: The scale \eqn{s_0}.
#' \item \code{s1}: A vector with the scales \eqn{s_1}.
#' \item \code{smax}: A vector with the scales \eqn{s_{max}}.
#' \item \code{smin}: A vector with the scales \eqn{s_{min}}.
#' \item \code{scalog}: A vector with the scalogram from which the scale indices are
#' computed.
#' \item \code{scalog_smax}: A vector with the maximum scalogram values \eqn{S(s_{max})}.
#' \item \code{scalog_smin}: A vector with the minimum scalogram values \eqn{S(s_{min})}.
#' \item \code{fourierfactor}: A factor for converting scales into periods.
#' }
#'
#' @examples
#'
#' dt <- 0.1
#' time <- seq(0, 50, dt)
#' signal <- c(sin(pi * time), sin(pi * time / 2))
#' si <- scale_index(signal = signal, dt = dt)
#'
#' # Another way, giving the scalogram instead of the signal:
#'
#' sc <- scalogram(signal = signal, dt = dt, energy_density = FALSE, makefigure = FALSE)
#' si <- scale_index(scalog = sc$scalog, scales = sc$scales, dt = dt)
#'
#' @section References:
#'
#' R. Benítez, V. J. Bolós, M. E. Ramírez. A wavelet-based tool for studying
#' non-periodicity. Comput. Math. Appl. 60 (2010), no. 3, 634-641.
#'
#' @export
#'

scale_index <-
  function(signal = NULL,
           scalog = NULL,
           dt = 1,
           scales = NULL,
           powerscales = TRUE,
           s1 = NULL,
           wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
           wparam = NULL,
           waverad = NULL,
           border_effects = c("BE", "INNER", "PER", "SYM"),
           makefigure = TRUE,
           figureperiod = TRUE,
           plot_scalog = FALSE,
           xlab = NULL,
           ylab = "Scale index",
           main = "Scale Index") {

  wname <- toupper(wname)
  wname <- match.arg(wname)

  fourierfactor <- fourier_factor(wname = wname, wparam = wparam)

  if (!is.null(scales)) {
    if (powerscales && length(scales) == 3) {
      scales <- pow2scales(scales)
    } else {
      if (is.unsorted(scales)) {
        warning("Scales were not sorted.")
        scales <- sort(scales)
      }
    }
    ns <- length(scales)
  }

  if (is.null(scalog)) {

    if (is.null(signal)) {
      stop("We need a signal or a scalogram (sc).")
    }

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

    if (is.null(scales)) {
      scmin <- 2 * dt / fourierfactor
      if (is.null(s1)) {
        scmax <- floor(nt / (2 * waverad)) * dt
      } else {
        scmax <- 2 * max(s1)
      }
      if (powerscales) {
        scales <- pow2scales(c(scmin, scmax, ceiling(256 / log2(scmax / scmin))))
      } else {
        scales <- seq(scmin, scmax, by = (scmax - scmin) / 256)
      }
      ns <- length(scales)
    }

  } else {

    if (is.null(scales)) {
      stop("Scales are needed.")
    }
    if (ns != length(scalog)) {
      stop("The number of scales does not match with length(scalog).")
    }

  }

  if (is.null(s1)) {
    if (scales[ns] < 2 * scales[1]) {
      stop("We need larger scales.")
    }
    index_s1 <- which(scales <= scales[ns] / 2)
    s1 <- scales[index_s1]
    ns1 <- length(s1)
    index_2s1 <- matrix(0, nrow = ns1, ncol = 1)
    for (i in 1:ns1) {
      index_2s1[i] <- max(which(scales <= 2 * s1[i]))
    }
  } else {
    if (is.unsorted(s1)) {
      s1 <- sort(s1)
    }
    ns1 <- length(s1)
    if (s1[1] < scales[1]){
      stop("s1 must be >= the minimum scale.")
    }
    if (2*s1[ns1] > scales[ns]) {
      stop("2*s1 must be <= the maximum scale.")
    }
    index_s1 <- matrix(0, nrow = ns1, ncol = 1)
    index_2s1 <- matrix(0, nrow = ns1, ncol = 1)
    for (i in 1:ns1) {
      index_s1[i] <- max(which(scales <= s1[i]))
      index_2s1[i] <- max(which(scales <= 2 * s1[i]))
    }
  }

  scales <- scales[1:index_2s1[ns1]]
  ns <- length(scales)
  if (ns < 4) {
    stop("We need more scales.")
  }

  if (is.null(scalog)) {

    sc <- scalogram(signal = signal, dt = dt,
                    scales = scales, powerscales = FALSE,
                    wname = wname, wparam = wparam, waverad = waverad,
                    border_effects = border_effects,
                    energy_density = FALSE,
                    makefigure = FALSE)

    scalog <- sc$scalog

  } else {
    scalog <- scalog[1:ns]
  }

  epsilon <- max(scalog) * 1e-6 # This is considered the "numerical zero".

  si <- matrix(0, nrow = ns1, ncol = 1)
  scalog_smin <- si
  scalog_smax <- si
  smax <- si
  smin <- si

  for (i in 1:ns1) {

    # Computation of smax

    scalog_smax[i] <- max(scalog[1:index_s1[i]])
    index_smax <- which.max(scalog[1:index_s1[i]])
    smax[i] <- scales[index_smax]

    # Computation of smin

    scalog_smin[i] <- min(scalog[index_smax:index_2s1[i]])
    index_smin <- which.min(scalog[index_smax:index_2s1[i]])
    smin[i] <- scales[index_smax + index_smin - 1]

    # Computation of SI

    if (scalog_smax[i] < epsilon) {
      si[i] <- 0
    } else {
      si[i] <- scalog_smin[i] / scalog_smax[i]
    }

  }

  if (plot_scalog) {

    if (figureperiod) {
      X <- fourierfactor * scales
      xlab_scalog <- "Period"
    } else {
      X <- scales
      xlab_scalog <- "Scale"
    }

    if (powerscales) {
      plot(log2(X), scalog, type = "l", xlab = xlab_scalog, ylab = "Scalogram", main = "Scalogram", xaxt = "n")
      axis(side = 1, at = floor(log2(X)), labels = 2^(floor(log2(X))))
    } else {
      plot(X, scalog, type = "l", xlab = xlab_scalog, ylab = "Scalogram", main = "Scalogram", xaxt = "n")
      axis(side = 1, at = X[1 + floor((0:8) * (ns - 1) / 8)])
    }

  }

  if (makefigure ) {

    if (ns1 > 1) {
      if (figureperiod) {
        X <- fourierfactor * s1
        if (is.null(xlab)) xlab <- expression('Period of s'[1])
      } else {
        X <- s1
        if (is.null(xlab)) xlab <- expression('s'[1])
      }

      if (powerscales) {
        plot(log2(X), si, type = "l", xlab = xlab, ylab = ylab, main = main, ylim = c(0, 1), xaxt = "n")
        axis(side = 1, at = floor(log2(X)), labels = 2^(floor(log2(X))))
      } else {
        plot(X, si, type = "l", xlab = xlab, ylab = ylab, main = main, ylim = c(0, 1), xaxt = "n")
        axis(side = 1, at = X[1 + floor((0:8) * (ns1 - 1) / 8)])
      }
    } else {
      warning("We can't plot a line with just one point.")
    }

  }

return(list(si = si,
            s0 = scales[1],
            s1 = s1,
            smax = smax,
            smin = smin,
            scalog = scalog,
            scalog_smax = scalog_smax,
            scalog_smin = scalog_smin,
            fourierfactor = fourierfactor)
       )
}
