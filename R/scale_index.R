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
#' @usage scale_index(signal,
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
#'                    xlab = NULL,
#'                    ylab = "Scale index",
#'                    main = "Scale Index")
#'
#' @param signal A vector containing the signal whose scale indices are wanted.
#' @param dt Numeric. The time step of the signals.
#' @param scales A vector containing the wavelet scales at wich
#' the scalogram is computed. This can be either a vector with all the scales, or
#' (if \code{powerscales} is TRUE) following Torrence and Compo 1998, a vector of three
#' elements with the minimum scale, the maximum scale and the number of suboctaves per
#' octave. If NULL, they are automatically computed.
#' @param powerscales Logical. If TRUE (default), construct power 2 scales from
#' \code{scales}. If \code{scales} is NULL, they are automatically computed.
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
#' @param xlab A string giving a custom X axis label. If NULL (default) the X label is
#' either "s1" or "Period of s1" depending on the value of \code{figureperiod}.
#' @param ylab A string giving a custom Y axis label.
#' @param main A string giving a custom main title for the figure.
#'
#' @return A list with the following fields:
#' \itemize{
#' \item \code{si}: A vector with the scale indices.
#' \item \code{s1}: A vector containing the scales \eqn{s_1}.
#' \item \code{smax}: A vector with the scales \eqn{s_{max}}.
#' \item \code{smin}: A vector with the scales \eqn{s_{min}}.
#' \item \code{scalog_smax}: A vector with the maximum scalogram values \eqn{S(s_{max})}.
#' \item \code{scalog_smin}: A vector with the minimum scalogram values \eqn{S(s_{min})}.
#' \item \code{fourierfactor}: A factor for converting scales into periods.
#' }
#'
#' @examples
#' dt <- 0.1
#' time <- seq(0, 50, dt)
#' signal <- c(sin(pi * time), sin(pi * time / 2))
#' si <- scale_index(signal = signal, dt = dt)
#'
#' @section References:
#'
#' R. Benítez, V. J. Bolós, M. E. Ramírez. A wavelet-based tool for studying
#' non-periodicity. Comput. Math. Appl. 60 (2010), no. 3, 634-641.
#'
#' @export
#'

scale_index <-
  function(signal,
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
           xlab = NULL,
           ylab = "Scale index",
           main = "Scale Index") {

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

  if (is.unsorted(s1)) {
    s1 <- sort(s1)
  }

  fourierfactor <- fourier_factor(wname = wname, wparam = wparam)

  if (is.null(scales)) {
    scmin <- 2 * dt / fourierfactor
    if (is.null(s1)) {
      scmax <- floor(nt / (2 * waverad)) * dt
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
    index_2s1 <- matrix(0, nrow = ns1, ncol = 1)
    for (i in 1:ns1) {
      index_2s1[i] <- max(which(scales <= 2 * s1[i]))
    }
  } else {
    s1null <- FALSE
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

  sc <- scalogram(signal = signal, dt = dt,
                  scales = scales, powerscales = FALSE,
                  wname = wname, wparam = wparam, waverad = waverad,
                  border_effects = border_effects,
                  energy_density = FALSE,
                  makefigure = FALSE)

  scalog <- sc$scalog
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

  if (makefigure ) {
    if (length(s1) > 1) {
      if (figureperiod) {
        X <- fourierfactor * s1
        if (is.null(xlab)) xlab <- expression('Period of s'[1])
      } else {
        X <- s1
        if (is.null(xlab)) xlab <- expression('s'[1])
      }
      plot(X, si, type = "l", xlab = xlab, ylab = ylab, main = main)
    } else {
      warning("We can't plot a line with just one point.")
    }
  }

return(list(si = si,
            s1 = s1,
            smax = smax,
            smin = smin,
            scalog_smax = scalog_smax,
            scalog_smin = scalog_smin,
            fourierfactor = fourierfactor)
       )
}
