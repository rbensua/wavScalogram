#' @title Scale index of a signal
#'
#' @description This function computes the scale index of a signal in the scale interval
#' \eqn{[s_0,s_1]}, for a given set of scale parameters \eqn{s_1} and taking \eqn{s_0} as
#' the minimum scale (see Benítez et al. 2010).
#'
#' The scale index of a signal in the scale interval \eqn{[s_0,s_1]} is given by the
#' quotient \deqn{\frac{S(s_{min})}{S(s_{max})},}{S(s_{min})/S(s_{max})} where \eqn{S} is
#' the scalogram, \eqn{s_{max} in [s_0,s_1]} is the smallest scale such that
#' \eqn{S(s)\le S(s_{max})} for all \eqn{s in [s_0,s_1]}, and
#' \eqn{s_{min} in [s_{max},2s_1]} is the smallest scale such that
#' \eqn{S(s_{min})\le S(s)} for all \eqn{s in [s_{max},2s_1]}.
#'
#' @usage scale_index(signal,
#'                    scales = NULL,
#'                    powerscales = TRUE,
#'                    s1 = NULL,
#'                    wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
#'                    wparam = NULL,
#'                    border_effects = c("BE", "INNER", "PER", "SYM"),
#'                    makefigure = FALSE)
#'
#' @param signal A vector containing the signal whose scale indices are wanted. The unit
#' of time is taken as the time difference between two consecutive data.
#' @param scales A vector containing the wavelet scales measured in units of time. This
#' can be either a vector with all the scales, or (if \code{powerscales = TRUE}) following
#' Torrence and Compo 1998, a vector of three elements with the minimum scale, the maximum
#' scale and the number of suboctaves per octave.
#' @param powerscales Logical. Construct power 2 scales.
#' @param s1 A vector containing the scales \eqn{s_1}. The scale indices are computed in
#' the intervals \eqn{[s_0,s_1]}, where \eqn{s_0} is the minimum scale.
#' @param wname A string, equal to "MORLET", "DOG", "PAUL", "HAAR" or "HAAR2". The
#' difference between "HAAR" and "HAAR2" is that "HAAR2" is more accurate but slower.
#' @param wparam Numeric. Parameters of the corresponding wavelet.
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
#' @param makefigure Logical. Plots a figure with the scale indices.
#'
#' @return A list with the following fields:
#'
#' \code{si}: A vector with the scale indices.
#'
#' \code{smax}: A vector with the scales \eqn{s_{max}}.
#'
#' \code{smin}: A vector with the scales \eqn{s_{min}}.
#'
#' @examples
#' time <- 1:500
#' signal <- c(sin(pi * time / 8), sin(pi * time / 16))
#' si <- scale_index(signal = signal, makefigure = TRUE)
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
           scales = NULL,
           powerscales = TRUE,
           s1 = NULL,
           wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
           wparam = NULL,
           border_effects = c("BE", "INNER", "PER", "SYM"),
           makefigure = FALSE) {

  wname <- toupper(wname)
  wname <- match.arg(wname)
  if (wname == "MORLET") {
    waverad <- 3
  } else if ((wname == "HAAR") || (wname == "HAAR2")) {
    waverad <- 0.5
  } else {
    waverad <- wavelet_radius(wname = wname, wparam = wparam)
    waverad <- waverad$left
  }

  border_effects <- toupper(border_effects)
  border_effects <- match.arg(border_effects)

  nt <- length(signal)

  if (is.unsorted(s1)) {
    s1 <- sort(s1)
  }

  if (is.null(scales)) {
    if (is.null(s1)) {
      scmax <- floor(nt / (2 * waverad))
    } else {
      scmax <- 2 * s1[length(s1)]
    }
    if (powerscales) {
      scales <- pow2scales(c(2, scmax, ceiling(256 / (log2(scmax) - 1))))
    } else {
      scales <- seq(2, scmax, by = scmax / 256)
    }
  } else {
    if (powerscales) {
      scales <- pow2scales(scales)
    } else {
      if (is.unsorted(scales))
        scales <- sort(scales)
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

  scalog <- scalogram(signal = signal, scales = scales, powerscales = FALSE, wname = wname, wparam = wparam,
                      border_effects = border_effects, energy_density = FALSE)

  scalog <- scalog$scalog
  epsilon <- max(scalog) * 1e-6 # This is considered the "numerical zero".

  smin <- matrix(0, nrow = ns1, ncol = 1)
  smax <- matrix(0, nrow = ns1, ncol = 1)
  si <- matrix(0, nrow = ns1, ncol = 1)
  for (i in 1:ns1) {

    # Computation of smax

    smax[i] <- max(scalog[1:index_s1[i]])
    index_smax <- which.max(scalog[1:index_s1[i]])

    # Computation of smin

    smin[i] <- min(scalog[index_smax:index_2s1[i]])

    # Computation of SI

    if (smax[i] < epsilon) {
      si[i] <- 0
    } else {
      si[i] <- smin[i] / smax[i]
    }

  }

  if (makefigure ) {
    if (length(s1) > 1){
    plot(s1, si, type = "l")
    } else{
      stop("We can't plot a line with a unique point.")
    }
  }

return(list(si = si,
            smax = smax,
            smin = smin)
       )
}
