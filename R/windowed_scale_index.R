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
#' radius \code{windowrad}, \eqn{s_{max} in [s_0,s_1]} is the smallest scale such that
#' \eqn{WS_{windowrad}(tc,s)\le WS_{windowrad}(tc,s_{max})} for all \eqn{s in [s_0,s_1]},
#' and \eqn{s_{min} in [s_{max},2s_1]} is the smallest scale such that
#' \eqn{WS_{windowrad}(tc,s_{min})\le WS_{windowrad}(tc,s)} for all
#' \eqn{s in [s_{max},2s_1]}.
#'
#' @usage windowed_scale_index(signal,
#'                             scales = NULL,
#'                             powerscales = TRUE,
#'                             s1 = NULL,
#'                             windowrad = NULL,
#'                             delta_t = NULL,
#'                             wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
#'                             wparam = NULL,
#'                             border_effects = c("BE", "INNER", "PER", "SYM"),
#'                             makefigure = FALSE)
#'
#' @param signal A vector containing the signal whose windowed scale indices are wanted.
#' The unit of time is taken as the time difference between two consecutive data.
#' @param scales A vector containing the wavelet scales measured in units of time. This
#' can be either a vector with all the scales, or (if \code{powerscales = TRUE})
#' following Torrence and Compo 1998, a vector of three elements with the minimum scale,
#' the maximum scale and the number of suboctaves per octave.
#' @param powerscales Logical. Construct power 2 scales.
#' @param s1 A vector containing the scales \eqn{s_1}. The windowed scale indices are
#' computed in the intervals \eqn{[s_0,s_1]}, where \eqn{s_0} is the minimum scale.
#' @param windowrad Numeric. Radius for the time windows.
#' @param delta_t Numeric. Increment of time for the construction of windows central
#' times.
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
#' @param makefigure Logical. Plots a figure with the windowed scale indices.
#'
#' @return A list with the following fields:
#'
#' \code{wsi}: A matrix of size \code{length(tcentral)}x\code{length(s1)} containing
#' the values of the corresponding windowed scale indices.
#'
#' \code{smax}: A vector of length \code{length(tcentral)} containing the scales
#' \eqn{s_{max}}.
#'
#' \code{smin}: A vector of length \code{length(tcentral)} containing the scales
#' \eqn{s_{min}}.
#'
#' \code{tcentral}: The vector of central times used in the computation of the \code{wsi}.
#'
#' \code{s1}: The vector of scales \eqn{s_1}.
#'
#' \code{coi_maxscale}: A vector of length \code{length(tcentral)} containing the values
#' of the maximum scale from which there are border effects.
#'
#' @importFrom fields image.plot
#' @importFrom colorRamps matlab.like
#'
#' @examples
#' time <- 1:500
#' signal <- c(sin(pi * time / 8), sin(pi * time / 16))
#' wsi <- windowed_scale_index(signal = signal, makefigure = TRUE)
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
           scales = NULL,
           powerscales = TRUE,
           s1 = NULL,
           windowrad = NULL,
           delta_t = NULL,
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

  if (is.null(scales)) {
    if (is.null(s1)) {
      scmax <- floor((nt - 2 * windowrad) / (2 * waverad))
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
    windowed_scalogram(signal = signal, scales = scales, powerscales = FALSE, windowrad = windowrad,
                       delta_t = delta_t, wname = wname, wparam = wparam, border_effects = border_effects,
                       energy_density = FALSE, makefigure = FALSE)

  nwsi <- length(wsc$tcentral)

  wsi <- matrix(NA, nrow = nwsi, ncol = ns1)
  smax <- matrix(0, nrow = nwsi, ncol = ns1)
  smin <- matrix(0, nrow = nwsi, ncol = ns1)
  epsilon <- max(wsc$wsc, na.rm = TRUE) * 1e-6 # This is considered the "numerical zero"
  for (i in 1:nwsi) {
    ni <- max(which(!is.na(wsc$wsc[i, index_2s1]))) # If border_effects = "INNER" there are NA in wsc.
    for (j in 1:ni) {
      smax[i, j] <- max(wsc$wsc[i, 1:index_s1[j]])
      if (smax[i, j] > epsilon) {
        index_smax <- which.max(wsc$wsc[i, 1:index_s1[j]])
        smin[i,j] <- min(wsc$wsc[i, index_smax:index_2s1[j]])
        wsi[i,j] <- smin[i, j] / smax[i, j]
      }
    }
  }

  # COI
  coi_maxscale <- wsc$coi_maxscale / 2

  if (makefigure) {
    wavPlot(
      Z = wsi,
      X = wsc$tcentral,
      Y = s1,
      Ylog = powerscales,
      coi = coi_maxscale,
      Xname = "Time",
      Yname = "s1",
      Zname = "WSI"
    )
  }


  return(list(
    wsi = wsi,
    smax = smax,
    smin = smin,
    tcentral = wsc$tcentral,
    s1 = s1,
    coi_maxscale = coi_maxscale
  ))
}
