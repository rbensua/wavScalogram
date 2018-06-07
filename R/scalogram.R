#' @title Scalogram of a signal
#'
#' @description This function computes the normalized scalogram of a signal for the scales
#' given. It is important to note that the notion of scalogram here is analogous
#' to the spectrum of the Fourier transform. It gives the contribution of each scale to
#' the total energy of the signal. For each scale \eqn{s}, it is defined as the square
#' root of the integral of the squared modulus of the wavelet transform w.r.t the time
#' variable \eqn{t}, i.e.
#'
#'   \deqn{S(s):= (\int_{-\infty}^{+\infty}|Wf(t,s)|^2 dt)^{1/2}.}{%
#'   S(s):=(integral_{-\infty}^{+\infty}|Wf(t,s)|^2 dt)^{1/2}.}
#'
#' "Normalized" means that the scalogram is divided by the square root of the number of
#' times, for comparison purposes between different values of the parameter
#' \code{border_effects}.
#'
#'
#' @usage scalogram(signal,
#'                  scales = NULL,
#'                  powerscales = TRUE,
#'                  wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
#'                  wparam = NULL,
#'                  border_effects = c("BE", "INNER", "PER", "SYM"),
#'                  energy_density = TRUE,
#'                  makefigure = FALSE)
#'
#' @param signal A vector containing the signal whose wavelet transform is wanted. The
#' unit of time is taken as the time difference between two consecutive data.
#' @param scales A vector with the wavelet scales measured in units of time. This can be
#' either a vector with all the scales, or (if \code{powerscales = TRUE}) following
#' Torrence and Compo 1998, a vector of three elements with the minimum scale, the maximum
#' scale and the number of suboctaves per octave.
#' @param powerscales Logical. Construct power 2 scales.
#' @param wname A string, equal to "MORLET", "DOG", "PAUL", "HAAR" or "HAAR2". The
#' difference between "HAAR" and "HAAR2" is that "HAAR2" is more accurate but slower.
#' @param wparam Numeric. Parameters of the corresponding wavelet.
#' @param border_effects String, equal to "BE", "INNER", "PER" or "SYM",
#' which indicates how to manage the border effects which arise usually when a convolution
#' is performed on finite-lenght signals.
#' \itemize{
#' \item "BE": With border effects, padding time series with zeroes.
#' \item "INNER": Normalized inner scalogram with security margin adapted for each
#'     different scale.
#' \item "PER": With border effects, using boundary wavelets (periodization of the
#' original time series).
#' \item "SYM": With border effects, using a symmetric catenation of the original time
#' series.
#' }
#' @param energy_density Logical. Divide the scalogram by sqrt(scales) for convert it into
#' energy density (defaults to TRUE).
#' @param makefigure Logical. Plots a figure with the scalogram.
#'
#' @return A list with the following fields:
#'
#' \code{scalog}: A vector of length \code{length(scales)}, containing the values of the
#' scalogram at each scale.
#'
#' \code{scales}: The vector of scales.
#'
#' @examples
#' time <- 1:500
#' signal <- c(sin(pi * time / 8), sin(pi * time / 16))
#' scalog <- scalogram(signal = signal, border_effects = "INNER", makefigure = TRUE)
#'
#' @section References:
#'
#' C. Torrence, G. P. Compo. A practical guide to wavelet analysis. B. Am. Meteorol. Soc.
#' 79 (1998), 61–78.
#'
#' V. J. Bolós, R. Benítez, R. Ferrer, R. Jammazi. The windowed scalogram difference: a
#' novel wavelet tool for comparing time series. Appl. Math. Comput., 312 (2017), 49-65.
#'
#' @export
#'

scalogram <-
  function(signal,
           scales = NULL,
           powerscales = TRUE,
           wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
           wparam = NULL,
           border_effects = c("BE", "INNER", "PER", "SYM"),
           energy_density = TRUE,
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

  if (is.null(scales)) {
    scmax <- floor(nt / (2 * waverad))
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

  scalog <- numeric(ns)

  if (border_effects == "BE") {

    coefs <- cwt_wst(scales = scales, signal = signal, wname = wname, wparam = wparam)
    scalog <- sqrt(colSums(abs(coefs) ^ 2))
    scalog <- scalog / sqrt(nt) # Normalization

  } else if(border_effects == "INNER") {

    coefs <- cwt_wst(scales = scales, signal = signal, wname = wname, wparam = wparam)
    scalog <- numeric(ns)
    for (i in 1:ns){
      nt_margin <- ceiling(waverad * scales[i])
      nt_end <- nt - nt_margin
      nt_ini <- nt_margin
      if (nt_end < nt_ini){
        print("We need more times for these scales")
        scalog <- scalog[1:(i - 1)]
        scales <- scales[1:(i - 1)]
        ns <- i-1
        break
      }
      scalog[i] <- sqrt(sum(abs(coefs[nt_ini:nt_end, i]) ^ 2))
      scalog[i] <- scalog[i] / sqrt(nt_end - nt_ini + 1) # Normalization
    }

  } else if(border_effects == "PER") {

    ndatcat <- ceiling(ceiling(waverad * scales[ns]) / nt) # Number of catenations = ndatcat*2+1
    ntcat <- (ndatcat * 2 + 1) * nt
    datcat <- numeric(ntcat)
    for (itcat in 1:ntcat) {
      datcat[itcat] = signal[itcat - floor((itcat - 1) / nt) * nt]
    }
    coefs <- cwt_wst(scales = scales, signal = datcat, wname = wname, wparam = wparam)
    nt_ini <- ndatcat * nt + 1
    nt_end <- nt_ini + nt - 1
    for (i in 1:ns) {
      scalog[i] <- sqrt(sum(abs(coefs[nt_ini:nt_end, i]) ^ 2))
    }
    scalog <- scalog / sqrt(nt) # Normalization

  } else {

    ndatcat <- ceiling(ceiling(waverad * scales[ns]) / nt) # Number of catenations = ndatcat*2+1
    dat_sym <- signal
    aux <- rev(signal)
    for (icat in 1:ndatcat) {
      dat_sym <- c(aux, dat_sym, aux)
      aux <- rev(aux)
    }
    coefs <- cwt_wst(scales = scales, signal = dat_sym, wname = wname, wparam = wparam)
    nt_ini <- ndatcat * nt + 1
    nt_end <- nt_ini + nt - 1
    for (i in 1:ns) {
      scalog[i] <- sqrt(sum(abs(coefs[nt_ini:nt_end, i]) ^ 2))
    }
    scalog <- scalog / sqrt(nt) # Normalization
  }

  if (energy_density) {
    scalog <- scalog / sqrt(scales)
  }

  if (makefigure) {
    if (powerscales) {
      plot(log2(scales), scalog, type = "l", xlab = "Scale", main = "Scalogram", xaxt = "n")
      axis(side = 1, at = floor(log2(scales)), labels = 2^(floor(log2(scales))))
    } else {
      plot(scales, scalog, type = "l", xlab = "Scale", main = "Scalogram", xaxt = "n")
      axis(side = 1, at = scales[1 + floor((0:8) * (ns - 1) / 8)])
    }
  }

return(list(
  scalog = scalog,
  scales = scales
  ))
}
