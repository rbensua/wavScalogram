#' @title Scalogram of a signal
#'
#' @description This function computes the normalized scalogram of a signal for the scales
#' given. It is important to note that the notion of scalogram here is analogous
#' to the spectrum of the Fourier transform. It gives the contribution of each scale to
#' the total energy of the signal. For each scale \eqn{s}, it is defined as the square
#' root of the integral of the squared modulus of the wavelet transform w.r.t. the time
#' variable \eqn{t}, i.e.
#'
#'   \deqn{S(s):= (\int_{-\infty}^{+\infty}|Wf(t,s)|^2 dt)^{1/2}.}{%
#'   S(s):=(integral_{-\infty}^{+\infty}|Wf(t,s)|^2 dt)^{1/2}.}
#'
#' "Normalized" means that the scalogram is divided by the square root of the number of
#' times, for comparison purposes between different values of the parameter
#' \code{border_effects}.
#'
#' @usage scalogram(signal,
#'                  dt = 1,
#'                  scales = NULL,
#'                  powerscales = TRUE,
#'                  wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
#'                  wparam = NULL,
#'                  waverad = NULL,
#'                  border_effects = c("BE", "INNER", "PER", "SYM"),
#'                  energy_density = TRUE,
#'                  makefigure = TRUE,
#'                  figureperiod = TRUE,
#'                  xlab = NULL,
#'                  ylab = "Scalogram",
#'                  main = "Scalogram")
#'
#' @param signal A vector containing the signal whose scalogram is wanted.
#' @param dt Numeric. The time step of the signal.
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
#' @param wname A string, equal to "MORLET", "DOG", "PAUL", "HAAR" or "HAAR2". The
#' difference between "HAAR" and "HAAR2" is that "HAAR2" is more accurate but slower.
#' @param wparam The corresponding nondimensional parameter for the wavelet function
#' (Morlet, DoG or Paul).
#' @param waverad Numeric. The radius of the wavelet used in the computations for the cone
#' of influence. If it is not specified, it is asumed to be \eqn{\sqrt{2}} for Morlet and DoG,
#' \eqn{1/\sqrt{2}} for Paul and 0.5 for Haar.
#' @param border_effects String, equal to "BE", "INNER", "PER" or "SYM", which indicates
#' how to manage the border effects which arise usually when a convolution is performed on
#' finite-lenght signals.
#' \itemize{
#' \item "BE": With border effects, padding time series with zeroes.
#' \item "INNER": Normalized inner scalogram with security margin adapted for each
#'     different scale.
#' \item "PER": With border effects, using boundary wavelets (periodization of the
#' original time series).
#' \item "SYM": With border effects, using a symmetric catenation of the original time
#' series.
#' }
#' @param energy_density Logical. If TRUE (default), divide the scalogram by the square
#' root of the scales for convert it into energy density.
#' @param makefigure Logical. If TRUE (default), a figure with the scalogram is plotted.
#' @param figureperiod Logical. If TRUE (default), periods are used in the figure instead
#' of scales.
#' @param xlab A string giving a custom X axis label. If NULL (default) the X label is
#' either "Scale" or "Period" depending on the value of \code{figureperiod}.
#' @param ylab A string giving a custom Y axis label.
#' @param main A string giving a custom main title for the figure.
#'
#' @return A list with the following fields:
#' \itemize{
#' \item \code{scalog}: A vector of length \code{length(scales)}, containing the values of
#' the scalogram at each scale.
#' \item \code{scales}: The vector of scales.
#' \item \code{energy}: If \code{energy_density} is TRUE, it is the \eqn{L^2} norm of
#' \code{scalog}.
#' \item \code{fourierfactor}: A factor for converting scales into periods.
#' }
#'
#' @examples
#' dt <- 0.1
#' time <- seq(0, 50, dt)
#' signal <- c(sin(pi * time), sin(pi * time / 2))
#' scalog <- scalogram(signal = signal, dt = dt, border_effects = "INNER")
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
           dt = 1,
           scales = NULL,
           powerscales = TRUE,
           wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
           wparam = NULL,
           waverad = NULL,
           border_effects = c("BE", "INNER", "PER", "SYM"),
           energy_density = TRUE,
           makefigure = TRUE,
           figureperiod = TRUE,
           xlab = NULL,
           ylab = "Scalogram",
           main = "Scalogram") {

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
  if (border_effects == "INNER") {
    border_effects_cwt <- "BE"
  } else {
    border_effects_cwt <- border_effects
  }

  nt <- length(signal)

  cwt <- cwt_wst(signal = signal,
                   dt = dt,
                   scales = scales,
                   powerscales = powerscales,
                   wname = wname,
                   wparam = wparam,
                   waverad = waverad,
                   border_effects = border_effects_cwt,
                   makefigure = FALSE)

  scalesdt <- cwt$scales
  scales <- scalesdt / dt
  ns <- length(scales)
  coefs <- matrix(cwt$coefs, nrow = nt, ncol = ns)

  scalog <- numeric(ns)

  if (border_effects == "INNER") {

    for (i in 1:ns) {
      nt_margin <- ceiling(waverad * scales[i])
      nt_end <- nt - nt_margin
      nt_ini <- nt_margin
      if (nt_end < nt_ini){
        print("We need more times for these scales")
        scalog <- scalog[1:(i - 1)]
        scales <- scales[1:(i - 1)]
        scalesdt <- scalesdt[1:(i - 1)]
        ns <- i - 1
        break
      }
      scalog[i] <- sqrt(sum(abs(coefs[nt_ini:nt_end, i]) ^ 2))
      scalog[i] <- scalog[i] / sqrt(nt_end - nt_ini + 1) # Normalization
    }

  } else {

    scalog <- sqrt(colSums(abs(coefs) ^ 2))
    scalog <- scalog / sqrt(nt) # Normalization

  }

  energy <- NA
  if (energy_density) {
    scalog <- scalog / sqrt(scalesdt)
    energy <- sqrt(sum(scalog ^ 2))
  }

  fourierfactor <- cwt$fourierfactor

  if (makefigure) {

    if (figureperiod) {
      X <- fourierfactor * scalesdt
      if (is.null(xlab)) xlab <- "Period"
    } else {
      X <- scalesdt
      if (is.null(xlab)) xlab <- "Scale"
    }

    if (powerscales) {
      plot(log2(X), scalog, type = "l", xlab = xlab, ylab = ylab, main = main, xaxt = "n")
      axis(side = 1, at = floor(log2(X)), labels = 2^(floor(log2(X))))
    } else {
      plot(X, scalog, type = "l", xlab = xlab, ylab = ylab, main = main, xaxt = "n")
      axis(side = 1, at = X[1 + floor((0:8) * (ns - 1) / 8)])
    }

  }

return(list(
  scalog = scalog,
  scales = scalesdt,
  energy = energy,
  fourierfactor = fourierfactor
  ))

}
