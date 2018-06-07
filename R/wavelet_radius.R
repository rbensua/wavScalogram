#' @title Wavelet radius
#'
#' @description
#' This function computes an approximation of the effective radius of a mother wavelet.
#'
#' @usage wavelet_radius(wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
#'                       wparam = NULL,
#'                       perc = .0025,
#'                       scale = 100,
#'                       n = 1000,
#'                       makefigure = FALSE)
#'
#' @param wname A string, equal to "MORLET", "DOG", "PAUL", "HAAR" or "HAAR2". The
#' difference between "HAAR" and "HAAR2" is that "HAAR2" is more accurate but slower.
#' @param wparam Numeric. Parameters of the corresponding wavelet.
#' @param perc Numeric. The wavelet radius is computed so that the area covered is at
#' least the 100*(1-\code{perc})\% of the total area of the mother wavelet.
#' @param scale Numeric. Scale of the wavelet used in the computations. It only affects
#' the accuracy.
#' @param n Numeric. The computations use a time series of length \eqn{2n+1}.
#' @param makefigure Logical. Plots a figure with the real part of the mother wavelet and
#' its modulus.
#'
#' @return A list with the following fields:
#'
#' \code{left}: The radius on the left.
#'
#' \code{right}: The radius on the right.
#'
#' @examples
#' waverad <- wavelet_radius(wname = "MORLET", makefigure = TRUE)
#'
#' @export
#'

wavelet_radius <-
  function(wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
           wparam = NULL,
           perc = .0025,
           scale = 100,
           n = 1000,
           makefigure = FALSE) {

  wname <- toupper(wname)
  wname <- match.arg(wname)

  signal <- c(numeric(n), 1, numeric(n))
  nt <- 2 * n + 1

  coefs <- cwt_wst(scales = scale, signal = signal, wname = wname, wparam = wparam)
  abscoefs <- abs(coefs)

  area <- numeric(nt)
  area[1] <- abscoefs[1]
  for (i in 2:nt) {
    area[i] <- area[i - 1] + abscoefs[i]
  }

  kk <- perc * area[nt] / 2

  i_right <- max(which(area <= kk))
  right <- (n + 1 - i_right - (kk - area[i_right]) / abscoefs[i_right + 1]) / scale

  abscoefs_rev <- rev(abscoefs)
  area_rev <- numeric(nt)
  area_rev[1] <- abscoefs_rev[1]
  for (i in 2:nt) {
    area_rev[i] <- area_rev[i - 1] + abscoefs_rev[i]
  }

  i_left <- max(which(area_rev <= kk))
  left <- (n + 1 - i_left - (kk - area_rev[i_left]) / abscoefs_rev[i_left + 1]) / scale
  i_left <- 2 * n + 2 - i_left

  if (makefigure) {
    i_left_margin <- min(nt, i_left + floor((i_left - n - 1) / 4))
    i_right_margin <- max(1, i_right - floor((n + 1 - i_right) / 4))
    xleft <- - (i_left_margin - n - 1) / scale
    xright <- (n + 1 - i_right_margin) / scale
    x <- seq(from = xleft, to = xright, by = 1 / scale)
    ylim <- sqrt(scale) * range(abscoefs, Re(coefs))
    plot(x, sqrt(scale) * abscoefs[i_left_margin:i_right_margin], type = "l", ylim = ylim, lty = 3)
    lines(x, sqrt(scale) * Re(coefs[i_left_margin:i_right_margin]))
    segments(x0 = -right, x1 = -right, y0 = ylim[1], y1 = ylim[2])
    segments(x0 = left, x1 = left, y0 = ylim[1], y1 = ylim[2])
  }

  return(list(
    left = left,
    right = right
    ))
}
