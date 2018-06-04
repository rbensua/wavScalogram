#' @title Daughter wavelets
#'
#' @description
#' This function computes the daughter wavelets for three different families: "MORLET", "DOG" and "PAUL"
#' (Derivatives of the Gaussian). It is a translation from the Matlab(R) function published by Torrence
#' and Compo following (Torrence & Compo, 1998).
#'
#' This is an internal function used for the computation of the wavelet continuous
#' transform with the \code{cwt_wst} function
#'
#' @usage wave_bases(wname = c("MORLET", "PAUL", "DOG"), k, scale, wparam)
#'
#' @param wname A string, equal to "MORLET", "PAUL" or "DOG".
#' @param k A vector with the Fourier frequencies at which to calculate the wavelet.
#' @param scales A vector with the wavelet scales.
#' @param wparam The corresponding nondimensional parameter for the wavelet function.
#'
#' @return A matrix of size length(scale) x length(k), \code{daughter},
#' containing the fourier transform of the daughter wavelets for the mother "wname",
#' for the given scales (one row for each scale).
#'
#' @examples
#' k <- (1:2000)*2/pi/100
#' scales <- c(1,2)
#' daughter <- wave_bases(wname = "MORLET", k = k, scale = scales)
#' plot(k, daughter[2,], type = "l", xlab = "", ylab = "", main = "Morlet daughter wavelets")
#' lines(k,daughter[1,], col = "red")
#' legend("topright", legend = c("s = 1", "s = 2"), lty = 1, col = c("red", "black"))
#'
#'
#' @section References:
#' C. Torrence, G. P. Compo. A practical guide to wavelet analysis. B. Am. Meteorol. Soc.
#' 79 (1998), 61–78.
#' @export


wave_bases <-
  function(wname = c("MORLET", "PAUL", "DOG"),
           k,
           scales,
           wparam = NULL) {

  wname <- toupper(wname)
  wname <- match.arg(wname)

  n <- length(k)
  ns <- length(scales)

  if (wname == "MORLET") {

    if (is.null(wparam)) {
      wparam <- 6
    }
    k0 <- wparam
    expnt <- -(diag(x = scales, ncol = ns) %*% matrix(rep(k, ns), nrow = ns, byrow = T) - k0) ^ 2 / 2
    expnt <- sweep(expnt, MARGIN = 2, (k > 0), `*`)
    Norm <- sqrt(scales * k[2])  *pi ^ (-.25) * sqrt(n)
    daughter <- diag(x = Norm, ncol = ns) %*% exp(expnt)
    daughter <- sweep(daughter, MARGIN = 2, (k > 0), `*`)

  #  fourier_factor <- (4 * pi) / (k0 + sqrt(2 + k0 ^ 2))
  #  coi <- fourier_factor / sqrt(2)

  } else if (wname == "PAUL") {

    if (is.null(wparam)) {
      wparam <- 4
    }
    m <- wparam
    expnt <- -diag(x = scales, ncol = ns) %*% matrix(rep(k, ns), nrow = ns, byrow = T)
    expnt <- sweep(expnt, MARGIN = 2, (k > 0), `*`)
    Norm <- sqrt(scales * k[2]) * (2 ^ m / sqrt(m * prod(2:(2 * m - 1)))) * sqrt(n)
    daughter <- diag(x = Norm, ncol = ns) %*% expnt ^ m * exp(expnt)
    daughter <- sweep(daughter, MARGIN = 2, (k > 0), `*`)

  #  fourier_factor <- 4 * pi / (2 * m + 1)
  #  coi  <- fourier_factor * sqrt(2)

  } else if (wname == "DOG") {

    if (is.null(wparam)) {
      wparam <- 2
    }
    m <- wparam
    preexpnt <- (diag(x = scales, ncol = ns) %*% matrix(rep(k, ns), nrow = ns, byrow = T))
    expnt <- -preexpnt ^ 2 / 2
    Norm <- sqrt(scales * k[2] / gamma(m + 0.5)) * sqrt(n)
    daughter <- diag(x = -Norm * (1i ^ m), ncol = ns) %*% preexpnt ^ m * exp(expnt)

  #  fourier_factor <- 2 * pi * sqrt(2 / (2 * m + 1))
  #  coi <- fourier_factor / sqrt(2)
  }
#  return(list(daughter = daughter, fourier_factor = fourier_factor, coi = coi))
  return(daughter = daughter)

}
