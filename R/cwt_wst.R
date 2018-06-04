#' @title Continuous wavelet transform
#'
#' @description
#' This function computes the continuous wavelet transform for some families of wavelet bases:
#' "MORLET", "DOG", "PAUL" and "HAAR".
#' It is a translation from the Matlab(R) function published by Torrence
#' and Compo (Torrence & Compo, 1998).
#'
#'
#' @usage cwt_wst(scales, signal, wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"), wparam)
#' @param scales A vector with the wavelet scales.
#' @param signal A vector containing the sigal whose wavelet transform is wanted.
#' @param wname A string, equal to "MORLET", "DOG", "PAUL", "HAAR" or "HAAR2". The difference between "HAAR" and "HAAR2" is that "HAAR2" is more accurate but slower.
#' @param wparam The corresponding nondimensional parameter for the wavelet function.
#'
#' @return A matrix of size length(signal) x length(scales),
#' containing the CWT coefficients of the signal.
#'
#' @examples
#' t <- 0:1000
#' nt <- length(t)
#' signal <- sin(2 * pi * t / 100) * (t <= 500) + 2 * cos(20 * pi * t / 100) * (t > 500)
#' plot(t, signal, type = "l")
#' dt <- 1
#' # Defining the scales (Torrence and Compo's way)
#' s0 <- 2 * dt
#' Dj <- 8
#' waverad <- 3 # Morlet wavelet radius
#' J <- log2((t[nt] - t[1] - 2) / (2 * s0 * dt))
#' scales <- s0 * 2 ^ ((0:(J * Dj)) / Dj)
#' cwt <- cwt_wst(scales = scales, signal = signal, wname = "MORLET")
#' \dontrun{
#' par(mfrow = c(2, 1))
#' plot(t, signal, type = "l")
#' image(t, log2(scales), abs(t(cwt)), ylim = rev(range(log2(scales))))
#'}
#'
#' @section References:
#'
#' C. Torrence, G. P. Compo. A practical guide to wavelet analysis. B. Am. Meteorol. Soc.
#' 79 (1998), 61–78.
#'
#' @export
#'

cwt_wst <-
  function(scales,
           signal,
           wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
           wparam = NULL) {

  wname <- toupper(wname)
  wname <- match.arg(wname)
  precision <- max(12, floor(log2(max(scales))) + 4)
  nt <- length(signal)
  ns <- length(scales)
  coefs <- matrix(0, nrow = nt, ncol = ns)

  if (wname == "HAAR") {

    step <- 2 ^ (-precision)
    xval <- seq(from = 0, to = 1 + step, by = step)
    xmax <- max(xval)
    psi_integ <- xval * (xval < 1 / 2) + (1 - xval) * (xval >= 1 / 2)
    for (k in 1:ns) {
      a <- scales[k]
      j <- 1 + floor((0:a * xmax) / (a * step))
      if (length(j) == 1) {
        j <- c(1, 1)
      }
      f <- rev(psi_integ[j])
      coefs[, k] <- -sqrt(a) * core(diff(convolve(signal, f, type = "open")), nt)
    }

  } else if (wname == "HAAR2") {

    for (j in 1:ns) {
      kmatrix <-  matrix(0, nrow = nt, ncol = nt)
      klen <- floor((scales[j] - 1) / 2)
      kfrac <- (((scales[j] - 1) / 2) - klen)
      for (i in 1:nt) {
        kmin <- max(1, i - klen)
        kmax <- min(nt, i + klen)
        if (kmin <= (i - 1)) {
          kmatrix[kmin:(i - 1), i] <- 1
          #for (k in kmin:(i - 1)) {
          #  coefs[i, j] = coefs[i, j] + signal[k]
          #}
        }
        if (kmin >= 2) {
          kmatrix[kmin - 1, i] <- kfrac
          #coefs[i, j] <- coefs[i, j] + signal[kmin - 1] * kfrac
        }
        if (kmax >= (i + 1)) {
          kmatrix[(i + 1):kmax, i] <- -1
          #for (k in (i + 1):kmax) {
          #  coefs[i, j] = coefs[i, j] - signal[k]
          #}
        }
        if (kmax <= (nt - 1)) {
          kmatrix[kmax + 1, i] <- -kfrac
          #coefs[i, j] <- coefs[i, j] - signal[kmax + 1] * kfrac
        }
      }
      coefs[, j] <- signal %*% kmatrix / sqrt(scales[j])
      #coefs[, j] <- coefs[, j] / sqrt(scales[j])
    }

  } else {

    coefs <- coefs + 1i*coefs
    f <- fft(signal)
    k <- 1:trunc(nt / 2)
    k <- k * 2 * pi / nt
    k <- c(0, k, -k[trunc((nt - 1) / 2):1])
    wav2 <- wave_bases(wname = wname, k = k, scales = scales, wparam = wparam)
    coefs <- mvfft(t(sweep(wav2, MARGIN = 2, f, `*`)), inverse = TRUE) / length(f)
  }

  return(coefs)
}

#' @title Extracts the center of a vector
#'
#' @description
#' This function is an internal function which extracts from a vector \code{x},
#' the center of the vector of length \code{n}. It emulates the Matlab(R) function \code{wkeep}.
#' This function is used by the cwt_wst function when the HAAR wavelet is selected.
#'


core <- function(x, n) {
  if (n > length(x)) {
    stop("Error. The length of the vector x should be greater than n.")
  }
  if (n == length(x)) {
    extracted <- x
  } else {
    if (length(x) %% 2 == 0){ # Length of x is even
      center1 <- length(x) / 2
      center2 <- length(x) / 2 + 1
      if (n %% 2 == 1) { # n is odd
        n1 <- center1 - floor(n / 2)
        n2 <- center2 + floor(n / 2) - 1
      } else { # n is even
        n1 <- center1 - (n / 2 - 1)
        n2 <- center2 + (n / 2 - 1)
      }
    } else { # Length of x is odd
      center <- floor(length(x) / 2) + 1
      if (n %% 2 == 1) { # n is odd
        n1 <- center - floor(n / 2)
        n2 <- center + floor(n / 2)
      } else { # n is even
        n1 <- center - n / 2
        n2 <- center + (n / 2 - 1)
      }
    }
    extracted <- x[n1:n2]
  }
  return(extracted)
}

