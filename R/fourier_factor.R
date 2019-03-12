#' @title Fourier factor of a wavelet
#'
#' @description This function computes the Fourier factor of a wavelet, according to
#' Torrence and Compo (1998).
#'
#' @usage fourier_factor(wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
#'                       wparam = NULL)
#'
#' @param wname A string, equal to "MORLET", "DOG", "PAUL", "HAAR" or "HAAR2" that
#' determines the wavelet function.
#' @param wparam The corresponding nondimensional parameter for the wavelet function
#' (Morlet, DoG or Paul).
#'
#' @return The numeric value of the Fourier factor.
#'
#' @examples
#' ff <- fourier_factor(wname = "DOG", wparam = 6)
#'
#' @section References:
#'
#' C. Torrence, G. P. Compo. A practical guide to wavelet analysis. B. Am. Meteorol. Soc.
#' 79 (1998), 61–78.
#'
#' @export
#'

fourier_factor <-
  function(wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
           wparam = NULL) {

    wname <- toupper(wname)
    wname <- match.arg(wname)

    if (wname == "MORLET") {
      if (is.null(wparam)) {
        wparam <- 6
      }
      ff <- (4 * pi) / (wparam + sqrt(2 + wparam ^ 2))
    } else if (wname == "PAUL") {
      if (is.null(wparam)) {
        wparam <- 4
      }
      ff <- 4 * pi / (2 * wparam + 1)
    } else if (wname == "DOG") {
      if (is.null(wparam)) {
        wparam <- 2
      }
      ff <- 2 * pi * sqrt(2 / (2 * wparam + 1))
    } else {  # HAAR
      ff <- 1
    }

    return(ff)

}
