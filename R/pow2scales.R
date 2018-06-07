#' @title Power 2 scales
#'
#' @description This function constructs power 2 scales from a vector of three elements
#' with the minimum scale, the maximum scale and the number of suboctaves per octave
#' (following Torrence and Compo 1998).
#'
#' @usage pow2scales(scales)
#'
#' @param scales A vector of three elements with the minimum scale, the maximum scale and
#' the number of suboctaves per octave.
#'
#' @return A vector with all the scales.
#'
#' @examples
#' scales <- pow2scales(c(2,128,8))
#'
#' @section References:
#'
#' C. Torrence, G. P. Compo. A practical guide to wavelet analysis. B. Am. Meteorol. Soc.
#' 79 (1998), 61–78.
#'
#' @export
#'

pow2scales <- function(scales) {

  if(length(scales) == 3) {
    scmin <- scales[1]
    scmax <- scales[2]
    J1 <- log2(scmax / scmin)
    Dj <- scales[3]
    scales <- scmin * 2 ^ ((0:floor(J1 * Dj)) / Dj)
  } else {
    stop("For power scales, a vector of length three should be provided.")
  }

 return(scales)

}
