#' @title Windowed scalograms of a signal
#'
#' @description This function computes the normalized windowed scalograms of a signal for
#' the scales given. It is computed using time windows with radius \code{windowrad}
#' centered at a vector of central times with increment of time \code{delta_t}. It is
#' important to note that the notion of scalogram here is analogous to the spectrum of the
#' Fourier transform. It gives the contribution of each scale to the total energy of the
#' signal. For each scale \eqn{s} and central time \eqn{tc}, it is defined as the square
#' root of the integral of the squared modulus of the wavelet transform w.r.t the time
#' variable \eqn{t}, i.e.
#'
#'   \deqn{WS_{windowrad}(tc,s):=
#'   (\int_{tc-windowrad}^{tc+windowrad}|Wf(t,s)|^2 dt)^{1/2}.}{WS_{windowrad}(tc,s):=
#'   (integral_{tc-windowrad}^{tc+windowrad}|Wf(t,s)|^2 dt)^{1/2}.}
#'
#' "Normalized" means that the windowed scalograms are divided by the square root of the
#' length of the respective time windows in order to be comparable between them.
#'
#'
#' @usage windowed_scalogram(signal,
#'                           dt = 1,
#'                           scales = NULL,
#'                           powerscales = TRUE,
#'                           windowrad = NULL,
#'                           delta_t = NULL,
#'                           wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
#'                           wparam = NULL,
#'                           waverad = NULL,
#'                           border_effects = c("BE", "INNER", "PER", "SYM"),
#'                           energy_density = TRUE,
#'                           makefigure = TRUE,
#'                           time_values = NULL,
#'                           figureperiod = TRUE)
#'
#' @param signal A vector containing the signal whose windowed scalogram is wanted.
#' @param dt Numeric. The time step of the signal.
#' @param scales A vector containing the wavelet scales at wich
#' the windowed scalogram is computed. This can be either a vector with all the scales, or
#' (if \code{powerscales} is TRUE) following Torrence and Compo 1998, a vector of three
#' elements with the minimum scale, the maximum scale and the number of suboctaves per
#' octave. If NULL, they are automatically computed.
#' @param powerscales Logical. If TRUE (default), construct power 2 scales from
#' \code{scales}. If \code{scales} is NULL, they are automatically computed.
#' @param windowrad Integer. Time radius for the windows, measured in \code{dt}. By
#' default, it is set to \eqn{celing(length(signal1) / 20)}.
#' @param delta_t Integer. Increment of time for the construction of windows central
#' times, measured in \code{dt}. By default, it is set to
#' \eqn{celing(length(signal) / 256)}.
#' @param wname A string, equal to "MORLET", "DOG", "PAUL", "HAAR" or "HAAR2". The
#' difference between "HAAR" and "HAAR2" is that "HAAR2" is more accurate but slower.
#' @param wparam The corresponding nondimensional parameter for the wavelet function
#' (Morlet, DoG or Paul).
#' @param waverad Numeric. The radius of the wavelet used in the computations for the cone
#' of influence. If it is not specified, it is computed by \code{wavelet_radius} for DoG
#' and Paul wavelets. For Haar and Morlet it is assumed to be 1 and 3 respectively.
#' @param border_effects String, equal to "BE", "INNER", "PER" or "SYM", which indicates
#' how to manage the border effects which arise usually when a convolution is performed on
#' finite-lenght signals.
#' \itemize{
#' \item "BE": With border effects, padding time series with zeroes.
#' \item "INNER": Normalized inner scalogram with security margin adapted for each
#'     different scale. Although there are no border effects, it is shown as a regular COI
#'     the zone in which the length of \eqn{J(s)} (see Benítez et al. 2010) is smaller and
#'     it has to be normalized.
#' \item "PER": With border effects, using boundary wavelets (periodization of the
#' original time series).
#' \item "SYM": With border effects, using a symmetric catenation of the original time
#' series.
#' }
#' @param energy_density Logical. If TRUE (default), divide the scalograms by the square
#' root of the scales for convert them into energy density.
#' @param makefigure Logical. If TRUE (default), a figure with the scalograms is plotted.
#' @param time_values A numerical vector of length \code{length(signal)} containing custom
#' time values for the figure. If NULL (default), it will be computed starting at 0.
#' @param figureperiod Logical. If TRUE (default), periods are used in the figure instead
#' of scales.
#'
#' @return A list with the following fields:
#' \itemize{
#' \item \code{wsc}: A matrix of size \code{length(tcentral)} x \code{length(scales)}
#' containing the values of the windowed scalograms at each scale and at each time window.
#' \item \code{tcentral}: The vector of central times at which the windows are centered.
#' \item \code{scales}: The vector of the scales.
#' \item \code{windowrad}: Radius for the time windows, measured in \code{dt}.
#' \item \code{fourier_factor}: A factor for converting scales to periods.
#' \item \code{coi_maxscale}: A vector of length \code{length(tcentral)} containing the
#' values of the maximum scale from which there are border effects for the respective
#' central time.
#' }
#'
#' @examples
#' dt <- 0.1
#' time <- seq(0, 50, dt)
#' signal <- c(sin(pi * time), sin(pi * time / 2))
#' wscalog <- windowed_scalogram(signal = signal, dt = dt)
#'
#'
#' @section References:
#'
#' C. Torrence, G. P. Compo. A practical guide to wavelet analysis. B. Am. Meteorol. Soc.
#' 79 (1998), 61–78.
#'
#' V. J. Bolós, R. Benítez, R. Ferrer, R. Jammazi. The windowed scalogram difference: a
#' novel wavelet tool for comparing time series. Appl. Math. Comput., 312 (2017), 49-65.
#'
#' R. Benítez, V. J. Bolós, M. E. Ramírez. A wavelet-based tool for studying
#' non-periodicity. Comput. Math. Appl. 60 (2010), no. 3, 634-641.
#'
#' @export
#'


windowed_scalogram <-
  function(signal,
           dt = 1,
           scales = NULL,
           powerscales = TRUE,
           windowrad = NULL,
           delta_t = NULL,
           wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
           wparam = NULL,
           waverad = NULL,
           border_effects = c("BE", "INNER", "PER", "SYM"),
           energy_density = TRUE,
           makefigure = TRUE,
           time_values = NULL,
           figureperiod = TRUE) {

  #  require(zoo)
  #  require(Matrix)

    wname <- toupper(wname)
    wname <- match.arg(wname)

    if (is.null(waverad)) {
      if (wname == "MORLET") {
        waverad <- 3
      } else if ((wname == "HAAR") || (wname == "HAAR2")) {
        waverad <- 0.5
      } else {
        waverad <- wavelet_radius(wname = wname, wparam = wparam)
        waverad <- waverad$left
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

    if (is.null(delta_t)) {
      delta_t <- ceiling(nt / 256)
    }

    if (is.null(windowrad)) {
      windowrad <- ceiling(nt / 20)
    } else {
      windowrad <- min(windowrad, floor((nt - 1) / 2))
    }

    if (is.null(scales)) {
      scmax <- floor((nt - 2 * windowrad) / (2 * waverad))
      if (powerscales) {
        scales <- pow2scales(c(2, scmax, ceiling(256 / (log2(scmax) - 1))))
      } else {
        scales <- seq(2, scmax, by = scmax / 256)
      }
      scalesdt <- scales * dt
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
      scalesdt <- scales
      scales <- scales / dt
    }

    ns <- length(scales)

    cwt <- cwt_wst(signal = signal,
                   dt = dt,
                   scales = scalesdt,
                   powerscales = FALSE,
                   wname = wname,
                   wparam = wparam,
                   waverad = waverad,
                   border_effects = border_effects_cwt,
                   makefigure = FALSE)

    coefs <- cwt$coefs

    if (border_effects == "INNER") {

      wrs <- ceiling(waverad * scales)
      tcentral_ini <- max(1 + windowrad, 1 + wrs[1] - windowrad)
      tcentral_end <- min(nt - windowrad, nt - wrs[1] + windowrad)
      if (tcentral_ini > tcentral_end) {
        stop("We need a larger signal")
      }
      tcentral <- seq(from = tcentral_ini, to = tcentral_end, by = delta_t)
      ntcentral <- length(tcentral)
      wsc <- matrix(NA, nrow = ntcentral, ncol = ns)

      abscoefs2 <- abs(coefs) ^ 2

      # Regular version
      for (i in 1:ntcentral) {
        for (j in 1:ns) {
          t_ini <- max(tcentral[i] - windowrad, 1 + wrs[j])
          t_end <- min(tcentral[i] + windowrad, nt - wrs[j])
          if (t_ini <= t_end) {
            wsc[i, j] <- sqrt(abs(sum(abscoefs2[t_ini:t_end, j]))) # abs: sometimes wsc is negative due to numerical errors
            wsc[i, j] <- wsc[i, j] / sqrt(t_end - t_ini + 1) # Normalization
          }
        }
      }

      wsc <- as.matrix(wsc)

    } else {

      tcentral_ini <- 1 + windowrad
      tcentral_end <- nt - windowrad
      tcentral <- seq(from = tcentral_ini, to = tcentral_end, by = delta_t)
      ntcentral <- length(tcentral)
      wsc <- matrix(0, nrow = ntcentral, ncol = ns)

      abscoefs2 <- abs(coefs) ^ 2

      if (delta_t < windowrad) { # Fast version
        for (j in 1:ns) {
          wsc[1, j] <- sum(abscoefs2[1:(1 + 2 * windowrad), j])
          for (i in 2:ntcentral) {
            wsc[i, j] <- wsc[i-1, j] - sum(abscoefs2[(tcentral[i] - windowrad - delta_t):(tcentral[i] - windowrad - 1), j]) + sum(abscoefs2[(tcentral[i] + windowrad - delta_t + 1):(tcentral[i] + windowrad), j])
          }
        }
      } else { # Regular version
        for (i in 1:ntcentral) {
          for (j in 1:ns) {
            wsc[i, j] <- sum(abscoefs2[(tcentral[i] - windowrad):(tcentral[i] + windowrad), j])
          }
        }
      }
      wsc <- as.matrix(sqrt(abs(wsc))) # abs: sometimes wsc is negative due to numerical errors
      wsc <- wsc / sqrt(2 * windowrad + 1) # Normalization

    }

    # COI
    coi_maxscale <- numeric(ntcentral)
    for (i in 1:ntcentral) {
      coi_maxscale[i] <- dt * min(tcentral[i] - windowrad - 1, nt - tcentral[i] - windowrad) / waverad
    }

    tcentraldt <- tcentral * dt

    # Energy density
    if (energy_density) {
      wsc <- t(t(wsc) / sqrt(scalesdt))
    }

    fourier_factor <- cwt$fourier_factor

    # Make figure
    if (makefigure) {

      if (figureperiod) {
        Y <- fourier_factor * scalesdt
        coi <- fourier_factor * coi_maxscale
        Yname <- "Period"
      } else {
        Y <- scalesdt
        coi <- coi_maxscale
        Yname <- "Scale"
      }

      if (is.null(time_values)) {
        X <- tcentraldt
      } else {
        if (length(time_values) != nt) {
          warning("Invalid length of time_values vector. Changing to default.")
          X <- tcentraldt
        } else {
          X <- time_values[tcentral]
        }
      }

      wavPlot(
        Z = wsc,
        X = X,
        Y = Y,
        Ylog = powerscales,
        coi = coi,
        Xname = "Time",
        Yname = Yname,
        Zname = "Windowed Scalogram"
      )

    }

    return(list(
      wsc = wsc,
      tcentral = tcentraldt,
      scales = scalesdt,
      windowrad = windowrad,
      fourier_factor = fourier_factor,
      coi_maxscale = as.numeric(coi_maxscale)
    ))
  }
