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
#'                           scales = NULL,
#'                           powerscales = TRUE,
#'                           windowrad = NULL,
#'                           delta_t = NULL,
#'                           wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
#'                           wparam = NULL,
#'                           border_effects = c("BE", "INNER", "PER", "SYM"),
#'                           energy_density = TRUE,
#'                           makefigure = FALSE)
#'
#' @param signal A vector containing the signal whose wavelet transform is wanted. The
#' unit of time is taken as the time difference between two consecutive data.
#' @param scales A vector containing the wavelet scales measured in units of time. This
#' can be either a vector with all the scales, or (if \code{powerscales = TRUE}) following
#' Torrence and Compo 1998, a vector of three elements with the minimum scale, the maximum
#' scale and the number of suboctaves per octave.
#' @param powerscales Logical. Construct power 2 scales.
#' @param windowrad Numeric. Radius for the time windows.
#' @param delta_t Numeric. Increment of time for the construction of windows central
#' times.
#' @param wname A string, equal to "MORLET", "DOG", "PAUL", "HAAR" or "HAAR2". The
#' difference between "HAAR" and "HAAR2" is that "HAAR2" is more accurate but slower.
#' @param wparam Numeric. Parameters of the corresponding wavelet.
#' @param border_effects String, equal to "BE", "INNER", "PER" or "SYM",
#' which indicates how to manage the border effects which arise usually when a convolution
#' is performed on finite-lenght signals.
#' \itemize{
#' \item "BE": With border effects, padding time series with zeroes.
#' \item "INNER": Normalized inner scalogram with security margin adapted for each
#'     different scale. Altough there are no border effects, it is shown as a regular COI
#'     the zone in which the length of \eqn{J(s)} (see Benítez et al. 2010) is smaller and
#'     it has to be normalized.
#' \item "PER": With border effects, using boundary wavelets (periodization of the
#' original time series).
#' \item "SYM": With border effects, using a symmetric catenation of the original time
#' series.
#' }
#' @param energy_density Logical. Divide the scalograms by sqrt(scales) for convert them
#' into energy density (defaults to TRUE).
#' @param makefigure Logical. Plots a figure with the scalogram.
#'
#'
#' @return A list with the following fields:
#'
#' \code{wsc}: A matrix of size \code{length(tcentral)}x\code{length(scales)} containing
#' the values of the windowed scalograms at each scale and at each time window.
#'
#' \code{tcentral}: The vector of central times.
#'
#' \code{scales}: The vector of the scales.
#'
#' \code{coi_maxscale}: A vector of length \code{length(tcentral)} containing the values
#' of the maximum scale from which there are border effects for the respective windowed
#' scalogram.
#'
#' @examples
#' time <- 1:500
#' signal <- c(sin(pi * time / 8), sin(pi * time / 16))
#' wscalog <- windowed_scalogram(signal = signal, makefigure = TRUE)
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
           scales = NULL,
           powerscales = TRUE,
           windowrad = NULL,
           delta_t = NULL,
           wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
           wparam = NULL,
           border_effects = c("BE", "INNER", "PER", "SYM"),
           energy_density = TRUE,
           makefigure = FALSE) {

  #  require(zoo)
  #  require(Matrix)

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
    } else {
      if (powerscales) {
        scales <- pow2scales(scales)
      } else {
        if (is.unsorted(scales))
          scales <- sort(scales)
      }
    }

    ns <- length(scales)

    if (border_effects == "BE") {

      tcentral_ini <- 1 + windowrad
      tcentral_end <- nt - windowrad
      tcentral <- seq(from = tcentral_ini, to = tcentral_end, by = delta_t)
      ntcentral <- length(tcentral)
      wsc <- matrix(0, nrow = ntcentral, ncol = ns)

      coefs <- cwt_wst(scales = scales, signal = signal , wname = wname, wparam = wparam)
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

    } else if (border_effects == "INNER") {

      wrs <- ceiling(waverad * scales)
      tcentral_ini <- max(1 + windowrad, 1 + wrs[1] - windowrad)
      tcentral_end <- min(nt - windowrad, nt - wrs[1] + windowrad)
      if (tcentral_ini > tcentral_end) {
        stop("We need a larger signal")
      }
      tcentral <- seq(from = tcentral_ini, to = tcentral_end, by = delta_t)
      ntcentral <- length(tcentral)
      wsc <- matrix(NA, nrow = ntcentral, ncol = ns)

      coefs <- cwt_wst(scales = scales, signal = signal , wname = wname, wparam = wparam)
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

    } else if (border_effects == "PER") {

      tcentral_ini <- 1 + windowrad
      tcentral_end <- nt - windowrad
      tcentral <- seq(from = tcentral_ini, to = tcentral_end, by = delta_t)
      ntcentral <- length(tcentral)
      wsc <- matrix(0, nrow = ntcentral, ncol = ns)

      ndatoscat <- ceiling(ceiling(waverad * scales[ns]) / nt) # Number of catenations = ndatoscat*2+1
      ntcat <- (ndatoscat * 2 + 1) * nt
      tcatcentral_ini <- ndatoscat * nt + tcentral_ini
      tcatcentral_end <- ndatoscat * nt + tcentral_end
      tcatcentral <- seq(from = tcatcentral_ini, to = tcatcentral_end, by = delta_t)

      datoscat <- matrix(0, nrow = 1, ncol = ntcat)
      for (itcat in 1:ntcat) {
        datoscat[itcat] <- signal[itcat - floor((itcat - 1) / nt) * nt]
      }

      coefs <- cwt_wst(scales = scales, signal = datoscat, wname = wname, wparam = wparam)
      abscoefs2 <- abs(coefs) ^ 2

      if (delta_t < windowrad) { # Fast version
        for (j in 1:ns) {
          wsc[1, j] <- sum(abscoefs2[(tcatcentral_ini - windowrad):(tcatcentral_ini + windowrad), j])
          for (i in 2:ntcentral) {
            wsc[i, j] <- wsc[i - 1, j] - sum(abscoefs2[(tcatcentral[i] - windowrad - delta_t):(tcatcentral[i] - windowrad - 1), j]) + sum(abscoefs2[(tcatcentral[i] + windowrad - delta_t + 1):(tcatcentral[i] + windowrad), j])
          }
        }
      } else { # Regular version
        for (i in 1:ntcentral) {
          for (j in 1:ns) {
            wsc[i, j] <- sum(abscoefs2[(tcatcentral[i] - windowrad):(tcatcentral[i] + windowrad), j])
          }
        }
      }

      wsc <- as.matrix(sqrt(abs(wsc))) # abs: sometimes wsc is negative due to numerical errors
      wsc <- wsc / sqrt(2 * windowrad + 1) # Normalization

    } else {

      tcentral_ini <- 1 + windowrad
      tcentral_end <- nt - windowrad
      tcentral <- seq(from = tcentral_ini, to = tcentral_end, by = delta_t)
      ntcentral <- length(tcentral);
      wsc <- matrix(0, nrow = ntcentral, ncol = ns)

      ndatoscat <- ceiling(ceiling(waverad * scales[ns]) / nt) # Number of catenations = ndatoscat*2+1
      tcatcentral_ini <- ndatoscat * nt + tcentral_ini
      tcatcentral_end <- ndatoscat * nt + tcentral_end
      tcatcentral <- seq(from = tcatcentral_ini, to = tcatcentral_end, by = delta_t)

      datos_sim <- signal
      aux <- rev(signal)
      for (icat in 1:ndatoscat) {
        datos_sim <- c(aux, datos_sim, aux)
        aux <- rev(aux)
      }

      coefs <- cwt_wst(scales = scales, signal = datos_sim, wname = wname, wparam = wparam)
      abscoefs2 <- abs(coefs) ^ 2

      if (delta_t < windowrad) { # Fast version
        for (j in 1:ns) {
          wsc[1, j] <- sum(abscoefs2[(tcatcentral_ini - windowrad):(tcatcentral_ini + windowrad), j])
          for (i in 2:ntcentral) {
            wsc[i, j] <- wsc[i - 1, j] - sum(abscoefs2[(tcatcentral[i] - windowrad - delta_t):(tcatcentral[i] - windowrad - 1), j])  + sum(abscoefs2[(tcatcentral[i] + windowrad - delta_t + 1):(tcatcentral[i] + windowrad), j])
          }
        }
      } else { # Regular version
        for (i in 1:ntcentral) {
          for (j in 1:ns) {
            wsc[i, j] <- sum(abscoefs2[(tcatcentral[i] - windowrad):(tcatcentral[i] + windowrad), j])
          }
        }
      }

      wsc <- as.matrix(sqrt(abs(wsc))) # abs: sometimes wsc is negative due to numerical errors
      wsc <- wsc/sqrt(2*windowrad+1) # Normalization
    }

    # COI
    coi_maxscale <- numeric(ntcentral)
    for (i in 1:ntcentral) {
      coi_maxscale[i] <- min(tcentral[i] - windowrad - 1, nt - tcentral[i] - windowrad) / waverad
    }

    # Energy density
    if (energy_density) {
      wsc <- t(t(wsc) / sqrt(scales))
    }

    # Make figure
    if (makefigure) {
      wavPlot(
        Z = wsc,
        X = tcentral,
        Y = scales,
        Ylog = powerscales,
        coi = coi_maxscale,
        Xname = "Time",
        Yname = "Scale",
        Zname = "WSC"
      )
    }

    return(list(
      wsc = wsc,
      tcentral = tcentral,
      scales = scales,
      coi_maxscale = as.numeric(coi_maxscale)
    ))
  }
