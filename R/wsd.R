#' @title Windowed Scalogram Difference
#'
#' @description This function computes the Windowed Scalogram Difference of two signals.
#' The definition and details can be found in (Bolós et al. 2017).
#'
#' @usage wsd(signal1,
#'            signal2,
#'            dt = 1,
#'            scaleparam = NULL,
#'            windowrad = NULL,
#'            rdist = NULL,
#'            delta_t = NULL,
#'            normalize = c("NO", "ENERGY", "MAX", "SCALE"),
#'            refscale = NULL,
#'            wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
#'            wparam = NULL,
#'            waverad = NULL,
#'            border_effects = c("BE", "INNER", "PER", "SYM"),
#'            mc_nrand = 0,
#'            commutative = TRUE,
#'            wscnoise = 0.02,
#'            compensation = 0,
#'            energy_density = TRUE,
#'            parallel = FALSE,
#'            makefigure = TRUE,
#'            time_values = NULL,
#'            figureperiod = TRUE,
#'            xlab = "Time",
#'            ylab = NULL,
#'            main = "-log2(WSD)")
#'
#' @param signal1 A vector containing the first signal.
#' @param signal2 A vector containing the second signal (its length should be equal to
#' that of \code{signal1}).
#' @param dt Numeric. The time step of the signals.
#' @param scaleparam A vector of three elements with the minimum scale, the maximum scale
#' and the number of suboctaves per octave for constructing power 2 scales (following
#' Torrence and Compo 1998). If NULL, they are automatically computed.
#' @param windowrad Integer. Time radius for the windows, measured in \code{dt}. By
#' default, it is set to \eqn{ceiling(length(signal1) / 20)}.
#' @param rdist Integer. Log-scale radius for the windows measured in suboctaves. By
#' default, it is set to \eqn{ceiling(length(scales) / 20)}.
#' @param delta_t Integer. Increment of time for the construction of windows central
#' times, measured in \code{dt}. By default, it is set to
#' \eqn{ceiling(length(signal1) / 256)}.
#' @param normalize String, equal to "NO", "ENERGY", "MAX" or "SCALE". If "ENERGY", signals are
#' divided by their respective energies. If "MAX", each signal is divided by the maximum
#' value attained by its scalogram. In these two cases, \code{energy_density} must be TRUE.
#' Finally, if "SCALE", each signal is divided by their scalogram value at scale
#' \code{refscale}.
#' @param refscale Numeric. The reference scale for \code{normalize}.
#' @param wname A string, equal to "MORLET", "DOG", "PAUL", "HAAR" or "HAAR2". The
#' difference between "HAAR" and "HAAR2" is that "HAAR2" is more accurate but slower.
#' @param wparam The corresponding nondimensional parameter for the wavelet function
#' (Morlet, DoG or Paul).
#' @param waverad Numeric. The radius of the wavelet used in the computations for the cone
#' of influence. If it is not specified, it is asumed to be \eqn{\sqrt{2}} for Morlet and DoG,
#' \eqn{1/\sqrt{2}} for Paul and 0.5 for Haar.
#' @param border_effects String, equal to "BE", "INNER", "PER" or "SYM",
#' which indicates how to manage the border effects which arise usually when a convolution
#' is performed on finite-lenght signals.
#' \itemize{
#' \item "BE": With border effects, padding time series with zeroes.
#' \item "INNER": Normalized inner scalogram with security margin adapted for each
#' different scale.
#' \item "PER": With border effects, using boundary wavelets (periodization of the
#' original time series).
#' \item "SYM": With border effects, using a symmetric catenation of the original time
#' series.
#' }
#' @param mc_nrand Integer. Number of Montecarlo simulations to be performed in order to
#' determine the 95\% and 5\% significance contours.
#' @param commutative Logical. If TRUE (default) the commutative windowed scalogram
#' difference. Otherwise a non-commutative (but simpler) version is computed (see Bolós et
#' al. 2017).
#' @param wscnoise Numeric in \eqn{[0,1]}. If a (windowed) scalogram takes values close to
#' zero, some problems may appear because we are considering relative differences.
#' Specifically, we can get high relative differences that in fact are not relevant, or
#' even divisions by zero.
#'
#'   If we consider absolute differences this would not happen but, on the other hand,
#'   using absolute differences is not appropriate for scalogram values not close to zero.
#'
#'   So, the parameter \code{wscnoise} stablishes a threshold for the scalogram values
#'   above which a relative difference is computed, and below which a difference
#'   proportional to the absolute difference is computed (the proportionality factor is
#'   determined by requiring continuity).
#'
#'   Finally, \code{wscnoise} can be interpreted as the relative amplitude of the noise in
#'   the scalograms and is chosen in order to make a relative (\eqn{= 0}), absolute
#'   (\eqn{= 1}) or mix (in \eqn{(0,1)}) difference between scalograms. Default value is
#'   set to \eqn{0.02}.
#' @param compensation Numeric. It is an alternative to \code{wscnoise} for
#' preventing numerical errors or non-relevant high relative differences when scalogram
#' values are close to zero (see Bolós et al. 2017). It should be a non-negative
#' relatively small value.
#' @param energy_density Logical. If TRUE (default), divide the scalograms by the square
#' root of the scales for convert them into energy density. Note that it does not affect
#' the results if \code{wscnoise} \eqn{= 0}.
#' @param parallel Logical. If TRUE, it uses function \code{parApply} from package
#' \code{parallel} for the Montecarlo simulations. When FALSE (default) it uses the normal
#' \code{apply} function.
#' @param makefigure Logical. If TRUE (default), a figure with the WSD is plotted.
#' @param time_values A numerical vector of length \code{length(signal)} containing custom
#' time values for the figure. If NULL (default), it will be computed starting at 0.
#' @param figureperiod Logical. If TRUE (default), periods are used in the figure instead
#' of scales.
#' @param xlab A string giving a custom X axis label.
#' @param ylab A string giving a custom Y axis label. If NULL (default) the Y label is
#' either "Scale" or "Period" depending on the value of \code{figureperiod}.
#' @param main A string giving a custom main title for the figure.
#'
#' @importFrom parallel parApply detectCores makeCluster stopCluster
#' @importFrom abind abind
#'
#' @return A list with the following fields:
#' \itemize{
#' \item \code{wsd}: A matrix of size \code{length(tcentral)} x \code{length(scales)}
#' containing the values of the windowed scalogram differences at each scale and at each
#' time window.
#' \item \code{tcentral}: The vector of central times used in the computations of the
#' windowed scalograms.
#' \item \code{scales}: The vector of scales.
#' \item \code{windowrad}: Radius for the time windows of the windowed scalograms,
#' measured in \code{dt}.
#' \item \code{rdist}: The log-scale radius for the windows measured in suboctaves.
#' \item \code{signif95}: A logical matrix of size \code{length(tcentral)} x
#' \code{length(scales)}. If TRUE, the corresponding point of the \code{wsd} matrix is in
#' the 95\% significance.
#' \item \code{signif05}: A logical matrix of size \code{length(tcentral)} x
#' \code{length(scales)}. If TRUE, the corresponding point of the \code{wsd} matrix is in
#' the 5\% significance.
#' \item \code{fourierfactor}: A factor for converting scales into periods.
#' \item \code{coi_maxscale}: A vector of length \code{length(tcentral)} containing the
#' values of the maximum scale from which there are border effects for the respective
#' central time.
#' }
#'
#' @examples
#'
#' nt <- 1500
#' time <- 1:nt
#' sd_noise <-  0.2 #% In Bolós et al. 2017 Figure 1, sd_noise = 1.
#' signal1 <- rnorm(n = nt, mean = 0, sd = sd_noise) + sin(time / 10)
#' signal2 <- rnorm(n = nt, mean = 0, sd = sd_noise) + sin(time / 10)
#' signal2[500:1000] = signal2[500:1000] + sin((500:1000) / 2)
#' \dontrun{
#' wsd <- wsd(signal1 = signal1, signal2 = signal2)
#' }
#'
#' @section References:
#' C. Torrence, G. P. Compo. A practical guide to wavelet analysis. B. Am. Meteorol. Soc.
#' 79 (1998), 61–78.
#'
#' V. J. Bolós, R. Benítez, R. Ferrer, R. Jammazi. The windowed scalogram difference: a
#' novel wavelet tool for comparing time series. Appl. Math. Comput., 312 (2017), 49-65.
#'
#' @export
#'

wsd <-
  function(signal1,
           signal2,
           dt = 1,
           scaleparam = NULL,
           windowrad = NULL,
           rdist = NULL,
           delta_t = NULL,
           normalize = c("NO", "ENERGY", "MAX", "SCALE"),
           refscale = NULL,
           wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
           wparam = NULL,
           waverad = NULL,
           border_effects = c("BE", "INNER", "PER", "SYM"),
           mc_nrand = 0,
           commutative = TRUE,
           wscnoise = 0.02,
           compensation = 0,
           energy_density = TRUE,
           parallel = FALSE,
           makefigure = TRUE,
           time_values = NULL,
           figureperiod = TRUE,
           xlab = "Time",
           ylab = NULL,
           main = "-log2(WSD)") {

  #  require(abind)

    wname <- toupper(wname)
    wname <- match.arg(wname)

    border_effects <- toupper(border_effects)
    border_effects <- match.arg(border_effects)

    normalize <- toupper(normalize)
    normalize <- match.arg(normalize)

    if (is.null(waverad)) {
      if ((wname == "MORLET") || (wname == "DOG")) {
        waverad <- sqrt(2)
      } else if (wname == "PAUL") {
        waverad <- 1 / sqrt(2)
      } else { # HAAR
        waverad <- 0.5
      }
    }

    nt <- length(signal1)
    if (nt != length(signal2)) {
      stop("Signals must have the same length.")
    }

    if (is.null(delta_t)) {
      delta_t <- ceiling(nt / 256)
    }

    if (is.null(windowrad)) {
      windowrad <- ceiling(nt / 20)
    } else if (windowrad > floor((nt - 1) / 2)) {
      windowrad <- floor((nt - 1) / 2)
      print("Windowrad re-adjusted.")
    }

    fourierfactor <- fourier_factor(wname = wname, wparam = wparam)

    if (is.null(scaleparam)) {
      scmin <- 2 / fourierfactor
      scmax <- floor((nt - 2 * windowrad) / (2 * waverad))
      scaleparam <- c(scmin * dt, scmax * dt, ceiling(256 / log2(scmax / scmin)))
    }

    scalesdt <- pow2scales(scaleparam)
    scales <- scalesdt / dt
    ns <- length(scales)

    if (is.null(rdist)) {
      rdist <- ceiling(ns / 20)
    }

    ### Normalization & Scalograms

    if (normalize != "NO") {

      if ((!energy_density) && (normalize != "SCALE")) {
        warning("For energy or max normalization, energy_density is set to TRUE.")
        energy_density <- TRUE
      }

      if (normalize == "SCALE") {
        if (is.null(refscale)) {
          stop("No refscale parameter specified for normalization.")
        } else {
          scales_aux <- refscale
        }
      } else {
        scales_aux <- scalesdt
      }

      scalog1 <-
        scalogram(signal = signal1, dt = dt,
                  scales = scales_aux, powerscales = FALSE,
                  border_effects = border_effects,
                  energy_density = energy_density,
                  wname = wname, wparam = wparam, waverad = waverad,
                  makefigure = FALSE)
      scalog2 <-
        scalogram(signal = signal2, dt = dt,
                  scales = scales_aux, powerscales = FALSE,
                  border_effects = border_effects,
                  energy_density = energy_density,
                  wname = wname, wparam = wparam, waverad = waverad,
                  makefigure = FALSE)

      # Normalization

      if (normalize == "ENERGY") {
        signal1 <- signal1 / scalog1$energy
        signal2 <- signal2 / scalog2$energy
      } else { # MAX or SCALE
        signal1 <- signal1 / max(scalog1$scalog)
        signal2 <- signal2 / max(scalog2$scalog)
      }

    }

    ### Windowed scalograms

    wsc1 <-
      windowed_scalogram(
        signal = signal1, dt = dt,
        scales = scalesdt, powerscales = FALSE,
        border_effects = border_effects,
        energy_density = energy_density,
        windowrad = windowrad,
        delta_t = delta_t,
        wname = wname, wparam = wparam, waverad = waverad,
        makefigure = FALSE)
    wsc2 <-
      windowed_scalogram(
        signal = signal2, dt = dt,
        scales = scalesdt, powerscales = FALSE,
        border_effects = border_effects,
        energy_density = energy_density,
        windowrad = windowrad,
        delta_t = delta_t,
        wname = wname, wparam = wparam, waverad = waverad,
        makefigure = FALSE)

    nwsc <- length(wsc1$tcentral)

    # Compensation

    fc <- 1 - compensation / max(wsc1$wsc, wsc2$wsc, na.rm = TRUE)
    wsc1$wsc <- compensation + fc * wsc1$wsc
    wsc2$wsc <- compensation + fc * wsc2$wsc

    ### Relative distances between scalograms

    wsd <- matrix(NA, nrow = nwsc, ncol = ns)

    noise1 <- wscnoise * max(wsc1$wsc, na.rm = TRUE)

    if (commutative) {
      noise2 <- wscnoise * max(wsc2$wsc, na.rm = TRUE)
      A <- 0.25 * ((wsc1$wsc - wsc2$wsc) / pmax(wsc1$wsc, noise1) + (wsc1$wsc - wsc2$wsc) / pmax(wsc2$wsc, noise2)) ^ 2
    } else {
      A <- ((wsc1$wsc - wsc2$wsc) / pmax(wsc1$wsc, noise1)) ^ 2
    }

    for (i in 1:nwsc) {
      maxsct <- max(which(!is.na(A[i, ])))
      for (j in 1:maxsct) {
        kmin <- max(1, j - rdist)
        kmax <- min(maxsct, j + rdist)
        aux <- (2 * rdist + 1)/(kmax - kmin + 1)
        wsd[i, j] <- sum(A[i, kmin:kmax]) * aux
      }
    }

    wsd <- sqrt(wsd)

    ### Montecarlo

    if (mc_nrand > 0) {
      mean1 <- mean(signal1)
      mean2 <- mean(signal2)
      sd1 <- sd(signal1)
      sd2 <- sd(signal2)
      distrnd <- array(0, dim = c(nwsc, ns, mc_nrand))

      print("Montecarlo...")
      rndSignals1 <- matrix(rnorm(mc_nrand * nt, mean = mean1, sd = sd1), nrow = mc_nrand, ncol = nt)
      rndSignals2 <- matrix(rnorm(mc_nrand * nt, mean = mean2, sd = sd2), nrow = mc_nrand, ncol = nt)
      rndSignals <- array(c(rndSignals1, rndSignals2), dim = c(mc_nrand, nt, 2))

      if(parallel) {
      #  require(parallel)
        no_cores <- detectCores() - 1

        # Initiate cluster
        cl1 <- makeCluster(no_cores)
        parApply(
          cl = cl1,
          rndSignals,
          MARGIN = 1,
          FUN = function(x) {
            aux_rnd_dist <-
              wsd(
                signal1 = x[, 1],
                signal2 = x[, 2],
                dt = dt,
                scaleparam = scaleparam,
                delta_t = delta_t,
                windowrad = windowrad,
                rdist = rdist,
                normalize = normalize,
                wname = wname, wparam = wparam, waverad = waverad,
                border_effects = border_effects,
                mc_nrand = 0,
                commutative = commutative,
                compensation = compensation,
                wscnoise = wscnoise,
                energy_density = energy_density,
                parallel = parallel,
                makefigure = FALSE)
          }
        ) -> kk
      stopCluster(cl1)

      } else {

        apply(
          rndSignals,
          MARGIN = 1,
          FUN = function(x) {
            aux_rnd_dist <-
              wsd(
                signal1 = x[, 1],
                signal2 = x[, 2],
                dt = dt,
                scaleparam = scaleparam,
                delta_t = delta_t,
                windowrad = windowrad,
                rdist = rdist,
                normalize = normalize,
                wname = wname, wparam = wparam, waverad = waverad,
                border_effects = border_effects,
                mc_nrand = 0,
                commutative = commutative,
                compensation = compensation,
                wscnoise = wscnoise,
                energy_density = energy_density,
                parallel = parallel,
                makefigure = FALSE)
          }
        ) -> kk
      }

      lapply(
        kk,
        function(x) {
          abind(x$wsd, along = 3)
        }
      ) -> kk2
      do.call(abind, kk2) -> distrnd

      # Significance

      print("Computing significance contours...")
      kk95 <- matrix(0, nrow = nwsc, ncol = ns)
      kk05 <- matrix(0, nrow = nwsc, ncol = ns)

      Aux <- -log2(distrnd)

      kk95 <- apply(
        X = Aux,
        MARGIN = c(1, 2),
        FUN = function(x)
          quantile(x, .95, na.rm = TRUE)
      )
      kk05 <-
        apply(
          X = Aux,
          MARGIN = c(1, 2),
          FUN = function(x)
            quantile(x, .05, na.rm = TRUE)
        )

      wsdsig95 <- ((-log2(wsd) - kk95) > 0)
      wsdsig05 <- ((-log2(wsd) - kk05) > 0)

    } else {
      wsdsig95 <- NULL
      wsdsig05 <- NULL
    }

    ### COI WSC

    coi_wsc <- wsc1$coi_maxscale

    ### Make figure

    tcentraldt <- wsc1$tcentral

    if (makefigure) {

      if (figureperiod) {
        Y <- fourierfactor * scalesdt
        coi <- fourierfactor * coi_wsc
        if (is.null(ylab)) ylab <- "Period"
      } else {
        Y <- scalesdt
        coi <- coi_wsc
        if (is.null(ylab)) ylab <- "Scale"
      }

      if (is.null(time_values)) {
        X <- tcentraldt
      } else {
        if (length(time_values) != nt) {
          warning("Invalid length of time_values vector. Changing to default.")
          X <- tcentraldt
        } else {
          X <- time_values[floor(tcentraldt / dt)]
        }
      }

      wavPlot(
        Z = -log2(wsd),
        X = X,
        Y = Y,
        Ylog = TRUE,
        coi = coi,
        rdist = rdist,
        sig95 = wsdsig95,
        sig05 = wsdsig05,
        Xname = xlab,
        Yname = ylab,
        Zname = main
      )

    }

    return(list(
      wsd = wsd,
      tcentral = tcentraldt,
      scales = scalesdt,
      windowrad = windowrad,
      rdist = rdist,
      signif95 = wsdsig95,
      signif05 = wsdsig05,
      fourierfactor = fourierfactor,
      coi_maxscale = coi_wsc
    ))

  }



