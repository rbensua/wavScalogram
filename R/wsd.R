#' @title Windowed Scalogram Difference
#'
#' @description This function computes the Windowed Scalogram Difference of two signals.
#' The definition and details can be found in (Bolós et al. 2017).
#'
#' @usage wsd(signal1,
#'            signal2,
#'            scaleparam = NULL,
#'            delta_t = NULL,
#'            windowrad = NULL,
#'            rdist = NULL,
#'            normalize = FALSE,
#'            wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
#'            wparam = NULL,
#'            border_effects = c("BE", "INNER", "PER", "SYM"),
#'            mc_nrand = 0,
#'            commutative = TRUE,
#'            wscnoise = 0.02,
#'            compensation = 0,
#'            energy_density = TRUE,
#'            parallel = FALSE,
#'            makefigure = TRUE)
#'
#' @param signal1 A vector containing the first signal.
#' @param signal2 A vector containing the second signal (its length should be equal to
#' that of \code{signal1}).
#' @param scaleparam A vector of three elements with the minimum scale, the maximum scale
#' and the number of suboctaves per octave for constructing power 2 scales (following
#' Torrence and Compo 1998), measured in units of time.
#' @param delta_t Numeric. Increment of time for the construction of windows central
#' times.
#' @param windowrad Numeric. Time radius for the windows.
#' @param rdist Numeric. Log-scale radius for the windows measured in suboctaves.
#' @param normalize Logical. Set to TRUE if the signals use different units.
#' @param wname A string, equal to "MORLET", "DOG", "PAUL", "HAAR" or "HAAR2". The
#' difference between "HAAR" and "HAAR2" is that "HAAR2" is more accurate but slower.
#' @param wparam Numeric. Parameters of the corresponding wavelet.
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
#' @param compensation Numeric in \eqn{[0,1]}. It is an alternative to \code{wscnoise} for
#' preventing numerical errors or non-relevant high relative differences when scalogram
#' values are close to zero (see Bolós et al. 2017).
#' @param energy_density Logical. Divide the scalograms by \eqn{\sqrt(scales)} for convert
#' them into energy density (defaults to TRUE). Note that it does not affect the results
#' if \code{wscnoise} \eqn{= 0}.
#' @param parallel Logical. If TRUE (default) uses function \code{parApply} from package
#' \code{parallel} for the Montecarlo simulations. When FALSE is uses the normal
#' \code{apply} function.
#' @param makefigure Logical. Plots a figure with the WSD.
#'
#' @importFrom parallel parApply detectCores makeCluster stopCluster
#' @importFrom abind abind
#'
#' @return A list with the following fields:
#'
#' \code{wsd}: A matrix of size \code{length(tcentral)}x\code{length(scales)} containing
#' the values of the windowed scalogram differences at each scale and at each time window.
#'
#' \code{t}: The vector of central times used in the computations of the windowed
#' scalograms.
#'
#' \code{scales}: The vector of scales used in the computations.
#'
#' \code{coi}: A vector of length \code{length(tcentral)} containing the values of the
#' maximum scale from which there are border effects for the windowed scalograms.
#'
#' \code{rdist}: The log-scale radius for the windows measured in suboctaves.
#'
#' \code{signif95}: A logical matrix of size
#' \code{length(tcentral)}x\code{length(scales)}. If TRUE, the corresponding point of the
#' \code{wsd} matrix is in the 95\% significance.
#'
#' \code{signif05}: A logical matrix of size
#' \code{length(tcentral)}x\code{length(scales)}. If TRUE, the corresponding point of the
#' \code{wsd} matrix is in the 5\% significance.
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
           scaleparam = NULL,
           delta_t = NULL,
           windowrad = NULL,
           rdist = NULL,
           normalize = FALSE,
           wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"),
           wparam = NULL,
           border_effects = c("BE", "INNER", "PER", "SYM"),
           mc_nrand = 0,
           commutative = TRUE,
           wscnoise = 0.02,
           compensation = 0,
           energy_density = TRUE,
           parallel = FALSE,
           makefigure = TRUE) {

  #  require(abind)

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

    if (is.null(scaleparam)) {
      scmax <- floor((nt - 2 * windowrad) / (2 * waverad))
      scaleparam <- c(2, scmax, ceiling(256 / (log2(scmax) - 1)))
    }

    scales <- pow2scales(scaleparam)
    ns <- length(scales)

    if (is.null(rdist)) {
      rdist <- floor(ns / 20)
    }

    ### Scalograms

    scalog1 <-
      scalogram(signal = signal1,
                scales = scales, powerscales = FALSE,
                border_effects = border_effects,
                energy_density = energy_density,
                wname = wname, wparam = wparam)
    scalog2 <-
      scalogram(signal = signal2,
                scales = scales, powerscales = FALSE,
                border_effects = border_effects,
                energy_density = energy_density,
                wname = wname, wparam = wparam)

    ### Windowed scalograms

    wsc1 <-
      windowed_scalogram(
        signal = signal1,
        scales = scales, powerscales = FALSE,
        border_effects = border_effects,
        energy_density = energy_density,
        windowrad = windowrad,
        delta_t = delta_t,
        wname = wname, wparam = wparam
      )
    wsc2 <-
      windowed_scalogram(
        signal = signal2,
        scales = scales, powerscales = FALSE,
        border_effects = border_effects,
        energy_density = energy_density,
        windowrad = windowrad,
        delta_t = delta_t,
        wname = wname, wparam = wparam
      )

    nwsc <- length(wsc1$tcentral)

    # Normalization

    norm1 <- sqrt(sum(scalog1$scalog ^ 2))
    if (normalize) {
      norm2 <- sqrt(sum(scalog2$scalog ^ 2))
    } else {
      norm2 <- norm1
    }
    scalog1$scalog <- scalog1$scalog / norm1
    scalog2$scalog <- scalog2$scalog / norm2
    wsc1$wsc <- wsc1$wsc / norm1
    wsc2$wsc <- wsc2$wsc / norm2

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
                scaleparam = scaleparam,
                delta_t = delta_t,
                windowrad = windowrad,
                rdist = rdist,
                normalize = normalize,
                wname = wname, wparam = wparam,
                border_effects = border_effects,
                mc_nrand = 0,
                commutative = commutative,
                compensation = compensation,
                wscnoise = wscnoise,
                energy_density = energy_density,
                makefigure = FALSE,
                parallel = parallel
              )
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
                scaleparam = scaleparam,
                delta_t = delta_t,
                windowrad = windowrad,
                rdist = rdist,
                normalize = normalize,
                wname = wname, wparam = wparam,
                border_effects = border_effects,
                mc_nrand = 0,
                commutative = commutative,
                compensation = compensation,
                wscnoise = wscnoise,
                energy_density = energy_density,
                makefigure = FALSE,
                parallel = parallel
              )
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

    if (makefigure) {
      wavPlot(
        Z = -log2(wsd),
        X = wsc1$tcentral,
        Y = scales,
        Ylog = TRUE,
        coi = coi_wsc,
        rdist = rdist,
        sig95 = wsdsig95,
        sig05 = wsdsig05,
        Xname = "Time",
        Yname = "Scale",
        Zname = "-log2(WSD)"
      )
    }

    return(list(
      wsd = wsd,
      t = wsc1$tcentral,
      scales = scales,
      coi = wsc1$coi_maxscale,
      rdist = rdist,
      signif95 = wsdsig95,
      signif05 = wsdsig05
    ))

  }



