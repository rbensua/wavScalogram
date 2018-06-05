#' @title Windowed Scalogram Difference
#'
#' @description This function computes the Windowed Scalogram Diference of two signals. The definition and details
#' can be found in (Bolos et. al, AMC 2017)
#'
#' @usage wsd(signal1, signal2, scaleparam = NULL, delta_t = NULL, windowrad = NULL, rdist = NULL, normalize = TRUE,
#' wname = c("MORLET", "HAAR"), border.effects = c("BE","INNER", "PER","SYM"), mc_nrand = NULL, compensated = TRUE,
#' commutative = TRUE, energy_density = TRUE, parallel = FALSE, makefigure = TRUE)
#'
#' @param signal1 A vector containing the first signal.
#' @param signal2 A vector containing the second signal (its length should be equal to that of \code{signal1})
#' @param scaleparam A vector of three elements with the minimum scale, the maximum scale and the number
#'          of suboctaves per octave for constructing power 2 scales (following Torrence and Compo 1998),
#'          measured in units of time.
#' @param delta_t Numeric. Increment of time for the construction of windows central times.
#' @param windowrad Numeric. Time radius for the windows.
#' @param rdist Numeric. Log-scale radius for the windows measured in suboctaves.
#' @param normalize Logical. Set to TRUE if the signals use different units.
#' @param wname A string, equal to 'MORLET' or 'HAAR'.
#' @param boder.effects String, equal to "BE", "INNER", "PER" or "SYM",
#' which indicates how to manage the border effects which arise usually when a convolution is
#' performed on finite-lenght signals.
#' \itemize{
#' \item "BE": With border effects, padding time series with zeroes.
#' \item "INNER": Normalized inner scalogram with security margin adapted for each different scale.
#' \item "PER": With border effects, using boundary wavelets (periodization of the
#' original time series). Recommended for "short" series.
#' \item "SYM": With border effects, using a symmetric catenation of the original time
#' series. Also recommended for "short" series.
#' }
#' @param compensated Logical. For preventing numerical errors in the computation of relative distances.
#' @param mc_nrand Integer. Number of Monte Carlo simulations to be performed in order to determine the 95\% significance contours.
#' @param commutative Logical. If TRUE (default) the commutative windowed scalogram difference. Otherwise a non-commutative (but simpler) version is computed
#' (see AMC 2017).
#' @param energy_density Logical. Divide the scalograms by sqrt(scales) for convert them into
#' energy density (defaults to TRUE).
#' @param parallel Logical. If TRUE (default) uses function \code{parApply} from package \code{parallel} for the Monte Carlo simulations.
#' When FALSE is uses the normal \code{apply} function.
#' @param makefigure Logical. Plots a figure with the WSD
#'
#' @importFrom  colorRamps matlab.like
#' @importFrom  fields image.plot
#' @importFrom  parallel parApply
#' @importFrom  abind abind
#' @importFrom reshape melt
#' @importFrom scales trans_new log_breaks
#' @import ggplot2
#' @import dplyr
#' @return A list of a matrix of size \code{length(scales)}x\code{length(tcentral)} containing
#' the values of the windowed scalograms at each scale and at each time window, a vector
#' of length \code{length(scales)} containing the scales, and a vector of length \code{length(tcentral)}
#' containing the values of the maxim scale from which there are border effects for the respective
#' windowed scalogram.
#'
#' @examples
#' nrand <- 10  # MonteCarlo
#' windowrad <- 50
#' rdist <- 5 # Scale radius for distance (height)
#' nt <- 1500
#' t <- 1:nt
#' sd_noise <-  0.2 #% In the paper Figure 1, sd_noise = 1
#' signal1 <- rnorm(n = nt, mean = 0, sd = sd_noise)+sin((1:nt)/10)
#' signal2 <- rnorm(n = nt, mean = 0, sd = sd_noise)+sin((1:nt)/10)
#' signal2[500:1000] = signal2[500:1000]+sin((500:1000)/2)
#'
#' waverad <- 3 # For Morlet
#' s0 <- 2
#' smax <- (nt-1-2*windowrad)/(2*waverad)
#' Dj <- 12
#' scaleparam <- c(s0, smax, Dj)
#' \dontrun{
#' wsd <-wsd(signal1 = signal1, signal2 = signal2, scaleparam = scaleparam, delta_t = 10,
#'        windowrad = windowrad, rdist = rdist, mc_nrand = nrand,
#'        wname = "MORLET", makefigure =  TRUE)
#' }
#' @section References:
#' C. Torrence, G. P. Compo. A practical guide to wavelet analysis. B. Am. Meteorol. Soc.
#' 79 (1998), 61–78.
#'
#' V. J. Bolós, R. Benítez, R. Ferrer, R. Jammazi. The windowed scalogram difference: a novel wavelet tool for comparing time series.
#' Appl. Math. Comput., 312 (2017), 49-65.
#'
#' @export
#'

wsd_old <-
  function(signal1,
           signal2,
           scaleparam = NULL,
           delta_t = NULL,
           windowrad = NULL,
           rdist = NULL,
           normalize = FALSE,
           wname = c("MORLET", "HAAR"),
           border.effects = c("BE","INNER", "PER","SYM"),
           mc_nrand = NULL,
           compensated = FALSE,
           commutative = TRUE,
           energy_density = TRUE,
           parallel = FALSE,
           makefigure = TRUE) {

    require(abind)

    wname <- toupper(wname)
    wname <- match.arg(wname)
    if (wname == "MORLET") {
      waverad <- 3
    } else if (wname == "HAAR") {
      waverad <- 0.5 # Antes ponía 3/4, por qué??
    }

    border.effects <- toupper(border.effects)
    border.effects <- match.arg(border.effects)

    if (is.null(mc_nrand)) {
      mc_nrand <- 10
    }

    nt <- length(signal1)
    if (nt != length(signal2)) {
      stop("Signals must have the same length.")
    }

    if (is.null(delta_t)) {
      delta_t <- ceiling(nt/256)
    }

    if (is.null(windowrad)) {
      windowrad <- ceiling(nt/20)
    } else if (windowrad > floor((nt-1)/2)) {
      windowrad <- floor((nt-1)/2)
      print("Windowrad re-adjusted.")
    }

    if (is.null(scaleparam)) {
      scmax <- floor((nt-2*windowrad)/(2*waverad))
      scaleparam <- c(2, scmax, ceiling(256/(log2(scmax)-1)))
    }

    scales <- pow2scales(scaleparam)
    ns <- length(scales)

    if (is.null(rdist)) {
      rdist <- floor(ns/20)
    }

    ### Scalograms

    scalog1 <-
      scalogram_old(signal = signal1,
                scales = scales, powerscales = FALSE,
                border.effects = border.effects,
                wname = wname)
    scalog2 <-
      scalogram_old(signal = signal2,
                scales = scales, powerscales = FALSE,
                border.effects = border.effects,
                wname = wname)

    norm1 <- sqrt(sum(scalog1$scalog ^ 2))
    if (normalize) {
      norm2 <- sqrt(sum(scalog2$scalog ^ 2))
    } else{
      norm2 <- norm1
    }
    scalog1 <- scalog1$scalog / norm1
    scalog2 <- scalog2$scalog / norm2

    ### Windowed scalograms

    wsc1 <-
      windowed_scalogram_old2(
        signal = signal1,
        scales = scales, powerscales = FALSE,
        border.effects = border.effects,
        windowrad = windowrad,
        delta_t = delta_t,
        wname = wname
      )
    wsc2 <-
      windowed_scalogram_old2(
        signal = signal2,
        scales = scales, powerscales = FALSE,
        border.effects = border.effects,
        windowrad = windowrad,
        delta_t = delta_t,
        wname = wname
      )
    if (energy_density) {
      wsc1$wsc <- wsc1$wsc / sqrt(scales)
      wsc2$wsc <- wsc2$wsc / sqrt(scales)
    }
    nwsc <- length(wsc1$tcentral)

    # Normalization
    wsc1$wsc <- wsc1$wsc / norm1
    wsc2$wsc <- wsc2$wsc / norm2

    # Compensation

    if (compensated) {
      fc <- 100 / min(mean(scalog1), mean(scalog2))
      wsc1$wsc <- 1 + fc * wsc1$wsc
      wsc2$wsc <- 1 + fc * wsc2$wsc
    }

    # Relative distances between scalograms
    wsd <- matrix(0, nrow = ns, ncol = nwsc)

    kmin <- pmax(1, (1:ns) - rdist)
    kmax <- pmin(ns, (1:ns) + rdist)
    idxlist <- mapply(
      kmin,
      kmax,
      FUN = function(i, j)
        i:j
    )
    if (commutative) {
      A <- ((wsc1$wsc ^ 2 - wsc2$wsc ^ 2) / (2 * wsc2$wsc * wsc1$wsc)) ^ 2
    } else{
      A <- ((wsc1$wsc ^ 2 - wsc2$wsc ^ 2) / wsc1$wsc) ^ 2
    }
    wsd <- apply(A, MARGIN = 2, function(r) {
      lapply(
        idxlist,
        FUN = function(x) {
          sum(r[x]) * (2 * rdist + 1) / (diff(range(x)) + 1)
        }
      )
    })
    wsd <- matrix(unlist(wsd), nrow = nwsc, byrow = TRUE)
    wsd <- t(sqrt(wsd))

    ## Montecarlo
    if (mc_nrand > 0) {
      mean1 <- mean(signal1)
      mean2 <- mean(signal2)
      sd1 <- sd(signal1)
      sd2 <- sd(signal2)
      distrnd <- array(0, dim = c(ns, nwsc, mc_nrand))

      print("Montecarlo...")
      rndSignals1 <- matrix(rnorm(mc_nrand * nt, mean = mean1, sd = sd1), nrow = mc_nrand, ncol = nt)
      rndSignals2 <- matrix(rnorm(mc_nrand * nt, mean = mean2, sd = sd2), nrow = mc_nrand, ncol = nt)
      rndSignals <- array(c(rndSignals1, rndSignals2), dim = c(mc_nrand, nt, 2))

      if(parallel){
        require(parallel)
        no_cores <- detectCores() - 1

        # Initiate cluster
        cl1 <- makeCluster(no_cores)
        parApply(cl = cl1, rndSignals, MARGIN = 1, FUN = function(x){

          aux_rnd_dist <-
            wsd_old(
              signal1 = x[,1],
              signal2 = x[,2],
              scaleparam = scaleparam,
              delta_t = delta_t,
              windowrad = windowrad,
              rdist = rdist,
              normalize = normalize,
              wname = wname,
              border.effects = border.effects,
              mc_nrand = 0,
              compensated = compensated,
              commutative = commutative,
              energy_density = energy_density,
              makefigure = FALSE,
              parallel = parallel
            )

        }) -> kk
        stopCluster(cl1)
      }else{
      apply(rndSignals, MARGIN = 1, FUN = function(x){

        aux_rnd_dist <-
          wsd_old(
            signal1 = x[,1],
            signal2 = x[,2],
            scaleparam = scaleparam,
            delta_t = delta_t,
            windowrad = windowrad,
            rdist = rdist,
            normalize = normalize,
            wname = wname,
            border.effects = border.effects,
            mc_nrand = 0,
            compensated = compensated,
            commutative = commutative,
            energy_density = energy_density,
            makefigure = FALSE,
            parallel = parallel
          )

      }) -> kk
      }
      lapply(kk, function(x){abind(x$wsd, along = 3)}) -> kk2
      do.call(abind,kk2) -> distrnd


      # Significance
      print("Computing significance contours...")
      kk95 <- matrix(0, nrow = ns, ncol = nwsc)
      kk05 <- matrix(0, nrow = ns, ncol = nwsc)

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

      wsdsig95 <- -log2(wsd) - kk95
      wsdsig05 <- -log2(wsd) - kk05
      sig95 = t(wsdsig95)
      sig05 = t(wsdsig05)
    } else {
      wsdsig95 <- NULL
      wsdsig05 <- NULL
      sig95 = NULL
      sig05 = NULL
    }
    # COI WSC
    if (border.effects %in% c("BE","PER","SYM")) {
      coi_wsc <- wsc1$coi_maxscale
    } else{
      coi_wsc <- NULL
    }
    if (makefigure){
      #require(colorRamps)
      #require(fields)

    wavPlot(
      Z = -log2(t(wsd)),
      X = wsc1$tcentral,
      Y = scales,
      Ylog = TRUE,
      coi = coi_wsc,
      sig95 = sig95,
      sig05 = sig05,
      Xname = "Time",
      Yname = "Scale",
      Zname = "WSD"
    )
      #wsdPlot(t = wsc1$tcentral, scales = scales, wsd = wsd, coi = coi_wsc,
      #       sig95 = wsdsig95, sig05 = wsdsig05, rdist = rdist)
    }
    return(list(
      wsd = wsd,
      t = wsc1$tcentral,
      scales = scales,
      coi = wsc1$coi_maxscale,
      signif95 = wsdsig95,
      signif05 = wsdsig05
    ))


  }


wsdPlot <- function(wsd, t, scales, sig95, sig05, coi, rdist,...) {
#  require(colorRamps)
#  require(fields)

  COLORS <- matlab.like(50)
  nesc <- length(scales)
  nwsc <- length(t)
  XY <- expand.grid(scales,t)
  Z <- melt(wsd)
  DF <- data.frame(Scales = XY$Var1, Time = XY$Var2, WSD = Z$value)

  middle <- min(-log2(DF$WSD))
  mask <- matrix(ncol = ncol(wsd), nrow = nrow(wsd))
  coimask <- mask
  mask[1:(1+rdist),] <- middle
  mask[(nesc-rdist):nesc,] <- middle
  DF$mask <- melt(mask)$value
  DF %>% ggplot(aes(x = Time, y = Scales)) + geom_raster(aes(fill = -log2(WSD))) +
    scale_y_continuous(expand = c(0,0), trans = reverselog_trans(2)) +
    geom_raster(aes(fill = mask), alpha = 0.3) +
    scale_fill_gradientn(colours = COLORS, na.value = 'transparent') +
    scale_x_continuous(expand = c(0,0)) +
    geom_hline(yintercept = scales[length(scales) - rdist]) +
    geom_hline(yintercept = scales[1 + rdist])  -> p

  if(!is.null(coi)){
  for(j in 1:dim(mask)[2]){
    coimask[,j] <- ifelse(scales > coi[j], middle, NA)
  }
  DF$coimask  <- melt(coimask)$value
  coi[coi < min(scales)] <- min(scales)
  coi[coi > max(scales)] <- max(scales)
  p <- p + geom_raster(aes(fill = coimask), alpha = 0.3) +
    scale_fill_gradientn(colours = COLORS, na.value = 'transparent') +
    geom_line(aes(x = time, y = scales),
                     data = data.frame(time = t,
                                       scales = as.numeric(coi)))
  }


  if(!is.null(sig95)){
      sig95 <- data.frame(Scales = DF$Scales , Time = DF$Time, sig95 = melt(sig95)$value)

  sig05 <- data.frame(Scales = DF$Scales , Time = DF$Time, sig05 = melt(sig05)$value)
  p <- p + geom_contour(data = sig95, aes(x = Time, y = Scales, z = sig95),
               breaks = 0, col = "black", lwd = 1) +
    geom_contour(data = sig05, aes(x = Time, y = Scales, z = sig05),
                 breaks = 0, col = "white", lwd = 1)

  }
  print(p)
  # Window margins (rdist)







#
#
#   image.plot(
#     t,
#     log2(scales),
#     -log2(t(wsd[nesc:1, ])),
#     col = COLORS,
#     yaxt = "n",
#     xlab = "Time",
#     ylab = "Scale",
#     ylim = c(min(log2(scales)), max(log2(scales))),
#     main = expression(-log[2](WSD))
#   )
#
#   # Window margins (rdist)
#   abline(h = log2(scales[nesc - rdist]))
#   abline(h = log2(scales[rdist+1]))
#   rdist_mask <- matrix(NA, nrow = nwsc, ncol = nesc)
#   rdist_mask[,1:(1+rdist)] <- 1
#   rdist_mask[,(nesc-rdist):nesc] <- 1
#   pal <- colorRampPalette(c(rgb(1, 1, 1), rgb(0, 0, 0)))
#   COLORS3 <- add.alpha(pal(20), 0.5)
#   image(t,
#         log2(scales),
#         rdist_mask[, nesc:1],
#         col = COLORS3,
#         add = T)
#
#   # COI
#   if (!is.null(coi)) {
#     coi_matrix <- matrix(NA, nrow = nwsc, ncol = nesc)
#     for (i in 1:nwsc) {
#       coi_matrix[i, scales > coi[i]] <- 1
#     }
#     pal <- colorRampPalette(c(rgb(1, 1, 1), rgb(0, 0, 0)))
#     COLORS2 <- add.alpha(pal(20), 0.5)
#     image(t,
#           log2(scales),
#           coi_matrix[, nesc:1],
#           col = COLORS2,
#           add = T)
#
#     log2coi <- (max(log2(scales)) - log2(coi) + min(log2(scales)))
#     segments(
#       x0 = t[-nwsc],
#       x1 = t[-1],
#       y0 = log2coi[-nwsc],
#       y1 = log2coi[-1],
#       ylim = c(max(log2coi), min(log2(scales)))
#     )
#   }
#
#
#
#   if (!is.null(sig95)) {
#     contour(
#       t,
#       log2(scales),
#       t(sig95[nesc:1, ]),
#       add = T,
#       levels = c(0),
#       lwd = 3,
#       labcex = 0.01,
#       labels = ""
#     )
#     contour(
#       t,
#       log2(scales),
#       t(sig05[nesc:1, ]),
#       add = T,
#       levels = c(0),
#       col = "white",
#       lwd = 3,
#       labels = "",
#       labcex = 0.01
#     )
#   }
#
#   axis(
#     side = 2,
#     at = seq(from = log2(max(scales)), to = floor(log2(min(scales))) , by = -1),
#     labels = 2 ^ (seq(
#       to = floor(log2(max(scales))), from = floor(log2(min(scales)))
#     ))
#   )
}


add.alpha <- function(COLORS, ALPHA){
  if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
  RGB <- col2rgb(COLORS, alpha=TRUE)
  RGB[4,] <- round(RGB[4,]*ALPHA)
  NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
  return(NEW.COLORS)
}





