diagnostic.summary <- function(codaMCMClist, HDImass = 0.95, gelman_diag = TRUE,
                               gelman_diag_multivariate = TRUE) {
  parameterNames = varnames(codaMCMClist)
  mcmcMat = as.matrix(codaMCMClist,chains=TRUE)
  summaryInfo = NULL
  for ( parName in parameterNames ) {
    summaryInfo = rbind( summaryInfo , summarizePost( mcmcMat[,parName], credMass = HDImass ) )
    thisRowName = parName
    rownames(summaryInfo)[NROW(summaryInfo)] = thisRowName
  }
  summaryInfo_df <- as.data.frame(summaryInfo)
  if(gelman_diag == TRUE) {
    psrf_df <- as.data.frame((gelman.diag(codaMCMClist, multivariate = gelman_diag_multivariate))$psrf)
    colnames(psrf_df) <- c("PSRF Point est.", "PSRF Upper C.I.")
    diagnostic_summary <- cbind(psrf_df, summaryInfo_df)
  }  else {
    diagnostic_summary <- summaryInfo_df
  }
}

# simplified version of a similar function in Kruschke 2015
#' @rdname amtl_bayes_helper
#' @export
summarizePost = function( paramSampleVec , credMass=0.95 ) {
  meanParam = mean( paramSampleVec )
  medianParam = median( paramSampleVec )
  dres = density( paramSampleVec )
  modeParam = dres$x[which.max(dres$y)]
  mcmcEffSz = round( effectiveSize( paramSampleVec ) , 1 )
  names(mcmcEffSz) = NULL
  MCSE = sd(paramSampleVec)/sqrt(mcmcEffSz)
  hdiLim = HDIofMCMC( paramSampleVec , credMass=credMass )
  return( c( Mean=meanParam , Median=medianParam , Mode=modeParam ,
             ESS=mcmcEffSz , MCSE = MCSE,
             HDImass=credMass , HDIlow=hdiLim[1] , HDIhigh=hdiLim[2]) )
}

# simplified version of a similar function in Kruschke 2015
#' @rdname amtl_bayes_helper
#' @export
HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

# from Kruschke 2015
#' @rdname amtl_bayes_helper
#' @export
gammaShRaFromMeanSD = function( mean , sd ) {
  if ( mean <=0 ) stop("mean must be > 0")
  if ( sd <=0 ) stop("sd must be > 0")
  shape = mean^2/sd^2
  rate = mean/sd^2
  return( list( shape=shape , rate=rate ) )
}

gammaShRaFromModeSD = function( mode , sd ) {
  if ( mode <=0 ) stop("mode must be > 0")
  if ( sd <=0 ) stop("sd must be > 0")
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

# from stackoverflow, helps to find initial values of parameters alpha and beta for
# estimation via MLE
#' @rdname amtl_bayes_helper
#' @export
beta_mom <- function(x) {

  m_x <- mean(x, na.rm = TRUE)
  s_x <- sd(x, na.rm = TRUE)

  alpha <- m_x*((m_x*(1 - m_x)/s_x^2) - 1)
  beta <- (1 - m_x)*((m_x*(1 - m_x)/s_x^2) - 1)

  return(list(alpha = alpha, beta = beta))
}


# computing age from logistic regression if slope and intercept of tooth loss are known
#' @rdname amtl_bayes_helper
#' @export
tooth.loss.recon <- function(x, amtl, tp, amtl_intercept, amtl_m, minimum_age = 15, trunc = FALSE) {
  amtl <- x[,"amtl"]
  tp <- x[,"tp"]
  tp_perc <- amtl/tp * 100
  tooth_loss_recon <- data.frame(age_approx = NA)
  for (i in 1:length(tp_perc)) {
    if (tp_perc[i] == 0) {
      age_approx <- ((log(1.5/101.5) - amtl_intercept) / amtl_m)
    } else if (tp_perc[i] == 100) {
      age_approx <- ((log(101.5/1.5) - amtl_intercept) / amtl_m)
    }
    else {
      amtl_log <- log(tp_perc[i]/(100 - tp_perc[i]))
      age_approx <- ((amtl_log  - amtl_intercept) / amtl_m)
    }
    if (trunc == TRUE) {
      if (age_approx < 15) {
        age_approx <- 15
      } else if (age_approx > 99) {
        age_approx <- 99
      } else {
        age_approx <- age_approx
      }
    }
    tooth_loss_recon <- rbind(tooth_loss_recon, age_approx)
  }
  tooth_loss_recon <- tooth_loss_recon[-1,]
  output <- cbind(x, "age_recon" = tooth_loss_recon)
  return(output)
}

# summarizing and comparing age estimates to true age
#' @rdname amtl_bayes_helper
#' @export
age.estimate.summarize <- function(x,
                                   known_age,
                                   chosen_mean,
                                   age_identifier = "age.s",
                                   HDIhigh = "HDIhigh",
                                   HDIlow = "HDIlow",
                                   HDImass = 0.95,
                                   gelman_diag = TRUE,
                                   gelman_diag_multivariate = F) {
  age_identifier_grep <- ifelse(age_identifier == "age.s", "^age.s\\[", "^age.s_c")
  x_mcmc_list <- as.mcmc.list(x)
  x_diag <- diagnostic.summary(x_mcmc_list,
                               HDImass = HDImass,
                               gelman_diag = gelman_diag,
                               gelman_diag_multivariate = gelman_diag_multivariate)
  x_diag_red <- x_diag[grep(age_identifier_grep,rownames(x_diag)),]
  x_mcmcMat = as.matrix(x_mcmc_list, chains=TRUE)

  estimated_age <- x_diag_red[,chosen_mean]
  Bias_model  <-  lm((known_age - estimated_age) ~ known_age)
  Bias <- Bias_model$coefficients[2]

  corrPearson <- cor.test(estimated_age,known_age, method="pearson")

  ages <- x_mcmcMat[,grep(age_identifier_grep,colnames(x_mcmcMat))]
  tmnlp_res <- tmnlp(known_age, ages)
  crps_res <- mean(scoringRules::crps_sample(known_age, t(ages)) )

  age_estimation <- data.frame(Mean_estimated = mean(estimated_age),
                               Offset = mean(estimated_age - known_age),
                               corrPearson = corrPearson$estimate,
                               corr_p = corrPearson$p.value,
                               Bias = Bias_model$coefficients[2],
                               Inaccuracy = mean(abs(estimated_age - known_age)),
                               RMSE = Metrics::rmse(known_age, estimated_age),
                               TMNLP = tmnlp_res,
                               CRPS = crps_res)
  HDIhigh <- x_diag_red[,HDIhigh]
  HDIlow <- x_diag_red[,HDIlow]
  age_estimation$Coverage <- sum(ifelse(known_age >= HDIlow & known_age <=
                                          HDIhigh, 1, 0)) / length(known_age) * 100
  age_estimation$HDI_Diff_median <- stats::quantile(HDIhigh - HDIlow, probs = c(0.5))
  age_estimation$HDI_Diff_quant_025 <- stats::quantile(HDIhigh - HDIlow, probs = c( 0.025))
  age_estimation$HDI_Diff_quant_975 <- stats::quantile(HDIhigh - HDIlow, probs = c(0.975))
  return(age_estimation)
}

# generating "the mean value of the negative log posterior evaluated at the known ages of the test observations"
#' @rdname amtl_bayes_helper
#' @export
# idea deriving from https://github.com/ElaineYChu/fs_mcp_us/blob/f5a54a74959442fd68647c7cba983cac62ed23cc/mcp_model_performance.R#L162C12-L162C46
tmnlp <- function(x, mcmcMat) { # x = true age, mcmcMat = output from coda MCMC
  val_vec <- c()
  x_length <- length(x)
  for (i in 1: x_length) {
    age_dens <- density(mcmcMat[i,])
    age_pred <- median(age_dens$y[round(age_dens$x) == x[i]])
    val_vec <- c(val_vec, log(age_pred))
  }
  tmnlp <- round(-sum(val_vec, na.rm=T) / x_length, 4)
  return(tmnlp)
}

# Continuous Ranked Probability Score, code adapted from
# https://towardsdatascience.com/crps-a-scoring-function-for-bayesian-machine-learning-models-dd55a7a337a8
# https://github.com/itamarfaran/public-sandbox/blob/master/crps.ipynb
#' @rdname amtl_bayes_helper
#' @export
crps_nrg <- function(y_true, y_pred) {
  num_samples <- nrow(y_pred)
  absolute_error <- colMeans(abs(sweep(y_pred, 2, y_true)))
  y_pred <- apply(y_pred, 2, sort)
  diff <- y_pred[-1, ] - y_pred[-num_samples, ]
  weight <- matrix(NA, nrow = num_samples - 1, ncol = ncol(y_pred))
  for (i in 1:(num_samples - 1)) {
    weight[i, ] <- i * (num_samples - i)
  }
  per_obs_crps <- absolute_error - colSums(diff * weight) / num_samples^2
  return(mean(per_obs_crps))
}

# Plot function from Kruschke 2015
#' @rdname amtl_bayes_helper
#' @export
DbdaDensPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL , credMass=0.95 ) {
  nChain = length(codaObject) # or nchain(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  hdiLims = NULL
  credMass_label = paste0(credMass, ' HDI')
  for ( cIdx in 1:nChain ) {
    densInfo = density(codaObject[,c(parName)][[cIdx]])
    xMat = cbind(xMat,densInfo$x)
    yMat = cbind(yMat,densInfo$y)
    hdiLims = cbind(hdiLims,HDIofMCMC(codaObject[,c(parName)][[cIdx]] , credMass=credMass))
  }
  matplot( xMat , yMat , type="l" , col=plColors ,
           main="" , xlab=parName , ylab="Density" )
  abline(h=0)
  points( hdiLims[1,] , rep(0,nChain) , col=plColors , pch="|" )
  points( hdiLims[2,] , rep(0,nChain) , col=plColors , pch="|" )
  text( mean(hdiLims) , 0 , credMass_label , adj=c(0.5,-0.2) )
}

# Function for the creation of summed age densities, with sampling to reduce size
#' @rdname amtl_bayes_helper
#' @export
coda.object.sum <- function(codaObject, sampled = TRUE) {
  nChain = length(codaObject) # or nchain(codaObject)
  parameterNames = varnames(codaObject) # get all parameter names for reference
  parameterNames <- grep("^age.s", parameterNames,value=TRUE)
  coda_object_simplified <- NULL
  coda_object_chains <- NULL
  for ( parName in parameterNames ) {
    onecolumn_chains <- NULL
    for ( cIdx in 1:nChain ) {
      onecolumn <- as.data.frame(codaObject[,c(parName)][[cIdx]])
      onecolumn_chains <- rbind(onecolumn_chains, onecolumn)
    }
    coda_object_chains[c(parName)] <- as.data.frame(onecolumn_chains)
  }
  coda_object_simplified <- melt(coda_object_chains)
  if (sampled == TRUE & nrow(coda_object_simplified) > 1000000 ) {
    coda_object_simplified <- dplyr::sample_n(coda_object_simplified, 1000000) # to avoid too large datasets
  }
  return(coda_object_simplified)
}

# Gompertz survival function, after Pflaumer 2011
#' @rdname amtl_bayes_helper
#' @export
gomp_lx <- function(x, a, b) {
  lx <- exp(a/b - a/b * exp(b * x))
  return(lx)
}

# Gompertz function for estimating life expectancy at age x, after Pflaumer 2011
#' @rdname amtl_bayes_helper
#' @export
gomp_ex <- function(x, a, b) {
  ex <- (-1) * (0.577221566 + log(a/b) + b * x - a/b * exp(b * x) +
                  (a/b * exp(b * x))^2 / 4) / b /exp(a/b - a/b * exp(b * x))
  return(ex)
}

# mean years lived after age x, following Sasaki/Kondo 2016, by integrating area below curve
#' @rdname amtl_bayes_helper
#' @export
gomp_mean_age <- function(x, a, b, lower = 0, upper = Inf) {
  f3 <- function(x) {
    gomp_lx(x, a, b)
  }
  mean_age <- round(integrate(f3, lower = lower, upper = upper)$value, 1)
  return(mean_age)
}

# simulate a population with tooth loss, choose a link function
amtl.sim.link <- function( y = 1000,
                           j = -4.4289,
                           q = 0.065,
                           M = 40,
                           b_min = 0.025,
                           b_max = 0.1,
                           link = "probit") { # options "weibull", "probit" or "logit"
  M_1 <- 0
  M_2 <- 0
  while ( (M < M_1 | M > M_2 )) {
    b_ <- runif(n = 1, min = b_min, max = b_max)
    a_ <- exp(rnorm(1, (-66.77 * (b_ - 0.0718) - 7.119), 0.0823))
    M_ <- 1 / b_ * log (b_/a_) + 15
    M_1 <- M_ - 1
    M_2 <- M_ + 1
  }
  tooth_loss <- data.frame(age = NA, amtl = NA, tp = NA)
  for (i in 1:y) {
    x <- round(flexsurv::rgompertz(1, b_, a_) ) + 15
    if (link == "logit") {
      amtl <- 32 - round(32 / (exp(q * rnorm(1, x, 25) + j) + 1) )
    } else if (link == "logit")  {
      amtl <- round(32 * (pnorm(q * rnorm(1, x, 25) + j)) )
    } else if (link == "weibull") {
      #amtl <- round(32 * pweibull(exp(-1.99782674 + 0.02717507 * x), exp(0.67653669)) )
      amtl <- round(32 * pweibull( rnorm(1, x, 25), j , q) )
    } else if (link == "exponential") {
      t <- rnorm(1, x, 25)
      if(x>75) {
        t <- 75
      }
      amtl <- round(32 * exp(j + q * t))
    }
    if(amtl > 32) {
      amtl <- 32
    } else if(amtl < 0) {
      amtl <- 0
    }
    tp <- 32
    tooth_loss <- rbind(tooth_loss, c(x, amtl, tp))
  }
  tooth_loss <- tooth_loss[-1,]
  return(tooth_loss)

}

# this function generates starting values for the Gompertz distribution
# if the starting age is not 15
gomp.a0 <- function(
    sampling = 100000,
    b_min = 0.02,
    b_max = 0.1,
    minimum_age = 15) {

  # we do not want too much overhead so no computation if the default age of 15 is true
  if (minimum_age == 15) {
    fit_coeff <- c(-66.77, -2.324914, 0.0823)
  } else {
    null_age <- minimum_age - 15

    ind_df <- data.frame(b = runif(n = sampling, min = b_min, max = b_max)) %>%
      mutate(a = exp(rnorm(n(), (-66.77 * (b - 0.0718) - 7.119), sqrt(0.0823) ))) %>%
      mutate(a0 = a * exp(b * null_age))

    fit <- lm(log(a0) ~ b, data = ind_df)
    rse <- sum(fit$residuals**2)/fit$df.residual # without squaring
    fit_coeff <- c(fit$coefficients[2], fit$coefficients[1], rse )
    fit_coeff <- unname(fit_coeff)
  }
  return(fit_coeff)
}

# generate output for Deviation Information Criterion (DIC)
get.dic <- function(x) {
  dic_result <- extract.runjags(x, what = "dic")
  ped_result <- extract.runjags(x, what = "ped")
  dic_df <- data.frame("D_mean" = round(sum(dic_result$deviance), 2),
                       "pD" = round(sum(dic_result$penalty), 2),
                       "DIC" = round(sum(dic_result$deviance) + sum(dic_result$penalty), 2),
                       "popt" = round(sum(ped_result$penalty), 2),
                       "PED" = round(sum(ped_result$deviance) + sum(ped_result$penalty), 2)
  )
  return(dic_df)
}

amtl.diagnostics.max <- function(x) {
  result <- data.frame(PSRF_max = max(x$`PSRF Point est.`[which(x$MCSE > 0)]),
                       PSRF_upper_max = max(x$`PSRF Upper C.I.`[which(x$MCSE > 0)]),
                       ESS_min = min(x$ESS[which(x$MCSE > 0)]
                       ))
  return(result)
}

# function for plotting 2 x 2 plots
bay.ta.plot <- function(x,
                        age_identifier = "age.s",
                        known_age, mean_choice
) {
  x$chosen_mean <- x[,mean_choice]
  age_identifier_grep <- ifelse(age_identifier == "age.s", "^age.s\\[", "^age.s_c")
  x_red <- x[grep(age_identifier_grep,rownames(x)),]
  x_red$known_age <- known_age
  known_age_density <- density(x_red$known_age, bw = 5)
  known_age_density_df <- data.frame(x = known_age_density$x, y = known_age_density$y)

  x_length <- length(known_age)
  x_ordered <- x_red[order(x_red$known_age),]
  x_ordered$id <- c(1:x_length)

  alpha <- x['a','chosen_mean']
  beta <- x['b','chosen_mean']

  plot1 <- ggplot(x_ordered, aes(x = id)) +
    geom_errorbar(aes(ymin=HDIlow, ymax=HDIhigh, color= known_age > HDIlow &
                        known_age < HDIhigh), lwd = 0.5 ,width = 0) +
    geom_point(aes(y=known_age), shape = 3, colour = "black") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_colour_manual(name = 'Age in range',
                        values = setNames(c('chartreuse4','coral2'),c(T, F))) +
    xlab("\nIndividuals ordered by known age-at-death") +
    ylab("HDIlow to HDIhigh, known age-at-death\n") +
    theme_light()

  plot2 <- ggplot() +
    geom_line(data = known_age_density_df,
              aes(x = x, y = y, col = "density of actual ages\n(bw = 5)\n")) +
    xlim(15, 100) +
    geom_function(fun = function(x) flexsurv::dgompertz(x - 15, beta, alpha),
                  aes(col = "Gompertz parameters\nfrom estimates")) +
    theme_light() +
    xlab("age-at-death") +
    ylab("Density\n") +
    scale_colour_manual(values = c("red","black")) +
    theme( legend.title = element_blank(),
           legend.spacing.y = unit(1.0, 'cm')) +
    guides(fill = guide_legend(byrow = F))

  plot3 <- ggplot (x_ordered, aes(x = known_age, y = chosen_mean)) +
    geom_point(shape = 21) + xlim(15,95) + ylim(15,95) +
    geom_smooth(method='lm', formula= y~x) + theme_light() +
    xlab("\nknown age-at-death") +
    ylab("Mode of estimated age-at-death\n") +
    geom_abline(slope = 1, intercept = 0, linetype = 3)

  plot4 <- ggplot(x_ordered, aes(x=known_age, y = ( known_age - chosen_mean))) +
    geom_point(shape = 21) + geom_smooth(method='lm', formula= y~x) +
    geom_hline(yintercept = 0, linetype = 3) + theme_light() +
    xlab("\nknown age-at-death") + ylab("Residual\n")

  plot_result <- ggpubr::ggarrange( plot1, plot2, plot3, plot4,nrow = 2, ncol = 2)

  return(plot_result)
}
