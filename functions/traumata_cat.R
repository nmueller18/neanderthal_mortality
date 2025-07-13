traumata.cat <- function( x,
                           age_beg,
                           age_end,
                           traumata,
                           category,
                           parameters = c( "a", "b","M"),
                           minimum_age = 15,
                           numSavedSteps = numSavedSteps,
                           thinSteps = thinSteps,
                           runjagsMethod="rjags",
                           nChains=3,
                           adaptSteps = 1000,
                           burnInSteps = 2000,
                           silent.jags = FALSE,
                           silent.runjags = FALSE) {
  age_beg <- x[[age_beg]]
  age_end <- x[[age_end]]
  traumata <- x[[traumata]]
  category <- x[[category]]
  Ntotal <- length(age_beg) # number of individuals
  ones <- rep(1,Ntotal)
  C <- 100000
  c = as.numeric(as.factor(category))
  cN = max(c)

  # Generate values for Gompertz alpha if minimum age is not 15
  gomp_a0 <- gomp.a0(minimum_age = minimum_age)

  # Generate inits-list with custom function,
  # including RNGs and seed for JAGS for reproducible results
  initsList <- function(){
    RNG_list <- c("base::Wichmann-Hill",
                  "base::Marsaglia-Multicarry",
                  "base::Super-Duper",
                  "base::Mersenne-Twister")
    init_list <- list(
      .RNG.name = sample(RNG_list, 1),
      .RNG.seed = sample(1:1e+06, 1),
      b = runif(cN, 0.02, 0.05),
      log_a = runif(cN, 4, 9) # equals 0.0001-0.018
    )
    return(init_list)
  }


  # Specify the data in a list, for later shipment to JAGS:
  # Be careful to specify only parameters which are actually used, otherwise you may confuse JAGS
  dataList = list(
    Ntotal = Ntotal ,
    c = c,
    cN = cN,
    C = C, # JAGS does not warn if too small!
    ones = ones,
    minimum_age = minimum_age,
    age_end = age_end,
    age_beg = age_beg,
    gomp_a0_m = gomp_a0[1],
    gomp_a0_ic = gomp_a0[2],
    gomp_a0_var = gomp_a0[3],
    traumata = traumata
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
    for ( i in 1:Ntotal ) {
      age[i] ~ dunif(age_beg[i] - minimum_age, age_end[i] - minimum_age)
      age.s[i] <- age[i] + minimum_age
      spy[i] <- a[c[i]] * exp(b[c[i]] * age[i]) *
          exp(-a[c[i]]/b[c[i]] * (exp(b[c[i]] * age[i]) - 1)) / C # implementing Gompertz probability density
      ones[i] ~ dbern( spy[i]  )
      traumata[i] ~ dbern( p[i]  )
      p[i] <- 1 / (1 + exp(-logit_p[i])) # location + scale parametrization
      logit_p[i] <- (age.s[i] - x0[c[i]]) / scale[c[i]]
    }
    for (j in 1:cN) {
      x0[j] ~ dnorm( 75 , 1/20^2 )
      scale[j] ~ dnorm( 0 , 1/10^2 ) T(0,)
    }
    for (j in 1:cN) {
      b[j]  ~ dgamma(0.01, 0.01) T(0.02, 0.1)# a must not be null  #dunif(0.02, 0.1)
      log_a_M[j] <- (gomp_a0_m * b[j] + gomp_a0_ic) * (-1)
      log_a[j]  ~ dgamma(log_a_M[j]^2 / gomp_a0_var, log_a_M[j] / gomp_a0_var) #T(log(0.02) * (-1),)
      a[j] <- exp(log_a[j] * (-1))
      M[j] <- 1 / b[j] * log (b[j]/a[j]) + minimum_age
    }

  }
  " # close quote for modelString

  runjags.options(silent.jags = silent.jags,
                  silent.runjags = silent.runjags)

  # RUN THE CHAINS
  runJagsOut <- run.jags( method = runjagsMethod ,
                          model = modelString ,
                          monitor = parameters ,
                          data = dataList ,
                          inits = initsList ,
                          n.chains = nChains ,
                          adapt = adaptSteps ,
                          burnin = burnInSteps ,
                          sample = ceiling(numSavedSteps/nChains) ,
                          thin = thinSteps ,
                          summarise = TRUE ,
                          plots = TRUE )
  codaSamples = as.mcmc.list( runJagsOut )
  return(codaSamples)
}
