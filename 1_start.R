require(pacman) || install.packages("pacman")
pacman::p_load(coda, dplyr, ggplot2, Metrics, ggpubr, rjags, runjags, flexsurv,
               scoringRules, MASS, parallel, nimble, tidybayes, posterior, tidyr,
               mortAAR)

options(dplyr.summarise.inform = FALSE)
source("./functions/bayes_helper.R")
source("./functions/anthro_age.R")
source("./functions/anthro_age_cat.R")
source("./functions/traumata.R")
source("./functions/traumata_cat.R")
source("./2_process_data.R")

# N.B.! Next, specify whether to run new MCMC analyses or to re-load previously
# run and saved MCMC results. When runNewMCMC=TRUE then the MCMC chains are run
# anew. When runNewMCMC=!TRUE (NOT TRUE) then the MCMC chains from the
# previous run are used. Because all MCMC processes are seeded, they should give
# the same results every time they are run.
# By using previously saved chains, all processes (i.e., knit) are much faster.

runNewMCMC = !TRUE

# Specify filename prefix for saved files and create a folder if needed:
saveFileDir = "preprocessed_files"
dir.create(file.path(".", saveFileDir), showWarnings = FALSE )

## Danemark
# Mortality
source("./3a_danemark_mortality.R")
plot(dan_lt[1:3], display = "Dx", line_vis = "colour")
plot(dan_lt[4], display = "Dx", line_vis = "colour")
dan_sites_mortality_plot

#Traumata
source("./4a_danemark_traumata.R")
plot(dan_trauma_lt[1:2], display = "dx", line_vis = "colour")
dan_age_trauma_plot
dan_sites_trauma_plot


## Neanderthals and Homo sapiens
# Mortality
source("./3_mortality.R")
plot(nea_sap_lt[1:2], display = "dx")
nea_sap_mort_plot

# Traumata of Neanderthals and Homo sapiens
source("./4_traumata.R")
nea_sap %>% group_by(taxon, trauma_presence) %>%
  summarize(n = n(), .groups = "drop") %>%
  pivot_wider( names_from = trauma_presence, values_from = n )
nea_sap_traumata_plot
plot(nea_sap_trauma_lt[1:4], display = "dx", line_vis = "colour")
