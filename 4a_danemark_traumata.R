if ( runNewMCMC ) {
  set.seed(5754)
  dan_trauma_MCMC <- traumata(
                          x = dan,
                          traumata = "trauma_presence",
                          age_beg = "Low",
                          age_end = "High",
                          parameters = c( "a", "b","M","x0", "scale"),
                          thinSteps = 1, minimum_age = 15,
                          runjagsMethod = "parallel",
                          numSavedSteps = 200000)
  save(dan_trauma_MCMC, file = file.path(".", saveFileDir,
                                     "dan_trauma_MCMC.Rdata"))
} else {
  load(file = file.path(".", saveFileDir, "dan_trauma_MCMC.Rdata"))
}
dan_trauma_MCMC_list = as.mcmc.list( dan_trauma_MCMC )
dan_trauma_MCMC_diag <- diagnostic.summary(dan_trauma_MCMC_list,
                                       HDImass = 0.95)

if ( runNewMCMC ) {
  set.seed(34)
  dan_age_trauma_MCMC <- gomp.anthr_age.cat(
                              x = dan,
                              category = "trauma_presence",
                              age_beg = "Low",
                              age_end = "High",
                              parameters = c( "a", "b","M","x0", "scale"),
                              thinSteps = 1, minimum_age = 15,
                              runjagsMethod = "parallel",
                              numSavedSteps = 200000)
  save(dan_age_trauma_MCMC, file = file.path(".", saveFileDir,
                                         "dan_age_trauma_MCMC.Rdata"))
} else {
  load(file = file.path(".", saveFileDir, "dan_age_trauma_MCMC.Rdata"))
}
dan_age_trauma_MCMC_list = as.mcmc.list( dan_age_trauma_MCMC )
dan_age_trauma_MCMC_diag <- diagnostic.summary(dan_age_trauma_MCMC_list,
                                           HDImass = 0.95)
mcmc_array <- as_draws_array(dan_age_trauma_MCMC_list)
draws_df <- as_draws_df(mcmc_array)
dan_age_trauma_plot <- draws_df %>%
  dplyr::select(`M[1]`, `M[2]`) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = recode(parameter,
                            `M[1]` = "without trauma",
                            `M[2]` = "with trauma")) %>%
  ggplot(aes(x = value, y = parameter)) +
  stat_halfeye(.width = 0.95, point_interval = mode_hdi, fill = "skyblue") +
  labs(x = "Modal age (years)",y = NULL) + theme_light() + xlim(10,60)


if ( runNewMCMC ) {
  set.seed(233)
  dan_sites_trauma_MCMC <- traumata.cat(
                                      x = dan,
                                      traumata = "trauma_presence",
                                      age_beg = "Low",
                                      age_end = "High",
                                      category = "Site",
                                      parameters = c( "a", "b","M","x0", "scale"),
                                      thinSteps = 1, minimum_age = 15,
                                      runjagsMethod = "parallel",
                                      numSavedSteps = 150000)
  save(dan_sites_trauma_MCMC, file = file.path(".", saveFileDir,
                                             "dan_sites_trauma_MCMC.Rdata"))
} else {
  load(file = file.path(".", saveFileDir, "dan_sites_trauma_MCMC.Rdata"))
}
dan_sites_trauma_MCMC_list = as.mcmc.list( dan_sites_trauma_MCMC )
dan_sites_trauma_MCMC_diag <- diagnostic.summary(dan_sites_trauma_MCMC_list,
                                               HDImass = 0.95)

mcmc_array <- as_draws_array(dan_sites_trauma_MCMC_list)
draws_df <- as_draws_df(mcmc_array)
dan_sites_trauma_plot <- draws_df %>%
  dplyr::select(`x0[1]`, `x0[2]`, `x0[3]`) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = recode(parameter,
                            `x0[1]` = "Sortebrødre",
                            `x0[2]` = "St. Mikkel",
                            `x0[3]` = "Tirup")) %>%
  ggplot(aes(x = value, y = parameter)) +
  stat_halfeye(.width = 0.95, point_interval = mode_hdi, fill = "skyblue") +
  labs(x = "Logistic location of traumata (years)",y = NULL) + theme_light() +
  xlim(50,120)
