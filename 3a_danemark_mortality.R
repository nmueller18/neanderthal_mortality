if ( runNewMCMC ) {
  set.seed(9896)
  dan_anthr_MCMC <- gomp.anthr_age(x = dan,
                                   age_beg = "Low",
                                   age_end = "High",
                                   parameters = c( "a", "b","M"),
                                   minimum_age = 15,
                                   runjagsMethod = "parallel",
                                   thinSteps = 1,
                                   numSavedSteps = 150000)
  save(dan_anthr_MCMC, file = file.path(".", saveFileDir,
                                        "dan_anthr_MCMC.Rdata"))
} else {
  load(file = file.path(".", saveFileDir, "dan_anthr_MCMC.Rdata"))
}

dan_anthr_MCMC_list = as.mcmc.list( dan_anthr_MCMC )
dan_anthr_MCMC_diag <- diagnostic.summary(dan_anthr_MCMC_list,
                                          HDImass = 0.95)

if ( runNewMCMC ) {
  set.seed(711)
  dan_sites_anthr_MCMC <- gomp.anthr_age.cat(x = dan,
                                            age_beg ="Low",
                                           age_end = "High",
                                           category = "Site",
                                           thinSteps = 1, minimum_age = 15,
                                           runjagsMethod = "parallel",
                                           numSavedSteps = 150000)
  save(dan_sites_anthr_MCMC, file = file.path(".", saveFileDir,
                                            "dan_sites_anthr_MCMC.Rdata"))
} else {
  load(file = file.path(".", saveFileDir, "dan_sites_anthr_MCMC.Rdata"))
}

dan_sites_anthr_MCMC_list = as.mcmc.list( dan_sites_anthr_MCMC )
dan_sites_anthr_MCMC_diag <- diagnostic.summary(dan_sites_anthr_MCMC_list,
                                              HDImass = 0.95)

mcmc_array <- as_draws_array(dan_sites_anthr_MCMC_list)
draws_df <- as_draws_df(mcmc_array)
dan_sites_mortality_plot <- draws_df %>%
  dplyr::select(`M[1]`, `M[2]`, `M[3]`) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = recode(parameter,
                            `M[1]` = "Sortebrødre",
                            `M[2]` = "St. Mikkel",
                            `M[3]` = "Tirup")) %>%
  ggplot(aes(x = value, y = parameter)) +
  stat_halfeye(.width = 0.95, point_interval = mode_hdi, fill = "skyblue") +
  labs(x = "Modal age-at-death (years)",y = NULL) + theme_light()
