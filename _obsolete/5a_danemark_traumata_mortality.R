if ( runNewMCMC ) {
  set.seed(401)
  dan_sites_mort_trauma_MCMC <- gomp.anthr_age.cat.hierach(
    x = dan,
    category = "trauma_presence",
    age_beg = "Low",
    age_end = "High",
    parameters = c("b_all", "a", "b","M"),
    thinSteps = 1, minimum_age = 12,
    runjagsMethod = "parallel",
    numSavedSteps = 200000)
  save(dan_sites_mort_trauma_MCMC, file = file.path(".", saveFileDir,
                                                   "dan_sites_mort_trauma_MCMC.Rdata"))
} else {
  load(file = file.path(".", saveFileDir, "dan_sites_mort_trauma_MCMC.Rdata"))
}

dan_sites_mort_trauma_MCMC_list = as.mcmc.list( dan_sites_mort_trauma_MCMC )

dan_sites_mort_trauma_MCMC_diag <- diagnostic.summary(dan_sites_mort_trauma_MCMC_list,
                                                     HDImass = 0.95)

mcmc_array <- as_draws_array(dan_sites_mort_trauma_MCMC_list)
draws_df <- as_draws_df(mcmc_array)
draws_df %>%
  dplyr::select(`M[1]`, `M[2]`) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = recode(parameter,
                            `M[1]` = "without trauma",
                            `M[2]` = "with trauma")) %>%
  ggplot(aes(x = value, y = parameter)) +
  stat_halfeye(.width = 0.95, point_interval = mode_hdi, fill = "skyblue") +
  labs(x = "Modal age-at-death (years)",y = NULL) + theme_light() +
  xlim(0,70)
