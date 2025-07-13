if ( runNewMCMC ) {
  set.seed(711)
  nea_sap_anthr_MCMC <- gomp.anthr_age.cat(x = nea_sap,
                                   age_beg = "age_min",
                                   age_end = "age_max",
                                   category = "taxon",
                                   thinSteps = 1, minimum_age = 12,
                                   runjagsMethod = "parallel",
                                   numSavedSteps = 100000)
  save(nea_sap_anthr_MCMC, file = file.path(".", saveFileDir,
                                        "nea_sap_anthr_MCMC.Rdata"))
} else {
  load(file = file.path(".", saveFileDir, "nea_sap_anthr_MCMC.Rdata"))
}
nea_sap_anthr_MCMC_list = as.mcmc.list( nea_sap_anthr_MCMC )
nea_sap_anthr_MCMC_diag <- diagnostic.summary(nea_sap_anthr_MCMC_list,
                                          HDImass = 0.95)
mcmc_array <- as_draws_array(nea_sap_anthr_MCMC_list)
draws_df <- as_draws_df(mcmc_array)
nea_sap_mort_plot <- draws_df %>%
  dplyr::select(`M[1]`, `M[2]`) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = recode(parameter,
                            `M[1]` = "Homo sapiens",
                            `M[2]` = "Neanderthals")) %>%
  ggplot(aes(x = value, y = parameter)) +
  stat_halfeye(.width = 0.95, point_interval = mode_hdi, fill = "skyblue") +
  labs(x = "Modal age-at-death (years)",y = NULL) + theme_light() + xlim(-20,40)
