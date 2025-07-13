if ( runNewMCMC ) {
  set.seed(233)
  nea_sap_trauma_MCMC <- traumata.cat(
                          x = nea_sap,
                          traumata = "trauma_presence",
                          age_beg = "age_min",
                          age_end = "age_max",
                          category = "taxon",
                          parameters = c( "a", "b","M","x0", "scale"),
                          thinSteps = 1, minimum_age = 12,
                          runjagsMethod = "parallel",
                          numSavedSteps = 100000)
  save(nea_sap_trauma_MCMC, file = file.path(".", saveFileDir,
                                     "nea_sap_trauma_MCMC.Rdata"))
} else {
  load(file = file.path(".", saveFileDir, "nea_sap_trauma_MCMC.Rdata"))
}
nea_sap_trauma_MCMC_list = as.mcmc.list( nea_sap_trauma_MCMC )
nea_sap_trauma_MCMC_diag <- diagnostic.summary(nea_sap_trauma_MCMC_list,
                                              HDImass = 0.95)
mcmc_array <- as_draws_array(nea_sap_trauma_MCMC_list)
draws_df <- as_draws_df(mcmc_array)
nea_sap_traumata_plot <- draws_df %>%
  dplyr::select(`x0[1]`, `x0[2]`) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = recode(parameter,
                            `x0[1]` = "Homo sapiens",
                            `x0[2]` = "Neanderthals")) %>%
  ggplot(aes(x = value, y = parameter)) +
  stat_halfeye(.width = 0.95, point_interval = mode_hdi, fill = "skyblue") +
  labs(x = "Logistic location of traumata (years)",y = NULL) + theme_light() +
  xlim(50,120)
