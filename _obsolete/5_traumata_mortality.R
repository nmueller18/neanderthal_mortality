if ( runNewMCMC ) {
  set.seed(401)
  nea_sap_taxon_trauma_MCMC <- gomp.anthr_age.cat.hierach(x = nea_sap,
                                           age_beg = "age_min",
                                           age_end = "age_max",
                                           category = "taxon_trauma",
                                           parameters = c("b_all", "a", "b","M"),
                                           thinSteps = 1, minimum_age = 12,
                                           runjagsMethod = "parallel",
                                           numSavedSteps = 200000)
  save(nea_sap_taxon_trauma_MCMC, file = file.path(".", saveFileDir,
                                            "nea_sap_taxon_trauma_MCMC.Rdata"))
} else {
  load(file = file.path(".", saveFileDir, "nea_sap_taxon_trauma_MCMC.Rdata"))
}
nea_sap_taxon_trauma_MCMC_list = as.mcmc.list( nea_sap_taxon_trauma_MCMC )
nea_sap_taxon_trauma_MCMC_diag <- diagnostic.summary(nea_sap_taxon_trauma_MCMC_list,
                                              HDImass = 0.95)

mcmc_array <- as_draws_array(nea_sap_taxon_trauma_MCMC_list)
draws_df <- as_draws_df(mcmc_array)
nea_sap_mort_trauma_plot <- draws_df %>%
  dplyr::select(`M[1]`, `M[2]`, `M[3]`, `M[4]`) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  mutate(parameter = recode(parameter,
                            `M[1]` = "Homo sapiens",
                            `M[2]` = "H. s. with trauma",
                            `M[3]` = "Neanderthals",
                            `M[4]` = "N. with trauma")) %>%
  ggplot(aes(x = value, y = parameter)) +
  stat_halfeye(.width = 0.95, point_interval = mode_hdi, fill = "skyblue") +
  labs(x = "Modal age-at-death (years)",y = NULL) + theme_light()
