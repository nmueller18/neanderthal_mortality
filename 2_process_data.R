## boldsen et al.: head traumata in Danemark
dan <- read.table("./data/pnas.tsv", header = T, sep = "\t")
dan$trauma_presence <- ifelse(dan$Fracture == "No healed fracture", 0, 1)

dan_prep <- mortAAR::prep.life.table(dan, agebeg = "Low", ageend = "High",
                                     group = "Site", agerange = "excluded")
dan_lt <- mortAAR::life.table(dan_prep)
dan_trauma_prep <- mortAAR::prep.life.table(dan,agebeg = "Low", ageend = "High",
                                            group = "Fracture",
                                            agerange = "excluded")
dan_trauma_lt <- mortAAR::life.table(dan_trauma_prep)


# Neanderthal and Homo sapiens
nea_sap <- readxl::read_excel("./data/Data_NEA_mortality_17.06.2025.xlsx")
nea_sap$age_max <- ifelse(nea_sap$age_max == "∞", 100, nea_sap$age_max)
nea_sap$age_max <- as.numeric(nea_sap$age_max)
nea_sap$taxon_trauma <- paste0(nea_sap$taxon, "_", nea_sap$trauma_presence)

nea_sap_prep <- mortAAR::prep.life.table(nea_sap,agebeg = "age_min", ageend = "age_max",
                        group = "taxon", agerange = "included")
nea_sap_lt <- mortAAR::life.table(nea_sap_prep)

nea_sap_trauma_prep <- mortAAR::prep.life.table(nea_sap,agebeg = "age_min", ageend = "age_max",
                                         group = "taxon_trauma", agerange = "included")
nea_sap_trauma_lt <- mortAAR::life.table(nea_sap_trauma_prep)
