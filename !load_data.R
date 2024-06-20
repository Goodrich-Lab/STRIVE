# Set up data for analysis
# library(tidylog)

# Read Data ---------------------------------------------------------------

data <- readxl::read_xlsx(fs::path(
  dir_data,
  "STRIVE_PFAS_demo_May24.xlsx")) %>% 
  janitor::clean_names()

pfas_name <- colnames(data)[3:27]

# For the analysis, only focusing on the PFAS with <25% below LOD
pfas_name_analysis <- c("pf_hx_s","pfda","pfna","pfos",
                        "pf_hp_a","pfbs", "pfoa","pf_pe_a",
                        "pf_un_a","pf_hp_s","pf_do_a","pf_pe_s")

legacy <- c("pf_hx_s","pfda","pfna","pfos",
            "pf_hp_a","pfoa","pf_un_a","pf_hp_s",
            "pf_do_a")

emerging <- c("pfbs", "pf_pe_a","pf_pe_s")

covars <- c("source", "age_at_enrollment","sex", 
            "race_eth_label", "rural", "smoking")

covars_analysis <- covars[1:4]

# Imputation of PFAS: min(pfas concentration)/sqrt(2)
data_imputed <- data %>%
  mutate_at(.vars = pfas_name_analysis,
            .funs = ~ifelse(is.na(.),unique(sort(.))[2]/sqrt(2),.))

# Adding normalized PFAS
data_scaled <- data_imputed %>% 
  mutate_at(.vars = pfas_name_analysis,
            .funs = list(scld = ~scale(.)%>% as.vector(.)))%>%
  mutate(smoking = ifelse(is.na(smoking), "Unknown", smoking),
         rural = ifelse(is.na(rural), "Unknown", rural),
         race_eth_label = ifelse(is.na(race_eth_label), "Unknown", race_eth_label))


pfas_name_scld <- paste0(pfas_name_analysis,"_", "scld")
legacy_scld <- paste0(legacy,"_", "scld")
emerging_scld <- paste0(emerging, "_", "scld")
# 
