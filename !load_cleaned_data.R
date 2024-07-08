# Set up data for analysis
# library(tidylog)

# Read Data ---------------------------------------------------------------

data <- read_csv(fs::path(
  dir_data,
  "cleaned_data/STRIVE_cleaned_data.csv"))

pfas_name <- colnames(data)[3:27]

# For the analysis, only focusing on the PFAS with <25% below LOD
pfas_name_analysis <- c("pf_hx_s","pfda","pfna","pfos",
                        "pf_hp_a","pfbs", "pfoa","pf_pe_a",
                        "pf_un_a","pf_hp_s","pf_do_a","pf_pe_s",
                        "pf_hx_a", "pfba")

legacy <- c("pf_hx_s","pfda","pfna","pfos",
            "pf_hp_a","pfoa","pf_un_a","pf_hp_s",
            "pf_do_a")

legacy_cat <- paste0(legacy, "_median")

emerging <- c("pfbs", "pf_pe_a","pf_pe_s", "pfba", "pf_hx_a")

emerging_cat <- c("pfbs_median", "pf_pe_a_median","pf_pe_s_median", "pfba_detected", "pf_hx_a_detected")

# emerging <- c("pfbs", "pf_pe_a","pf_pe_s", "x9cl_pf3ons")

covars <- c("source", "age_at_enrollment","sex", 
            "race_eth_label", "rural", "smoking")

covars_analysis <- covars[1:4]


pfas_name_scld <- paste0(pfas_name_analysis,"_", "scld")
legacy_scld <- paste0(legacy,"_", "scld")
emerging_scld <- paste0(emerging, "_", "scld")
# 


