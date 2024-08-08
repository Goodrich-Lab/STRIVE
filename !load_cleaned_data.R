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

potential_conf <- c("source", "age_at_enrollment","sex", "bmi",
             "race_eth_label", "race_final_label","ethnicity",
             "rural", "smoking","trig_mg_d_l","sq_drink_alcohol","sq_average_drink_per_day","sq_self_hep_b",
             "sq_self_hep_c","supp_meds_tylenol","supp_meds_steroids","sq_water_well",
             "sq_water_tap_unfiltered","sq_water_house_filtration", "sq_water_faucet_filter",
             "sq_water_charcoal_filter","sq_water_bottled","sq_water_none","sq_water_other_type")


covars <- c("source", "age_at_enrollment", "sex","smoking", "sq_average_drink_per_day")

pfas_name_scld <- paste0(pfas_name_analysis,"_", "scld")
legacy_scld <- paste0(legacy,"_", "scld")
emerging_scld <- paste0(emerging, "_", "scld")
# 


