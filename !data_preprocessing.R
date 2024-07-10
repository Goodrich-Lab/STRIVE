# Set up data for analysis
# library(tidylog)
source(fs::path(here::here("!libraries.R")))
source(fs::path(here::here("!directories.R")))

# Read Data ---------------------------------------------------------------
data <- readxl::read_xlsx(fs::path(
  dir_data,
  "STRIVE_PFAS_demo_08july24.xlsx")) %>% 
  janitor::clean_names()%>%
  mutate(sex = ifelse(sex == 1, "Male", "Female"),
         rural = ifelse(rural == 1, "Living in rural area", "Living in Metro area"),
         smoking = ifelse(smoking == 1, "Smoke or use vape", "Don't smoke or use vape"),
         case_control = ifelse(case_control == 1, "With cirrhosis", "Healthy"),
         ## added by BS
         ethnicity = ifelse(ethnicity == 1, "Hispanic", "Not Hispanic"),
         sq_drink_alcohol = case_when(
           sq_drink_alcohol == 1 ~ "Yes, current drinker",
           sq_drink_alcohol == 2 ~ "No, former drinker (stopped)",
           sq_drink_alcohol == 3 ~ "No, never drinker",
           TRUE ~ "Unknown/Not Reported"), 
         sq_average_drink_per_day = case_when(
             sq_average_drink_per_day == 1 ~ "Less than 1 alcoholic drink per day",
             sq_average_drink_per_day == 2 ~ "1-2 alcoholic drinks per day",
             sq_average_drink_per_day == 3 ~ "3-4 alcoholic drinks per day",
             sq_average_drink_per_day == 4 ~ "More than 4 alcoholic drinks per day",
             TRUE ~ "Unknown/Not Reported"
           )
         ## to here
         )%>%
<<<<<<< HEAD
  mutate_at(.vars = c("sq_self_hep_b","sq_self_hep_c","supp_meds_tylenol",
                      "supp_meds_steroids","sq_water_well","sq_water_tap_unfiltered",
                      "sq_water_house_filtration","sq_water_faucet_filter",
                      "sq_water_charcoal_filter","sq_water_bottled",
                      "sq_water_none","sq_water_other_type","sq_water_dont_know"),
            .funs = ~ifelse(. == 1, "Yes", "No"))%>%
  mutate_at(.vars = c("rural", "smoking", "ethnicity", "sq_self_hep_b","sq_self_hep_c","supp_meds_tylenol",
                      "supp_meds_steroids","sq_water_well","sq_water_tap_unfiltered",
                      "sq_water_house_filtration","sq_water_faucet_filter",
                      "sq_water_charcoal_filter","sq_water_bottled",
                      "sq_water_none","sq_water_other_type","sq_water_dont_know"),
            .funs = ~ifelse(is.na(.), "Unknown/Not Reported", .))%>%
=======
  mutate(rural = ifelse(is.na(rural), "Unknown/Not Reported", rural),
         smoking = ifelse(is.na(smoking), "Unknown/Not Reported", smoking),
         ethnicity = ifelse(is.na(ethnicity), "Unknown/Not Reported", ethnicity))%>%
>>>>>>> efe0c442ef0eb11430bc2f1f4c41b0b9c5a1a9ba
  ## added or edited by BS
  mutate(race_eth_label = ifelse(is.na(race_eth_label), "Unknown/Not Reported", race_eth_label),
         race_final_label = ifelse(is.na(race_final_label), "Unknown/Not Reported", race_final_label),
         sq_average_drink_per_day = ifelse(is.na(sq_average_drink_per_day), "Unknown/Not Reported", sq_average_drink_per_day))
## to here
pfas_name <- colnames(data)[3:27]

# For the analysis, only focusing on the PFAS with <25% below LOD
pfas_name_analysis <- c("pf_hx_s","pfda","pfna","pfos",
                        "pf_hp_a","pfbs", "pfoa","pf_pe_a",
                        "pf_un_a","pf_hp_s","pf_do_a","pf_pe_s",
                        # "x9cl_pf3ons", 
                        "pf_hx_a", "pfba")

legacy <- c("pf_hx_s","pfda","pfna","pfos",
            "pf_hp_a","pfoa","pf_un_a","pf_hp_s",
            "pf_do_a")

legacy_cat <- paste0(legacy, "_median")

emerging <- c("pfbs", "pf_pe_a","pf_pe_s", "pfba", "pf_hx_a")

emerging_cat <- c("pfbs_median", "pf_pe_a_median","pf_pe_s_median", "pfba_detected", "pf_hx_a_detected")


# Imputation of PFAS: min(pfas concentration)/sqrt(2)
data_imputed <- data %>%
  mutate_at(.vars = pfas_name_analysis,
            .funs = ~ifelse(is.na(.),unique(sort(.))[1]/sqrt(2),.))%>%
  mutate_at(.vars = c("pfba", "pf_hx_a"),
            .funs = list(detected = ~ifelse(. == sort(.)[1], 0, 1)))%>%
  mutate_at(.vars = c(legacy, emerging[-c(4,5)]),
            .funs = list(median = ~ifelse(. <median(.), 0, 1)))


# Adding normalized PFAS
data_scaled <- data_imputed %>% 
  mutate_at(.vars = pfas_name_analysis,
            .funs = list(scld = ~scale(.)%>% as.vector(.)))

pfas_name_scld <- paste0(pfas_name_analysis,"_", "scld")
legacy_scld <- paste0(legacy,"_", "scld")
emerging_scld <- paste0(emerging, "_", "scld")

# Adding new cirrhosis outcome defined by AST/ALT (Cirrhosis: AST/ALT ratio > 1; Health: AST/ALT ratio <= 1)
data_scaled1 <- data_scaled %>%
  mutate(`AST/ALT` = ast_u_l/alt_u_l,
         cirrhosis = ifelse(`AST/ALT` > 1, "With cirrhosis", "Healthy"))

write_csv(data_scaled1, fs::path(dir_data, "cleaned_data/STRIVE_cleaned_data.csv"))

rm(data)
rm(data_imputed)
rm(data_scaled)
rm(data_scaled1)
