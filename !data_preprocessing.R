# Set up data for analysis
# library(tidylog)
source(fs::path(here::here("!libraries.R")))
source(fs::path(here::here("!directories.R")))
source(fs::path(here::here("fxn_baysian_below_lod_imputation.R")))

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
  ## added or edited by BS
  mutate(race_eth_label = ifelse(is.na(race_eth_label), "Unknown/Not Reported", race_eth_label),
         race_final_label = ifelse(is.na(race_final_label), "Unknown/Not Reported", race_final_label),
         sq_average_drink_per_day = ifelse(is.na(sq_average_drink_per_day), "Unknown/Not Reported", sq_average_drink_per_day))
## to here
pfas_name <- colnames(data)[3:27]

# For the analysis, only focusing on the PFAS with <25% below LOD----------
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

# Imputation of PFAS: min(pfas concentration)/sqrt(2)-----------
#data_imputed <- data %>%
 # mutate_at(.vars = pfas_name_analysis,
  #          .funs = ~ifelse(is.na(.),unique(sort(.))[1]/sqrt(2),.))%>%
#  mutate_at(.vars = c("pfba", "pf_hx_a"),
 #           .funs = list(detected = ~ifelse(. == sort(.)[1], 0, 1)))%>%
#  mutate_at(.vars = c(legacy, emerging[-c(4,5)]),
 #           .funs = list(median = ~ifelse(. <median(.), 0, 1)))

# Calculate imputed vars ---------------------- 
# Set seed for reproducibility
set.seed(1234)

data_imputed <- data %>%
  mutate(across(.cols= pfas_name_analysis, .fns = log))

data_imputed <- data_imputed %>%
  mutate(pf_hx_s = impute.X.BDL(X.obs=pf_hx_s, LOD=min(pf_hx_s, na.rm = TRUE)),
         pfda = impute.X.BDL(X.obs=pfda, LOD=min(pfda, na.rm = TRUE)),
         pfna = impute.X.BDL(X.obs=pfna, LOD=min(pfna, na.rm = TRUE)),
         pfos = impute.X.BDL(X.obs=pfos, LOD=min(pfos, na.rm = TRUE)),
         pf_hp_a = impute.X.BDL(X.obs=pf_hp_a, LOD=min(pf_hp_a, na.rm = TRUE)),
         pfbs = impute.X.BDL(X.obs=pfbs, LOD=min(pfbs, na.rm = TRUE)),
         pfoa = impute.X.BDL(X.obs=pfoa, LOD=min(pfoa, na.rm = TRUE)),
         pf_pe_a = impute.X.BDL(X.obs=pf_pe_a, LOD=min(pf_pe_a, na.rm = TRUE)),
         pf_un_a = impute.X.BDL(X.obs=pf_un_a, LOD=min(pf_un_a, na.rm = TRUE)),
         pf_hp_s = impute.X.BDL(X.obs=pf_hp_s, LOD=min(pf_hp_s, na.rm = TRUE)),
         pf_do_a = impute.X.BDL(X.obs=pf_do_a, LOD=min(pf_do_a, na.rm = TRUE)),
         pf_pe_s = impute.X.BDL(X.obs=pf_pe_s, LOD=min(pf_pe_s, na.rm = TRUE)),
         pfba = impute.X.BDL(X.obs=pfba, LOD=min(pfba, na.rm = TRUE)),
         pf_hx_a = impute.X.BDL(X.obs=pf_hx_a, LOD=min(pf_hx_a, na.rm = TRUE))
  )

data_imputed <- data_imputed %>%
  mutate(across(.cols = pfas_name_analysis, .fns = exp))


# Adding normalized and log2 transformed PFAS------------
data_scaled <- data_imputed %>% 
  mutate_at(.vars = pfas_name_analysis,
            .funs = list(scld = ~scale(.)%>% as.vector(.),
                         log2 = ~log2(.) %>% as.vector()))

pfas_name_scld <- paste0(pfas_name_analysis,"_", "scld")
legacy_scld <- paste0(legacy,"_", "scld")
emerging_scld <- paste0(emerging, "_", "scld")


rm(data)
rm(data_imputed)

# Replace Old Covars data with New received ones------- 
library(haven)
# load new data
covars_df <- read_csv(fs::path(
  dir_data,
  "2025 Newly updated data/strive_rqst_covs_bojung.csv"))%>%
  janitor::clean_names() 

# Variables needed to be replaced
vars_name <- c("source","age_at_enrollment","sex","race_final_label","smoking","rural",
               "sq_drink_alcohol","sq_average_drink_per_day","sq_self_hep_b",
               "sq_self_hep_c","sq_water_well","sq_water_tap_unfiltered",
               "sq_water_house_filtration","sq_water_faucet_filter",
               "sq_water_charcoal_filter","sq_water_bottled","sq_water_none",
               "sq_water_other_type","sq_water_dont_know")

data1 <- data_scaled %>% tidylog::left_join(covars_df)%>%
  dplyr::select(-vars_name) %>%
  rename(source = site,
         age_at_enrollment = slm1_age,
         sex = slm1_sex,
         race_final_label = slm1_race_ethnic,
         rural = slm1_rucc_metro ,
         smoking = s1_smoked_cigs,
         sq_drink_alcohol = s1_drink_alcohol,
         sq_average_drink_per_day = s1_average_drink_per_day,
         sq_self_hep_b = s1_self_hep_b,
         sq_self_hep_c = s1_self_hep_c,
         sq_water_well = s1_water_well,
         sq_water_tap_unfiltered= s1_water_tap_unfiltered,
         sq_water_house_filtration = s1_water_house_filtration,
         sq_water_faucet_filter = s1_water_faucet_filter,
         sq_water_charcoal_filter =s1_water_charcoal_filter,
         sq_water_bottled = s1_water_bottled,
         sq_water_none = s1_water_none,
         sq_water_other_type = s1_water_other,
         sq_water_dont_know= s1_water_dont_know)%>%
  mutate(sex = ifelse(sex == 1, "Male", "Female"),
         # rural = ifelse(rural == 1, "Living in rural area", "Living in Metro area"),
         rural = ifelse(rural <= 3, "Metro Counties", "Non-Metro Counties"),
         smoking = case_when(smoking == 1 ~ "Everyday", 
                             smoking == 2|smoking == 3 ~ "Some days",
                             smoking == 4 ~ "Not at all",
                             TRUE ~ "Unknow/Not Reported"),
         s1_use_vape = case_when(s1_use_vape == 1 ~ "Daily", 
                                 s1_use_vape == 2 ~ "At least weekly (but not daily)",
                                 s1_use_vape == 3 ~ "Less often than weekly",
                                 s1_use_vape == 4 ~ "Not at all",
                                 TRUE ~ "Unknow/Not Reported"),
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
         ),
         race_final_label = case_when(
           race_final_label == 1 ~ "American Indian or Alaska Native",
           race_final_label == 2 ~ "Asian",
           race_final_label == 3 ~ "Black or African American",
           race_final_label == 4 ~ "Hispanic or Latino",
           race_final_label == 6 ~ "White",
           race_final_label == 8 ~ "Unknown/Not Reported"
         )
         ## to here
  )%>%
  mutate_at(.vars = c("sq_self_hep_b","sq_self_hep_c","supp_meds_tylenol",
                      "supp_meds_steroids","sq_water_well","sq_water_tap_unfiltered",
                      "sq_water_house_filtration","sq_water_faucet_filter",
                      "sq_water_charcoal_filter","sq_water_bottled",
                      "sq_water_none","sq_water_dont_know",
                      "s1_meds_pain_reliever_current","s1_meds_pain_reliever_take"),
            .funs = ~ifelse(. >= 1, "Yes", "No"))%>%
  mutate_at(.vars = c("rural", "smoking", "ethnicity", "sq_self_hep_b","sq_self_hep_c","supp_meds_tylenol",
                      "supp_meds_steroids","sq_water_well","sq_water_tap_unfiltered",
                      "sq_water_house_filtration","sq_water_faucet_filter",
                      "sq_water_charcoal_filter","sq_water_bottled",
                      "sq_water_none","sq_water_dont_know",
                      "s1_meds_pain_reliever_current","s1_meds_pain_reliever_take"),
            .funs = ~ifelse(is.na(.), "Unknown/Not Reported", .))%>%
  mutate(sq_water_other_type = case_when(sq_water_other_type %in% 
                                        c("N/A", "na", "none", "None", "NONE") ~ "No",
                                        !is.na(sq_water_other_type)~"Yes",
                                        is.na(sq_water_other_type) ~ "Unknown/Not Reported")) %>%
  ## added or edited by BS
  mutate(race_eth_label = ifelse(is.na(race_eth_label), "Unknown/Not Reported", race_eth_label),
         race_final_label = ifelse(is.na(race_final_label), "Unknown/Not Reported", race_final_label),
         sq_average_drink_per_day = ifelse(is.na(sq_average_drink_per_day), "Unknown/Not Reported", sq_average_drink_per_day))


# Replace Old Alt data with new ones------

# load new alt data
metab_df <- read_csv(fs::path(
  dir_data,
  "2025 Newly updated data/strive_rqst_metab_Bojung.csv"))%>%
  janitor::clean_names()

# Adding new cirrhosis outcome defined by AST/ALT (Cirrhosis: AST/ALT ratio > 1; Health: AST/ALT ratio <= 1) and other alt variables-----------
data2 <- data1 %>% 
  dplyr::select(-intersect(colnames(metab_df), colnames(data1))[-1], -log_alt, -log_ast, -log_trig, -trig_mg_d_l)%>%
  tidylog::left_join(metab_df, by = "strive_id")%>%
  mutate(`AST/ALT` = ast_u_l/alt_u_l,
         cirrhosis = ifelse(`AST/ALT` > 1, "With cirrhosis", "Healthy"),
         log_alt = log(alt_u_l),
         log_ast = log(ast_u_l),
         log_trig = log(trig_mg_dl))%>%
  mutate(alt_cat1 = ifelse((alt_u_l <= 29 & sex == "Male")|
                             (alt_u_l <= 19 & sex == "Female"), "Norm", "High"),
         alt_cat2 = ifelse((alt_u_l <= 33 & sex == "Male")|
                             (alt_u_l <= 25 & sex == "Female"), "Norm", "High"),
         ast_cat1 = ifelse((ast_u_l <= 29 & sex == "Male")|
                             (ast_u_l <= 19 & sex == "Female"), "Norm", "High"),
         ast_cat2 = ifelse((ast_u_l <= 33 & sex == "Male")|
                             (ast_u_l <= 25 & sex == "Female"), "Norm", "High"))%>%
  mutate(
    tg_ast = ifelse(trig_flags == "High" & ast_cat2 == "High", "High", "Norm"),
    tg_ast_alt = ifelse(trig_flags == "High" & alt_cat2 == "High" & ast_cat2 == "High", "High", "Norm"),
    alt_ast = ifelse(alt_cat2 == "High" & ast_cat2 == "High", "High", "Norm"))%>%
   mutate(obesity_status = ifelse(slm1_bmi >= 30, "Obesity", "Non-obesity"))

write_csv(data2, fs::path(dir_data, "cleaned_data/STRIVE_cleaned_data_V2.csv"))

