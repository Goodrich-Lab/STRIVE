# Set up data for analysis
# library(tidylog)

# Read Data ---------------------------------------------------------------

data <- readxl::read_xlsx(fs::path(
  dir_data,
  "STRIVE_PFAS_demo_28jun24.xlsx")) %>% 
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
  mutate(rural = ifelse(is.na(rural), "Unknown/Not Reported", rural),
         smoking = ifelse(is.na(smoking), "Unknown/Not Reported", smoking))%>%
  drop_na(case_control)

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

# emerging <- c("pfbs", "pf_pe_a","pf_pe_s", "x9cl_pf3ons")

## added or edited by BS
covars <- c("source", "age_at_enrollment","sex", 
            "rural", "smoking","race_final_label", "sq_average_drink_per_day")

covars_analysis <- covars[1:7]
## to here

# Imputation of PFAS: min(pfas concentration)/sqrt(2)
data_imputed <- data %>%
  mutate_at(.vars = pfas_name_analysis,
            .funs = ~ifelse(is.na(.),unique(sort(.))[1]/sqrt(2),.))%>%
  mutate_at(.vars = c("pfba", "pf_hx_a"),
            .funs = list(detected = ~ifelse(. == sort(.)[1], 0, 1)))%>%
  # mutate_at(.vars = c(legacy, emerging[-c(4,5)]),
  #           .funs = list(median = ~ifelse(. < median(.), 0, 1)))
  mutate_at(.vars = c(legacy, emerging[-c(4,5)]),
            .funs = list(median = ~ifelse(. < quantile(.,0.90), 0, 1)))


# Adding normalized PFAS
data_scaled <- data_imputed %>% 
  mutate_at(.vars = pfas_name_analysis,
            .funs = list(scld = ~scale(.)%>% as.vector(.)))%>%
  ## added or edited by BS
  mutate(smoking = ifelse(is.na(smoking), "Unknown", smoking),
         rural = ifelse(is.na(rural), "Unknown", rural),
         race_eth_label = ifelse(is.na(race_eth_label), "Unknown", race_eth_label),
         race_final_label = ifelse(is.na(race_final_label), "Unknown", race_final_label),
         sq_average_drink_per_day = ifelse(is.na(sq_average_drink_per_day), "Unknown", sq_average_drink_per_day))
  ## to here

pfas_name_scld <- paste0(pfas_name_analysis,"_", "scld")
legacy_scld <- paste0(legacy,"_", "scld")
emerging_scld <- paste0(emerging, "_", "scld")
# 


