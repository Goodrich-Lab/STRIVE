fit_bwqs_model <- function(expos,
                           out,
                           covariates, 
                           data,
                           set_prior)
{
  data_bwqs <- data %>% 
    mutate(
      y = !!sym(out)
    ) %>% 
    data.frame
  
  expos_addx <- paste0(expos, "_lg2")
  
  fit_bwqs <- bwqs(
    formula = as.formula(paste0(
      "y",
      " ~ ",
      paste(covariates, collapse = " + "))),
    mix_name = expos_addx,
    data = data_bwqs,
    q = 4, iter = 10000,
    family = "binomial",
    prior = set_prior
  )
  
  output <- fit_bwqs$summary_fit %>% 
    as_tibble(., rownames = "coef_name")
  output <- output %>%
    dplyr::mutate(
      var_name = gsub("C_|W_X|W_", "", coef_name)
    )
  
  return(output)
}

# Clean PFAS names --------
clean_pfas_names <- function(pfas_names) {
  cleaned_names <- pfas_names %>%
    # Capitalize "pf" to "PF"
    gsub("^pf", "PF", .) %>%
    # Remove underscores and capitalize the following letter
    gsub("_([a-z])", "\\U\\1", ., perl = TRUE) %>%
    # Ensure fully uppercase formatting for known PFAS abbreviations
    gsub("PFba", "PFBA", .) %>%
    gsub("PFbs", "PFBS", .) %>%
    gsub("PFda", "PFDA", .) %>%
    gsub("PFna", "PFNA", .) %>%
    gsub("PFoa", "PFOA", .) %>%
    gsub("PFos", "PFOS", .) %>%
    gsub("x9clPf3ons", "9Cl-PF3ONS", .) %>%
    
    # Specific cleaning for known sums
    gsub("sumAllPfas", "Σ all PFAS", .) %>%
    gsub("sumEmergingPfas", "Σ Short-Chain", .) %>%
    gsub("sumLegacyPfas", "Σ Long-Chain", .)
  
  return(cleaned_names)
}

# Clean Terms ----------
clean_questionnaire_names <- function(term_names) {
  cleaned_terms <- term_names %>%
    # Replace underscores with spaces
    gsub("_", " ", .) %>%
    # Capitalize abbreviations and standardize labels
    gsub("age at enrollment", "Age at Enrollment", ., ignore.case = TRUE) %>%
    gsub("avg drink per day modified", "Alcohol per day: ", ., ignore.case = TRUE) %>%
    gsub("race eth label", "Race/Ethnicity", ., ignore.case = TRUE) %>%
    gsub("sq self hep b", "History of Hepatitis B", ., ignore.case = TRUE) %>%
    gsub("sq self hep c", "History of Hepatitis C", ., ignore.case = TRUE) %>%
    gsub("rural", "Rural vs. Urban", ., ignore.case = TRUE) %>%
    gsub("sex", "Sex", .) %>%
    gsub("smoking", "Smoking/Vaping Status (Yes/No)", ., ignore.case = TRUE) %>%
    gsub("source", "Study Site", ., ignore.case = TRUE) %>%
    gsub("supp meds steroids", "Steroid Medication Use", ., ignore.case = TRUE) %>%
    gsub("sq water bottled",          "Water habits: Uses Bottled Water", ., ignore.case = TRUE) %>%
    gsub("sq water charcoal filter",  "Water habits: Uses Charcoal Filter", ., ignore.case = TRUE) %>%
    gsub("sq water faucet filter",    "Water habits: Uses Faucet Filter", ., ignore.case = TRUE) %>%
    gsub("sq water house filtration", "Water habits: Uses House Filtration System", ., ignore.case = TRUE) %>%
    gsub("sq water none",             "Water habits: Does not Drink Water", ., ignore.case = TRUE) %>%
    gsub("sq water other type",       "Water habits: Uses Other Water Type", ., ignore.case = TRUE) %>%
    gsub("sq water tap unfiltered",   "Water habits: Uses Unfiltered Tap Water", ., ignore.case = TRUE) %>%
    gsub("sq water well",             "Water habits: Uses Well Water", ., ignore.case = TRUE) %>%
    gsub("sq drink alcohol",   "Alcohol Intake", ., ignore.case = TRUE) 
  
  return(cleaned_terms)
}



# Clean Terms ----------
clean_model_terms <- function(term_names) {
  cleaned_terms <- term_names %>%
    # Replace underscores with spaces
    gsub("_", " ", .) %>%
    # Capitalize abbreviations and standardize labels
    gsub("avg drink per day modified", "Alcohol per day: ", ., ignore.case = TRUE) %>%
    gsub("race eth label", "Race/Ethnicity: ", ., ignore.case = TRUE) %>%
    gsub("sq self hep b", "Self-Reported Hepatitis B: ", ., ignore.case = TRUE) %>%
    gsub("sq self hep c", "Self-Reported Hepatitis C: ", ., ignore.case = TRUE) %>%
    gsub("rural", "Rural: ", ., ignore.case = TRUE) %>%
    gsub("sexMale", "Sex: Male", .) %>%
    gsub("smoking", "Smoking Status: ", ., ignore.case = TRUE) %>%
    gsub("source", "Source: ", ., ignore.case = TRUE) %>%
    gsub("supp meds steroids", "Steroid Medications: ", ., ignore.case = TRUE) %>%
    gsub("sq water bottled", "Uses Bottled Water: ", ., ignore.case = TRUE) %>%
    gsub("sq water charcoal filter", "Uses Charcoal Filter: ", ., ignore.case = TRUE) %>%
    gsub("sq water faucet filter", "Uses Faucet Filter: ", ., ignore.case = TRUE) %>%
    gsub("sq water house filtration", "House Filtration System: ", ., ignore.case = TRUE) %>%
    gsub("sq water none", "No Filter Used: ", ., ignore.case = TRUE) %>%
    gsub("sq water other type", "Other Water Type: ", ., ignore.case = TRUE) %>%
    gsub("sq water tap unfiltered", "Uses Unfiltered Tap Water: ", ., ignore.case = TRUE) %>%
    gsub("sq water well", "Uses Well Water: ", ., ignore.case = TRUE) %>%
    # Clean categorical responses
    gsub("Unknown/Not Reported", "(Unknown)", .) %>%
    gsub("Yes", "(Yes)", .) %>%
    gsub("Living in rural area", "Lives in Rural Area", .) %>%
    gsub("Smoke or use vape", "Smokes or Vapes", .) %>%
    gsub("NHB", "Non-Hispanic Black", .) %>%
    gsub("NHO", "Non-Hispanic Other", .) %>%
    gsub("Hispanic", "Hispanic", .) %>% 
    gsub("age at enrollment", "Age at enrollment", .)
  
  return(cleaned_terms)
}
