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
