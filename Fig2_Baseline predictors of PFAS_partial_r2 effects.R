# Load necessary libraries
library(dplyr)
library(broom)
library(tidyr)
library(ggplot2)
library(forcats)
library(cowplot)
library(ggh4x)
library(circlize)
library(emmeans)
library(ComplexHeatmap)

source(here::here("!directories.r"))
source(here::here("!functions.r"))
source(here::here("!load_cleaned_data.r"))

short_chain_color = "#1B9E77"
long_chain_color = "#D95F02"

# data$x9cl_pf3ons <- replace_na(data$x9cl_pf3ons,  0.1/sqrt(2))

potential_conf <- c('age_at_enrollment', 'sex', 'slm1_bmi', 'race_eth_label', 'rural', 'source',
'smoking', 'sq_drink_alcohol', 'sq_average_drink_per_day', 's1_meds_pain_reliever_take', 'sq_water_well', 
'sq_water_tap_unfiltered', 'sq_water_house_filtration', 
'sq_water_faucet_filter', 'sq_water_charcoal_filter',
'sq_water_bottled', 'sq_water_none')
  
potential_conf <- setdiff(potential_conf, c("slm1_bmi", "s1_meds_pain_reliever_take",'sq_water_dont_know'))
covar_names <- potential_conf
data$race_eth_label <- relevel(as.factor(data$race_eth_label), ref = "NHW")

# legacy <- c(legacy, "x9cl_pf3ons")

long_chain_pfas <-  legacy
short_chain_pfas <- setdiff(emerging, c("pfba", "pf_hx_a"))

# Identify numeric covariates and scale them
data <- data %>%
  mutate(across(all_of(covar_names), ~ if(is.numeric(.)) scale(.) |> as.numeric() else ., .names = "scaled_{col}"))

scaled_covar_names <- paste0("scaled_", covar_names)


## A. Finalize covariate names --------

# Calculate the sum of all legacy PFAS and all short_chain_pfas PFAS
data <- data %>%
  mutate(sum_legacy_pfas = rowSums(dplyr::select(., all_of(legacy)), na.rm = TRUE),
         sum_emerging_pfas = rowSums(dplyr::select(., all_of(short_chain_pfas)), na.rm = TRUE), 
         sum_all_pfas = rowSums(dplyr::select(., all_of(c(legacy, short_chain_pfas))), na.rm = TRUE))


# Add the summed PFAS groups to the classification
summed_pfas <- c("sum_legacy_pfas", "sum_emerging_pfas", "sum_all_pfas")

# Create a list of PFAS and their chain lengths
pfas_list <- c(legacy, short_chain_pfas, summed_pfas)
pfas_chain <- case_when(pfas_list %in% long_chain_pfas ~ "Long-chain", 
                        pfas_list %in% summed_pfas ~ "Summed PFAS", 
                        pfas_list %in% short_chain_pfas ~ "Short-chain")

# Create the new variable, defaulting to scaled_sq_drink_alcohol values
data$scaled_drink_combined <- as.character(data$scaled_sq_drink_alcohol)

# Replace "Yes, current drinker" with values from scaled_sq_average_drink_per_day
replace_idx <- data$scaled_sq_drink_alcohol == "Yes, current drinker"
data$scaled_drink_combined[replace_idx] <- as.character(data$scaled_sq_average_drink_per_day[replace_idx])

# Convert to factor with appropriate levels
data$scaled_drink_combined <- factor(data$scaled_drink_combined, levels = c(
  "No, never drinker",
  "No, former drinker (stopped)",
  "Less than 1 alcoholic drink per day",
  "1-2 alcoholic drinks per day",
  "3-4 alcoholic drinks per day",
  "Unknown/Not Reported"
))

# Check the new variable
table(data$scaled_drink_combined)


# Sets specific levels as the reference
data$scaled_sex <- factor(data$scaled_sex, levels = c("Female", "Male"))  
data$scaled_race_eth_label <- factor(data$scaled_race_eth_label)
data$scaled_rural <- factor(data$scaled_rural, levels = c("Non-Metro Counties", "Metro Counties","Unknown/Not Reported"))
data$scaled_source <- factor(data$scaled_source, levels= c("emory", "duke", "ncsu", "unc"))
data$scaled_smoking <- factor(data$scaled_smoking, levels=c("Not at all", "Everyday", "Some days", "Unknow/Not Reported"))
data$scaled_sq_drink_alcohol <- factor(data$scaled_sq_drink_alcohol, levels =c("No, never drinker", "No, former drinker (stopped)", "Yes, current drinker", "Unknown/Not Reported"))
data$scaled_sq_average_drink_per_day <- factor(data$scaled_sq_average_drink_per_day, levels=c("Less than 1 alcoholic drink per day","1-2 alcoholic drinks per day","3-4 alcoholic drinks per day","Unknown/Not Reported"))
data$scaled_drink_combined <- factor(data$scaled_drink_combined, levels=c("No, never drinker", "No, former drinker (stopped)","Less than 1 alcoholic drink per day","1-2 alcoholic drinks per day","3-4 alcoholic drinks per day","Unknown/Not Reported"))
data$scaled_sq_water_well <- factor(data$scaled_sq_water_well)
data$scaled_sq_water_tap_unfiltered <- factor(data$scaled_sq_water_tap_unfiltered)
data$scaled_sq_water_house_filtration <- factor(data$scaled_sq_water_house_filtration)
data$scaled_sq_water_faucet_filter <- factor(data$scaled_sq_water_faucet_filter)
data$scaled_sq_water_charcoal_filter <- factor(data$scaled_sq_water_charcoal_filter)
data$scaled_sq_water_bottled <- factor(data$scaled_sq_water_bottled)
data$scaled_sq_water_none <- factor(data$scaled_sq_water_none)

# Remove the two variables from the covariate list
scaled_covar_names <- setdiff(scaled_covar_names, c("scaled_sq_drink_alcohol", "scaled_sq_average_drink_per_day"))

# Optionally, add the new combined variable if needed
scaled_covar_names <- union(scaled_covar_names, "scaled_drink_combined")

# Check the updated list
print(scaled_covar_names)


# Initialize an empty dataframe to store results for each regression
results_df <- data.frame()

# Loop over each PFAS in the list and perform a linear regression ------
# pfas <- "pf_hp_a"
for (pfas in pfas_list) {
  # Define the regression formula with all covariates
  formula_str <- paste("scale(log2(", pfas, ")) ~ -1 + ", paste(scaled_covar_names, collapse = " + "))
  
  # Fit the linear model
  model <- lm(as.formula(formula_str), data = data)
  partial_r2 <- rsq::rsq.partial(model, type = "n")
  
  
  resout <- data.frame(variable = partial_r2$variable, 
                       partial_r2 = partial_r2$partial.rsq) |> 
    tidylog::left_join(tidy(car::Anova(model)) |> rename("variable" = "term"))
  
  
  # Extract results and format them
  model_results <- resout |>
    # filter(term != "(Intercept)") %>%  # Exclude the intercept
    mutate(pfas = pfas,
           chain_length = case_when(pfas %in% long_chain_pfas ~ "Long-chain", 
                                    pfas %in% summed_pfas ~ "Summed PFAS", 
                                    pfas %in% short_chain_pfas ~ "Short-chain"),
           R2 = summary(model)$r.squared,
           adj_R2 = summary(model)$adj.r.squared 
    )  # Add chain length
  
  # Combine results into the main results dataframe
  results_df <- bind_rows(results_df, model_results)
}

####
for (pfas in pfas_list) {
  # Define the regression formula with all covariates
  formula_str <- paste("scale(log2(", pfas, ")) ~ -1 + ", paste(scaled_covar_names, collapse = " + "))
  
  # Fit the linear model
  model <- lm(as.formula(formula_str), data = data)
  partial_r2 <- rsq::rsq.partial(model, type = "n")
  
  # Extract beta coefficients, SE, and confidence intervals
  coef_results <- broom::tidy(model, conf.int = TRUE) |> 
    rename(variable = term) # Rename term to variable to match other outputs
  
  # Merge with partial R-squared results
  resout <- data.frame(variable = partial_r2$variable, 
                       partial_r2 = partial_r2$partial.rsq) |> 
    tidylog::left_join(coef_results)
  
  # Add additional metadata and format results
  model_results <- resout |>
    mutate(pfas = pfas,
           chain_length = case_when(pfas %in% long_chain_pfas ~ "Long-chain", 
                                    pfas %in% summed_pfas ~ "Summed PFAS", 
                                    pfas %in% short_chain_pfas ~ "Short-chain"),
           R2 = summary(model)$r.squared,
           adj_R2 = summary(model)$adj.r.squared) 
  
  # Combine results into the main results dataframe
  results_df <- bind_rows(results_df, model_results)
}
write.csv(results_df, fs::path(dir_result, "results_df.csv"))

coef_df <- data.frame()  # Initialize an empty dataframe for coefficient results

for (pfas in pfas_list) {
  # Define the regression formula with all covariates
  formula_str <- paste("scale(log2(", pfas, ")) ~ -1 + ", paste(scaled_covar_names, collapse = " + "))
  
  # Fit the linear model
  model <- lm(as.formula(formula_str), data = data)
  
  # Extract coefficients, SE, and confidence intervals
  coef_results <- broom::tidy(model, conf.int = TRUE) |> 
    rename(variable = term) |> 
    mutate(pfas = pfas)  # Add PFAS name to each coefficient
  
  # Append results to coef_df
  coef_df <- bind_rows(coef_df, coef_results)
}

# Save the coefficient results separately
write.csv(coef_df, fs::path(dir_result, "coef_df.csv"), row.names = FALSE)

###
# Prepare heatmap data --------
# Create annotation DF
r2_values <- results_df %>%
  group_by(pfas) %>%
  slice_head() |>
  ungroup() |>
  tidylog::select(pfas, R2, adj_R2)

#naming
# Define a mapping of original PFAS names to cleaned names
pfas_name_map <- c(
  "pf_hx_s" = "PFHxS",
  "pfda" = "PFDA",
  "pfna" = "PFNA",
  "pfos" = "PFOS",
  "pf_hp_a" = "PFHpA",
  "pfoa" = "PFOA",
  "pf_un_a" = "PFUnA",
  "pf_hp_s" = "PFHpS",
  "pf_do_a" = "PFDoA",
  "pfbs" = "PFBS",
  "pf_pe_a" = "PFPeA",
  "pf_pe_s" = "PFPeS",
  "sum_legacy_pfas" = "Summed Legacy PFAS",
  "sum_emerging_pfas" = "Summed Emerging PFAS",
  "sum_all_pfas" = "Summed PFAS"
)


# Create a lookup table to assign long vs. short chain annotations
chain_annotations <- data.frame(
  pfas = pfas_list,
  Chain = case_when(pfas_list %in% long_chain_pfas ~ "Long-chain", 
                    pfas_list %in% summed_pfas ~ "Summed PFAS", 
                    pfas_list %in% short_chain_pfas ~ "Short-chain")) |>
  mutate(
    pfas_clean_names = clean_pfas_names(pfas),
    mol_weight = case_when(
      pfas_clean_names == "PFDA"    ~ 514.08,
      pfas_clean_names == "PFOS"    ~ 500.13,
      pfas_clean_names == "PFHxS"   ~ 400.12,
      pfas_clean_names == "PFNA"    ~ 464.08,
      pfas_clean_names == "PFBS"    ~ 300.10,
      pfas_clean_names == "PFHpA"   ~ 364.06,
      pfas_clean_names == "PFOA"    ~ 414.07,
      pfas_clean_names == "PFPeA"   ~ 264.05,
      pfas_clean_names == "PFUnA"   ~ 564.09,
      pfas_clean_names == "PFHpS"   ~ 450.11,
      pfas_clean_names == "PFDoA"   ~ 614.10,
      pfas_clean_names == "PFPeS"   ~ 350.10,
      pfas_clean_names == "PFHxA"   ~ 314.05, 
      pfas_clean_names == "PFBA"   ~ 214.04, 
      pfas_clean_names == "9Cl-PF3ONS"   ~ 570.67, 
      TRUE ~ NA_real_), 
    functional_group = case_when(str_detect(pfas, "Sum") ~ "Summed PFAS", 
                                 grepl("PF3ONS$", pfas) ~ "PFESA", 
                                 grepl("s$", pfas) ~ "Sulfonic Acid", 
                                 grepl("a$", pfas) ~ "Carboxylic Acid"))

# Combine R² values and chain length annotations for PFAS columns and reorder to match heatmap columns
pfas_annotations <- r2_values %>%
  tidylog::left_join(chain_annotations, by = "pfas") %>%
  mutate(pfas = clean_pfas_names(pfas_list), 
         Type = if_else(str_detect(Chain, "Summed"), "Summed PFAS", "Individual PFAS"), 
         Chain2 = case_when(str_detect(pfas_clean_names, "Long-Chain") ~ "Long-chain", 
                            str_detect(pfas_clean_names, "Short-Chain") ~ "Short-chain", 
                            TRUE ~ Chain)) |>
  arrange(pfas) |>
  column_to_rownames(var = "pfas")


# Prepare the data matrix for the heatmap
heatmap_data <- results_df %>%
  filter(!is.na(df)) |>
  dplyr::select(pfas, variable, partial_r2) %>%
  tidylog::pivot_wider(names_from = pfas, values_from = partial_r2) %>%
  column_to_rownames("variable") |>
  rename_with(.fn = clean_pfas_names) |>
  as.matrix()


# Prepare significance matrix (TRUE for significant, FALSE otherwise)
sig_matrix <- results_df %>%
  filter(!is.na(df)) |>
  mutate(significant = if_else(p.value < 0.05, "*", "")) %>%
  dplyr::select(pfas, variable, significant) %>%
  tidylog::pivot_wider(names_from = pfas, values_from = significant) %>%
  column_to_rownames("variable") |>
  rename_with(.fn = clean_pfas_names) |>
  as.matrix()

# Clean row names for heatmap
rownames(heatmap_data) <- str_remove(rownames(sig_matrix), "scaled_") |> 
  clean_questionnaire_names()

# Ensure that pfas_annotations and heatmap_data columns are in the same order
pfas_annotations <- pfas_annotations[colnames(heatmap_data), ]

# Ensure no negative values in heatmap_data for display purposes
heatmap_data[heatmap_data < 0] <- 0

### Column Annotations --------
# the annotations for chain length and R² values
col_annotations <- HeatmapAnnotation(
  `Molecular Weight`  = pfas_annotations$mol_weight,
  `Overall Adjusted R²` = pfas_annotations$adj_R2,
  `PFAS Category` = pfas_annotations$Chain2,
  col = list(
    `PFAS Category` = c("Long-chain" = long_chain_color, "Short-chain" = short_chain_color, "Summed PFAS" = "gray"),
    `Molecular Weight` = colorRamp2(c(min(pfas_annotations$mol_weight, na.rm = TRUE),
                                           max(pfas_annotations$mol_weight, na.rm = TRUE)),
                                         c("white", "blue")), 
    `Overall Adjusted R²` = colorRamp2(c(min(pfas_annotations$adj_R2, na.rm = TRUE),
                                          max(pfas_annotations$adj_R2, na.rm = TRUE)),
                                        c("white", "red")))
  # show_legend = FALSE
)


# Define color scale for the heatmap
heatmap_colors <- colorRamp2(c(0, max(heatmap_data, na.rm = TRUE)), c("white", "#E1BE6A"))
# heatmap_colors <- colorRamp2(c(0, max(pfas_annotations$R2, na.rm = TRUE)), c("white", "#E1BE6A"))


## Create the heatmap --------

ht1 <- Heatmap(heatmap_data,
               name = "Partial R²",
               col = heatmap_colors,
               row_names_side = "right",
               cluster_rows = TRUE,
               cluster_row_slices = FALSE, 
               cluster_columns = TRUE,  
               cluster_column_slices = FALSE,
               column_split = pfas_annotations$Type,  # Split columns by chain length
               # bottom_annotation = col_annotations[,c(1,3)],  # Add chain length and R2 as column annotations
               top_annotation = col_annotations[,c(3,1,2)],  # Add chain length and R2 as column annotations
               show_column_names = TRUE,
               show_row_names = TRUE,
               row_dend_side = "left",
               column_title = " ",
               row_title = "Baseline Lifestyle/Demographic Factors",
               row_names_max_width = unit(15, "cm"),  # Set width to prevent overlap
               cell_fun = function(j, i, x, y, width, height, fill) {
                 if (sig_matrix[i, j] == "*") {
                   # grid.text("*", x, y, gp = gpar(fontsize = 20, col = "black"))
                   grid.text("*", x, y - unit(1.5, "mm"), gp = gpar(fontsize = 20, col = "black"))
                   
                 }
               }
)

ht1
# save heatmap 
png(fs::path(dir_figures, "Heatmap_adjusted_R2_joint_models.png"),
    width = 1300, height = 1000, res = 150)  # Adjust resolution and size as needed
(draw(ht1, padding = unit(c(1, 3, 0, 10), "mm"), merge_legends = TRUE))
dev.off()

# Plot overall model R2 vs. chain type ---------
# Test for differences by molecular weight- not sig
pfas_annotations %>%
  filter(!is.na(mol_weight),
         # pfas != "PFBA",
         pfas != "9Cl-PF3ONS") %>%
  lm(R2 ~ mol_weight, data = .) |>
  summary()

pfas_annotations |> 
  group_by(Chain) |> 
  summarise(mean_adj_r2=mean(adj_R2), 
            sd=sd(adj_R2))

# Test for differences by PFAS group- significant
temp <- pfas_annotations %>%
  filter(!is.na(mol_weight))  

# Run test
(ttestp = t.test(temp$adj_R2 ~ temp$Chain)$p.value)


## Boxplot of R2 vs. chainlength ----
data_for_boxplot <- pfas_annotations |> 
  mutate(pfas_clean_names = str_remove(pfas_clean_names, " PFAS"), 
         summed_var = if_else(str_detect(Chain, "Summed"), "Summed", "Individual"))

yadj = 0.0
xadj = 0

(boxplot_of_r2 <- data_for_boxplot %>%
    filter(str_detect(pfas_clean_names, "all", negate = TRUE)) |>
    ggplot(aes(x = Chain2, y = adj_R2)) +
    # Boxplot: Remove summed PFAS for the boxplot
    geom_boxplot(
      data = data_for_boxplot %>% filter(Type != "Summed PFAS"),
      aes(fill = Chain2),
      position = position_nudge(x = -0.1-xadj),
      width = 0.07, 
      show.legend = FALSE) +  # Shift boxplot slightly to the left
    # Add labels on the right side of each point
    ggrepel::geom_label_repel(aes(label = pfas_clean_names), 
                              # position = position_nudge(x = -xadj), 
                              nudge_x = 0.15-xadj,  # Nudges labels to the right
                              # xlim=c(1+.3-xadj, 2+.3-xadj),
                              hjust = 0,  # Aligns labels to the left side
                              min.segment.length = unit(0, 'lines'),
                              direction = "y", 
                              max.overlaps = 10) +  # Keeps labels aligned horizontally
    geom_point(aes(fill = Chain2, size = summed_var, shape = summed_var), 
               color = "black", show.legend = FALSE, 
               position = position_nudge(x = -xadj)) +  
    # Add a line above the two boxplots
    annotate("segment", x = 1 -xadj, xend = 2 - xadj, y = yadj+0.35, yend = yadj+0.35, color = "black") +
    # Add small vertical ticks at the ends of the line
    annotate("segment", x = 1 -xadj, xend = 1 - xadj, y = yadj+0.34, yend = yadj+0.35, color = "black") +  # Left tick
    annotate("segment", x = 2 -xadj, xend = 2 - xadj, y = yadj+0.34, yend = yadj+0.35, color = "black") +  # Right tick
    # Add p-value annotation above the line
    annotate("text",    x = 1.5 -xadj, y = yadj+0.358, label = paste0("P = ",round(ttestp, 3))) +
    scale_size_manual(breaks = c("Summed", "Individual"), values = c(3, 2)) +
    scale_shape_manual(breaks = c("Summed", "Individual"), values = c(23, 21)) +
    scale_color_manual(breaks = c("Long-chain", "Short-chain"), 
                       values = c(long_chain_color, short_chain_color)) +
    scale_fill_manual(breaks = c("Long-chain", "Short-chain"), 
                       values = c(long_chain_color, short_chain_color)) +
    ylab(expression("Adjusted R"^2 ~ "for overall model")) +
    ylim(c(0, .45)) + 
    # xlab("Long vs. Short-chain") +
    xlab(NULL) +
    theme_cowplot() + 
    theme(axis.text.x = element_text(angle = 90, vjust = .5))
)

# Save the boxplot
ggsave(fs::path(dir_figures, "boxplot_of_adjusted_r2.jpeg"),
       boxplot_of_r2, width = 5.5, height = 7, dpi = 300)

## Scatter plot of R2 vs. mol weight ----
# pfas_annotations |>
#   filter(!is.na(mol_weight)) |>
#   ggplot( aes(x = mol_weight, y = R2)) +
#   geom_point(aes(color = Chain)) +
#   geom_smooth(method = "lm") +
#   ggrepel::geom_label_repel(aes(label = pfas_clean_names),
#                             min.segment.length = unit(0, 'lines')) +
#   theme_cowplot()

table(data$smoking)




