# =============================================================================
# ELSA DATA PREPROCESSING FOR CAUSAL MEDIATION ANALYSIS
# =============================================================================
# Project: Bayesian Causal Mediation Analysis with Ordinal Outcomes
# Data: English Longitudinal Study of Ageing (ELSA) Waves 2-7
# 
# FUNDING AND ACKNOWLEDGEMENTS
# This research was funded by the Economic and Social Research Council (ESRC) 
# through grant ES/W012588/1 (Investigators: C. R. Victor, J. Fisher, K. Yu, 
# K. Smith, M. Thomas). ESRC is part of UK Research and Innovation (UKRI). 
# The views expressed in this publication are those of the authors and do not 
# necessarily represent those of the ESRC or UKRI.
#
# ELSA is supported by the National Institute on Aging of the National 
# Institutes of Health under Award Number R01AG017644 and by the National 
# Institute for Health Research (NIHR) Policy Research Programme 
# (RP-PG-0310-1074_03). The views expressed are those of the authors and do 
# not necessarily represent those of the National Institutes of Health, the 
# NIHR, or the Department of Health and Social Care.
#
# DATA AVAILABILITY
# ELSA data are available through the UK Data Service:
# https://beta.ukdataservice.ac.uk/datacatalogue/series/series?id=200011
#
# CODE ATTRIBUTION
# Portions of this preprocessing script are adapted from:
# Mayerl, H., Stolz, E., & Freidl, W. (2023). Lonely and depressed in older 
# age: prospective associations and common vulnerabilities. Aging & Mental 
# Health, 27(3), 640-645. https://doi.org/10.1080/13607863.2022.2056138
# Original code available at: https://osf.io/jrhq7/overview
# =============================================================================

# Load required packages
library(foreign)
library(car)
library(dplyr)
library(tidyr)
library(mice)
library(stringr)

# Custom functions for scale calculation
mymean <- function(x, y=3) {
  if (sum(!is.na(x)) >= y) mean(x, na.rm = TRUE) else NA
}

mysum <- function(x, y=3) {
  if (sum(!is.na(x)) >= y) sum(x, na.rm = TRUE) else NA
}

# =============================================================================
# DATA IMPORT
# =============================================================================

# Import ELSA core data files (save as RDS for faster subsequent loading)
elsa_2 <- read.spss("./wave_2_core_data_v4.sav", to.data.frame = TRUE)
elsa_2_2 <- read.spss("./wave_2_ifs_derived_variables.sav", to.data.frame = TRUE)
elsa_3 <- read.spss("./wave_3_elsa_data_v4.sav", to.data.frame = TRUE)
elsa_4 <- read.spss("./wave_4_elsa_data_v3.sav", to.data.frame = TRUE)
elsa_5 <- read.spss("./wave_5_elsa_data_v4.sav", to.data.frame = TRUE)
elsa_6 <- read.spss("./wave_6_elsa_data_v2.sav", to.data.frame = TRUE)
elsa_7 <- read.spss("./wave_7_elsa_data.sav", to.data.frame = TRUE)

saveRDS(elsa_2, "elsa_2.rds")
saveRDS(elsa_2_2, "elsa_2_2.rds")
saveRDS(elsa_3, "elsa_3.rds")
saveRDS(elsa_4, "elsa_4.rds")
saveRDS(elsa_5, "elsa_5.rds")
saveRDS(elsa_6, "elsa_6.rds")
saveRDS(elsa_7, "elsa_7.rds")

# =============================================================================
# BASELINE DEMOGRAPHIC VARIABLES (WAVE 2)
# =============================================================================

# Age categories: 50-59 (0), 60-69 (1), 70-79 (2), 80+ (3)
elsa_2$age_y <- as.numeric(as.character(elsa_2$dhager))
elsa_3$age_y <- as.numeric(as.character(elsa_3$dhager))
elsa_2$age_gr <- cut(elsa_2$age_y, breaks = c(0, 50, 60, 70, 80, 91), right = FALSE)
elsa_3$age_gr <- cut(elsa_3$age_y, breaks = c(0, 50, 60, 70, 80, 91), right = FALSE)

# Gender
elsa_2$dhsex2 <- factor(elsa_2$DhSex)
elsa_3$dhsex2 <- factor(elsa_3$dhsex)

# Cohabitation status (descriptive purposes only)
elsa_2$couple2 <- factor(
  ifelse(elsa_2$couple == "Married" | elsa_2$couple == "Cohabit", 1, 0),
  labels = c("Single", "Cohabit"))
elsa_3$couple2 <- factor(
  ifelse(elsa_3$couple == "Married" | elsa_3$couple == "Cohabit", 1, 0),
  labels = c("Single", "Cohabit"))

# Marital status
elsa_2$marital <- car::recode(elsa_2$DiMar, 
                              "'Single, that is never married'='single'; 
   'Married, first and only marriage'='married';
   'Remarried, second or later marriage'='married'; 
   'Legally separated'='legallyseparated';
   'Divorced'='divorced'; 
   'Widowed'='widowed'; 
   else=NA")

elsa_2$defacto_marital <- factor(
  ifelse(elsa_2$couple2 == "Cohabit" | elsa_2$marital == "married", 
         "partnered", "not_partnered"),
  labels = c("not_partnered", "partnered"))

# Treatment variable: Living alone (single-person household)
elsa_2$livalone <- with(elsa_2, factor(
  ifelse(as.numeric(DhR) %in% 4:25 | as.numeric(DhR2) %in% 4:25 | 
           as.numeric(DhR3) %in% 4:25 | as.numeric(DhR4) %in% 4:25 | 
           as.numeric(DhR5) %in% 4:25 | as.numeric(DhR6) %in% 4:25 | 
           as.numeric(DhR7) %in% 4:25 | as.numeric(DhR8) %in% 4:25 | 
           as.numeric(DhR9) %in% 4:25 | as.numeric(DhR10) %in% 4:25 | 
           as.numeric(DhR11) %in% 4:25 | as.numeric(DhR12) %in% 4:25, 0, 1), 
  labels = c("no", "yes")))

# Educational qualification
elsa_2 <- merge(elsa_2, elsa_2_2[, c("idauniq", "edqual")], by="idauniq", all.x=TRUE)
elsa_2$edqual2 <- car::recode(elsa_2$edqual, 
                              "'nvq4/nvq5/degree or equiv'='high'; 
   'higher ed below degree'='high';
   'nvq3/gce a level equiv'='low'; 
   'nvq2/gce o level equiv'='low';
   'nvq1/cse other grade equiv'='low'; 
   'no qualification'='low'; 
   else=NA")

# Life satisfaction scale
names_sclife <- c("sclifea", "sclifeb", "sclifec", "sclifed", "sclifee")
elsa_2[, names_sclife] <- apply(elsa_2[, names_sclife], 2, function(x) {
  car::recode(x, "'Strongly disagree'=1; 'Disagree'=2; 'Slightly disagree'=3; 
               'Neither agree nor disagree'=4; 'Slightly agree'=5; 'Agree'=6; 
               'Strongly agree'=7; else=NA")
})
elsa_2$sclife_score_w2 <- apply(elsa_2[, names_sclife], 1, function(x) mysum(x, 5))

# Restrict to participants aged 50+ at baseline
elsa_2 <- subset(elsa_2, age_y >= 50)
elsa_3 <- elsa_3[elsa_3$idauniq %in% elsa_2$idauniq, ]
elsa_4 <- elsa_4[elsa_4$idauniq %in% elsa_2$idauniq, ]
elsa_5 <- elsa_5[elsa_5$idauniq %in% elsa_2$idauniq, ]
elsa_6 <- elsa_6[elsa_6$idauniq %in% elsa_2$idauniq, ]
elsa_7 <- elsa_7[elsa_7$idauniq %in% elsa_2$idauniq, ]

# =============================================================================
# OUTCOME VARIABLE: LONELINESS (UCLA 3-ITEM SCALE, RANGE 3-9)
# =============================================================================

names_lonel <- c("scfeela", "scfeelb", "scfeelc")

# Handle ELSA-specific missing values
handle_missing <- function(x) {
  ifelse(x %in% c(-9, -8, -2, -1), NA, x)
}

# Recode loneliness items (note ELSA coding inconsistencies)
elsa_2[, names_lonel[-3]] <- apply(elsa_2[, names_lonel[-3]], 2, function(x) {
  car::recode(x, "'Hardly ever or never'=1; 'Some of the time'=2; 
               'Often'=3; else=NA")
})
elsa_2$scfeelc <- car::recode(elsa_2$scfeelc, 
                              "'Hardly ever of never'=1; 'Some of the time'=2; 'Often'=3; else=NA", 
                              as.factor = FALSE)

elsa_3[, names_lonel[-3]] <- apply(elsa_3[, names_lonel[-3]], 2, function(x) {
  car::recode(x, "'Hardly ever or never'=1; 'Some of the time'=2; 
               'Often'=3; else=NA")
})
elsa_3$scfeelc <- car::recode(elsa_3$scfeelc, 
                              "'Hardly ever of never'=1; 'Some of the time'=2; 'Often'=3; else=NA", 
                              as.factor = FALSE)

elsa_4[, names_lonel[-3]] <- apply(elsa_4[, names_lonel[-3]], 2, function(x) {
  car::recode(x, "'Hardly ever or never'=1; 'Some of the time'=2; 
               'Often'=3; else=NA")
})
elsa_4$scfeelc <- car::recode(elsa_4$scfeelc, 
                              "'Hardly ever of never'=1; 'Some of the time'=2; 'Often'=3; else=NA", 
                              as.factor = FALSE)

elsa_5[, names_lonel] <- apply(elsa_5[, names_lonel], 2, function(x) {
  car::recode(x, "'Hardly ever or never'=1; 'Some of the time'=2; 
               'Often'=3; else=NA")
})

elsa_6[, names_lonel[-3]] <- apply(elsa_6[, names_lonel[-3]], 2, function(x) {
  car::recode(x, "'Hardly ever or never'=1; 'Some of the time'=2; 
               'Often'=3; else=NA")
})
elsa_6$scfeelc <- car::recode(elsa_6$scfeelc, 
                              "'Hardly ever of never'=1; 'Some of the time'=2; 'Often'=3; else=NA", 
                              as.factor = FALSE)

elsa_7[, names_lonel] <- apply(elsa_7[, names_lonel], 2, function(x) {
  car::recode(x, "'Hardly ever or never'=1; 'Some of the time'=2; 
               'Often'=3; else=NA")
})

# Calculate loneliness sum scores
for(wave in 2:7) {
  wave_data <- get(paste0("elsa_", wave))
  wave_data[[paste0("lonel_score_w", wave)]] <- 
    apply(wave_data[, names_lonel], 1, function(x) mysum(x, 3))
  assign(paste0("elsa_", wave), wave_data)
}

# =============================================================================
# MEDIATOR 1: DEPRESSION (CES-D EXCLUDING LONELINESS ITEM, RANGE 0-7)
# =============================================================================
# Following Mayerl et al. (2023), exclude loneliness item (PScedE) to avoid
# tautological overlap with outcome measure

names_cesd <- tolower(c("PScedA", "PScedB", "PScedC", "PScedG", "PScedH"))
names_cesd_rev <- tolower(c("PScedD", "PScedF"))
names_cesd_all <- c(names_cesd, names_cesd_rev)

# Standardise variable names to lowercase
names(elsa_2)[match(toupper(names_cesd_all), names(elsa_2))] <- names_cesd_all
names(elsa_6)[match(toupper(names_cesd_all), names(elsa_6))] <- names_cesd_all
names(elsa_7)[match(toupper(names_cesd_all), names(elsa_7))] <- names_cesd_all

# Recode CES-D items
for(wave in 2:7) {
  wave_data <- get(paste0("elsa_", wave))
  
  # Regular items
  wave_data[, names_cesd] <- apply(wave_data[, names_cesd], 2, function(x) {
    car::recode(x, "'Yes'=1; 'No'=0; else=NA", as.factor=FALSE)
  })
  
  # Reverse-coded items
  wave_data[, names_cesd_rev] <- apply(wave_data[, names_cesd_rev], 2, function(x) {
    car::recode(x, "'Yes'=0; 'No'=1; else=NA", as.factor=FALSE)
  })
  
  # Sum score
  wave_data[[paste0("cesd_score_w", wave)]] <- 
    apply(wave_data[, names_cesd_all], 1, function(x) mysum(x, 7))
  
  assign(paste0("elsa_", wave), wave_data)
}

# =============================================================================
# MEDIATOR 2: TRANSPORTATION MOBILITY
# =============================================================================

handle_transport_missing <- function(x) {
  valid_responses <- c("Yes", "No", "No, passenger only")
  ifelse(x %in% valid_responses, x, NA)
}

create_transport_mobility <- function(wave_data, car_var, drive_var) {
  car_clean <- handle_transport_missing(as.character(wave_data[[car_var]]))
  drive_clean <- handle_transport_missing(as.character(wave_data[[drive_var]]))
  
  dplyr::case_when(
    car_clean == "Yes" | drive_clean == "Yes" ~ 0,
    car_clean == "No" & (is.na(drive_clean) | drive_clean == "No, passenger only") ~ 1,
    TRUE ~ NA_real_
  )
}

# Apply to all waves (adjusting for variable name case differences)
elsa_2$transport_mobility_w2 <- create_transport_mobility(elsa_2, "SPCar", "SPCarA")
elsa_3$transport_mobility_w3 <- create_transport_mobility(elsa_3, "spcar", "spcara")
elsa_4$transport_mobility_w4 <- create_transport_mobility(elsa_4, "spcar", "spcara")
elsa_5$transport_mobility_w5 <- create_transport_mobility(elsa_5, "spcar", "spcara")
elsa_6$transport_mobility_w6 <- create_transport_mobility(elsa_6, "SPCar", "SPCarA")
elsa_7$transport_mobility_w7 <- create_transport_mobility(elsa_7, "SPCar", "SPCarA")

# =============================================================================
# MEDIATOR 3: MOBILITY LIMITATIONS
# =============================================================================

handle_mobility_missing <- function(x) {
  valid_responses <- c("...no difficulty,", "some difficulty,", "much difficulty,", 
                       "Unable to do this", "No difficulty", "Some difficulty", 
                       "Much difficulty")
  x <- as.character(x)
  ifelse(x %in% valid_responses, x, NA_character_)
}

create_mobility_limitations <- function(wave_data, mobility_var) {
  mobility_clean <- handle_mobility_missing(as.character(wave_data[[mobility_var]]))
  
  dplyr::case_when(
    mobility_clean %in% c("...no difficulty,", "No difficulty") ~ 0,
    mobility_clean %in% c("some difficulty,", "Some difficulty") ~ 1,
    mobility_clean %in% c("much difficulty,", "Much difficulty") ~ 2,
    mobility_clean == "Unable to do this" ~ 3,
    TRUE ~ NA_real_
  )
}

# Apply to all waves
elsa_2$mobility_limitations_w2 <- create_mobility_limitations(elsa_2, "HeFunc")
elsa_3$mobility_limitations_w3 <- create_mobility_limitations(elsa_3, "hefunc")
elsa_4$mobility_limitations_w4 <- create_mobility_limitations(elsa_4, "hefunc")
elsa_5$mobility_limitations_w5 <- create_mobility_limitations(elsa_5, "HeFunc")
elsa_6$mobility_limitations_w6 <- create_mobility_limitations(elsa_6, "HeFunc")
elsa_7$mobility_limitations_w7 <- create_mobility_limitations(elsa_7, "HeFunc")

# =============================================================================
# MERGE DATA ACROSS WAVES
# =============================================================================

# Loneliness outcome
data_lonel2_7 <- Reduce(function(x, y) merge(x, y, by="idauniq", all=TRUE), 
                        list(
                          elsa_2[, c("idauniq", paste0("lonel_score_w", 2), "age_y", "age_gr", 
                                     "dhsex2", "couple2", "marital", "defacto_marital", "livalone", 
                                     "Hehelf", "edqual2", "sclife_score_w2")],
                          elsa_3[, c("idauniq", paste0("lonel_score_w", 3))],
                          elsa_4[, c("idauniq", paste0("lonel_score_w", 4))],
                          elsa_5[, c("idauniq", paste0("lonel_score_w", 5))],
                          elsa_6[, c("idauniq", paste0("lonel_score_w", 6))],
                          elsa_7[, c("idauniq", paste0("lonel_score_w", 7))]
                        ))

# Add mediators
mediator_waves <- lapply(2:7, function(w) {
  wave_data <- get(paste0("elsa_", w))
  wave_data[, c("idauniq", 
                paste0("cesd_score_w", w),
                paste0("transport_mobility_w", w),
                paste0("mobility_limitations_w", w))]
})

data_complete <- Reduce(function(x, y) merge(x, y, by="idauniq", all=TRUE),
                        c(list(data_lonel2_7), mediator_waves))

# Sample selection: at least one observation for loneliness and one for each mediator
names_lonel_score <- paste0("lonel_score_w", 2:7)
names_cesd_score <- paste0("cesd_score_w", 2:7)

elsa_final <- data_complete[rowSums(!is.na(data_complete[, names_lonel_score])) >= 1, ]
elsa_final <- elsa_final[rowSums(!is.na(elsa_final[, names_cesd_score])) >= 1, ]

write.csv(elsa_final, "elsa_waves2_7_complete.csv", row.names = FALSE)

# =============================================================================
# VARIABLE PREPARATION FOR ANALYSIS
# =============================================================================

data <- read.csv("elsa_waves2_7_complete.csv")

# Select analysis variables
selected_vars <- c(
  "idauniq", "age_y", "age_gr", "dhsex2", "couple2", "marital", 
  "defacto_marital", "livalone", "Hehelf", "edqual2", "sclife_score_w2",
  paste0("lonel_score_w", 2:7),
  paste0("cesd_score_w", 2:7),
  paste0("transport_mobility_w", 2:7),
  paste0("mobility_limitations_w", 2:7)
)

analysis_data <- data[, selected_vars[selected_vars %in% names(data)]]

# Convert to ordinal factors
lonel_vars <- paste0("lonel_score_w", 2:7)
analysis_data[lonel_vars] <- lapply(analysis_data[lonel_vars], function(x) {
  factor(pmin(pmax(round(x), 3), 9), levels = 3:9, ordered = TRUE)
})

cesd_vars <- paste0("cesd_score_w", 2:7)
analysis_data[cesd_vars] <- lapply(analysis_data[cesd_vars], function(x) {
  factor(pmin(pmax(round(x), 0), 7), levels = 0:7, ordered = TRUE)
})

# =============================================================================
# MULTIPLE IMPUTATION
# =============================================================================

# Impute education based on age-specific probabilities
if(sum(is.na(analysis_data$edqual2)) > 0) {
  edu_probs <- analysis_data %>%
    filter(!is.na(edqual2)) %>%
    group_by(age_gr) %>%
    summarise(prob_high = mean(edqual2 == "high", na.rm = TRUE), .groups = 'drop')
  
  set.seed(500)
  for(age_group in unique(analysis_data$age_gr[!is.na(analysis_data$age_gr)])) {
    missing_mask <- is.na(analysis_data$edqual2) & analysis_data$age_gr == age_group
    n_missing <- sum(missing_mask)
    
    if(n_missing > 0) {
      prob_high <- edu_probs$prob_high[edu_probs$age_gr == age_group]
      if(length(prob_high) == 0) {
        prob_high <- mean(analysis_data$edqual2 == "high", na.rm = TRUE)
      }
      analysis_data$edqual2[missing_mask] <- sample(
        c("high", "low"), n_missing, prob = c(prob_high, 1 - prob_high), replace = TRUE
      )
    }
  }
}

# Recode baseline variables
analysis_data <- analysis_data %>%
  mutate(
    self_reported_health = dplyr::recode(Hehelf,
                                         "or, poor?" = "poor", "fair," = "poor", 
                                         "excellent," = "good", "very good," = "good", "good," = "good",
                                         .default = NA_character_),
    sclife = factor(ifelse(sclife_score_w2 > median(sclife_score_w2, na.rm = TRUE), 
                           "high", "low")),
    age_gr = case_when(
      age_gr == "[50,60)" ~ 0, age_gr == "[60,70)" ~ 1,
      age_gr == "[70,80)" ~ 2, age_gr == "[80,91)" ~ 3,
      TRUE ~ NA_real_
    ),
    dhsex2 = ifelse(dhsex2 == "Female", 0, 1),
    edqual2 = ifelse(edqual2 == "high", 0, 1),
    self_reported_health = ifelse(self_reported_health == "good", 0, 1),
    sclife = ifelse(sclife == "high", 0, 1),
    livalone = ifelse(livalone == "no", 0, 1)
  )

# Convert transport mobility to numeric
transport_vars <- paste0("transport_mobility_w", 2:7)
for(var in transport_vars) {
  if(var %in% names(analysis_data)) {
    analysis_data[[var]] <- ifelse(analysis_data[[var]] == "High_mobility", 0, 1)
  }
}

# Convert mobility limitations to numeric
mobility_vars <- paste0("mobility_limitations_w", 2:7)
for(var in mobility_vars) {
  if(var %in% names(analysis_data)) {
    analysis_data[[var]] <- case_when(
      analysis_data[[var]] == "No_limitations" ~ 0,
      analysis_data[[var]] == "Some_limitations" ~ 1,
      analysis_data[[var]] == "High_limitations" ~ 2,
      analysis_data[[var]] == "Unable" ~ 3,
      TRUE ~ NA_real_
    )
  }
}

# Convert all variables to integers
base_vars <- c("age_gr", "dhsex2", "livalone", "edqual2", 
               "self_reported_health", "sclife")
all_numeric_vars <- c(base_vars, transport_vars, mobility_vars, lonel_vars, cesd_vars)
existing_vars <- all_numeric_vars[all_numeric_vars %in% names(analysis_data)]
analysis_data[existing_vars] <- lapply(analysis_data[existing_vars], as.integer)

# Remove unnecessary columns
analysis_data <- analysis_data[, !(names(analysis_data) %in% 
                                     c("Hehelf", "sclife_score_w2"))]

# Perform single imputation
imputation_methods <- rep("pmm", ncol(analysis_data))
names(imputation_methods) <- names(analysis_data)

lonel_vars_in_data <- lonel_vars[lonel_vars %in% names(analysis_data)]
cesd_vars_in_data <- cesd_vars[cesd_vars %in% names(analysis_data)]
transport_vars_in_data <- transport_vars[transport_vars %in% names(analysis_data)]
mobility_vars_in_data <- mobility_vars[mobility_vars %in% names(analysis_data)]

imputation_methods[lonel_vars_in_data] <- "polr"
imputation_methods[cesd_vars_in_data] <- "polr"
imputation_methods[transport_vars_in_data] <- "logreg"
imputation_methods[mobility_vars_in_data] <- "polr"

if("edqual2" %in% names(analysis_data)) imputation_methods["edqual2"] <- "logreg"

if(sum(sapply(analysis_data, function(x) sum(is.na(x)))) > 0) {
  imputed <- mice(analysis_data, method = imputation_methods, m = 1, 
                  maxit = 100, seed = 500, printFlag = FALSE)
  completed_data <- complete(imputed, 1)
} else {
  completed_data <- analysis_data
}

write.csv(completed_data, "elsa_imputed_waves2_7.csv", row.names = FALSE)

# =============================================================================
# CONVERT TO LONGITUDINAL FORMAT
# =============================================================================

time_varying_pattern <- paste(c("lonel_score_w", "cesd_score_w", 
                                "transport_mobility_w", "mobility_limitations_w"), 
                              collapse = "|")

elsa_long <- completed_data %>%
  pivot_longer(
    cols = matches(time_varying_pattern),
    names_to = "variable_wave",
    values_to = "value"
  ) %>%
  mutate(
    wave = as.integer(str_extract(variable_wave, "\\d+$")),
    variable_type = case_when(
      str_detect(variable_wave, "lonel_score") ~ "loneliness",
      str_detect(variable_wave, "cesd_score") ~ "depression",
      str_detect(variable_wave, "transport_mobility") ~ "transport_mobility",
      str_detect(variable_wave, "mobility_limitations") ~ "mobility_limitations"
    )
  ) %>%
  select(-variable_wave) %>%
  pivot_wider(names_from = variable_type, values_from = value) %>%
  arrange(idauniq, wave) %>%
  filter(!(is.na(loneliness) & is.na(depression) & 
             is.na(transport_mobility) & is.na(mobility_limitations)))

write.csv(elsa_long, "elsa_longitudinal_analysis.csv", row.names = FALSE)

# =============================================================================
# END OF PREPROCESSING
# =============================================================================