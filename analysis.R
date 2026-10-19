#=========================================================================
# 01_analysis.R
#
# Menstrual Cycle Status and the Acute Subjective Psychedelic Experience
# Reproduces all analyses reported in the manuscript and supplement.
#
# Input: data/menstruation_psychedelics_minimal.csv (n = 169)
# Run from the repository root.
#=========================================================================

library(dplyr)
library(broom)
library(effectsize)
library(psych)
library(car)
library(logistf)

df <- read.csv("data/menstruation_psychedelics_minimal.csv",
               stringsAsFactors = FALSE)
stopifnot(nrow(df) == 169)

# -------------------------------------------------------------------------
# 1. Four-level cycle status (reference = menstruating during experience)
# -------------------------------------------------------------------------
CYCLE_LEVELS <- c("Menstruating during experience",
                  "Menstruating regularly, not during experience",
                  "Not menstruating regularly",
                  "Male participants")

df <- df %>%
  mutate(
    Cycle = case_when(
      Sex_Birth_Follow == 0     ~ 3,
      Cycle_Month_Follow == 0   ~ 2,
      Cycle_Intense_Follow == 1 ~ 0,
      Cycle_Intense_Follow == 0 ~ 1,
      TRUE ~ NA_real_
    ),
    Cycle = factor(Cycle, levels = 0:3, labels = CYCLE_LEVELS)
  )

# -------------------------------------------------------------------------
# 2. Sample characteristics (Table 1)
# -------------------------------------------------------------------------
table(df$Cycle, useNA = "always")
table(df$Pregnant_Follow, df$Cycle, useNA = "always")
table(df$Breastfeeding_Follow, df$Cycle, useNA = "always")
table(df$Dose_Follow, df$Cycle)
table(df$Age_Four_Baseline, df$Cycle)
table(df$Any_Harm_Follow, df$Cycle)

# Substance used for the most intense experience
substances <- c(Psilocybin  = "Psilocybin_Month_Intense_Follow",
                LSD         = "LSD_Month_Intense_Follow",
                DMT         = "DMT_Month_Intense_Follow",
                Mescaline   = "Mescaline_Month_Intense_Follow",
                Ayahuasca   = "Ayahuasca_Month_Intense_Follow",
                Peyote      = "Peyote_Month_Intense_Follow",
                `San Pedro` = "Pedro_Month_Intense_Follow")

df <- df %>%
  mutate(
    n_subst = rowSums(across(all_of(unname(substances))) == 1, na.rm = TRUE),
    Substance_Intense = case_when(
      n_subst > 1  ~ "Multiple",
      n_subst == 0 ~ NA_character_,
      TRUE ~ names(substances)[max.col(across(all_of(unname(substances))) == 1,
                                       ties.method = "first")]
    ),
    Substance_Intense = factor(Substance_Intense,
                               levels = c(names(substances), "Multiple"))
  )

tab_subst <- table(df$Substance_Intense, df$Cycle)
tab_subst <- cbind(tab_subst, Overall = rowSums(tab_subst))
tab_subst
round(100 * prop.table(tab_subst, margin = 2), 1)

# -------------------------------------------------------------------------
# 3. Internal consistency (Cronbach's alpha)
# -------------------------------------------------------------------------
edi_items <- df %>% dplyr::select(EDI_Follow_1:EDI_Follow_8)
meq_items <- df %>% dplyr::select(MEQ_Follow_1:MEQ_Follow_30)
ceq_items <- df %>% dplyr::select(CEQ_Follow_1:CEQ_Follow_26)

round(c(EDI = psych::alpha(edi_items)$total$raw_alpha,
        MEQ = psych::alpha(meq_items)$total$raw_alpha,
        CEQ = psych::alpha(ceq_items)$total$raw_alpha), 2)

subscale_items <- list(
  "MEQ: Mystical"      = paste0("MEQ_Follow_", c(4, 5, 6, 9, 14, 15, 16, 18, 20, 21, 23, 24, 25, 26, 28)),
  "MEQ: Positive Mood" = paste0("MEQ_Follow_", c(2, 8, 12, 17, 27, 30)),
  "MEQ: Transcendence" = paste0("MEQ_Follow_", c(1, 7, 11, 13, 19, 22)),
  "MEQ: Ineffability"  = paste0("MEQ_Follow_", c(3, 10, 29)),
  "CEQ: Fear"          = paste0("CEQ_Follow_", c(4, 7, 14, 21, 26)),
  "CEQ: Grief"         = paste0("CEQ_Follow_", c(2, 6, 9, 11, 23, 25)),
  "CEQ: Physical"      = paste0("CEQ_Follow_", c(3, 5, 15, 17, 18)),
  "CEQ: Insanity"      = paste0("CEQ_Follow_", c(8, 13, 19)),
  "CEQ: Isolation"     = paste0("CEQ_Follow_", c(1, 10, 24)),
  "CEQ: Death"         = paste0("CEQ_Follow_", c(16, 20)),
  "CEQ: Paranoia"      = paste0("CEQ_Follow_", c(12, 22))
)

subscale_alpha_table <- bind_rows(lapply(names(subscale_items), function(label) {
  a <- psych::alpha(df[, subscale_items[[label]]])
  data.frame(subscale = label,
             n_items  = length(subscale_items[[label]]),
             alpha    = round(a$total$raw_alpha, 2))
}))
print(subscale_alpha_table)

# -------------------------------------------------------------------------
# 4. Covariates
# -------------------------------------------------------------------------
df <- df %>%
  mutate(Dose_Follow       = factor(Dose_Follow),
         Age_Four_Baseline = factor(Age_Four_Baseline))

covar_formula <- "Dose_Follow + Age_Four_Baseline"

# Helper: cycle-contrast rows with unstandardized and standardized
# (outcome z-scored) coefficients
clean_lm <- function(model, dv_name, pattern = "^Cycle") {
  raw <- broom::tidy(model) %>% dplyr::filter(grepl(pattern, term))
  std <- effectsize::standardize_parameters(model, ci = 0.95, method = "posthoc") %>%
    as.data.frame() %>%
    dplyr::filter(grepl(pattern, Parameter)) %>%
    dplyr::rename(term = Parameter, beta = Std_Coefficient,
                  ci_low = CI_low, ci_high = CI_high)
  raw %>%
    dplyr::left_join(std, by = "term") %>%
    dplyr::mutate(outcome = dv_name) %>%
    dplyr::select(outcome, term, estimate, std.error, statistic, p.value,
                  beta, ci_low, ci_high)
}

# -------------------------------------------------------------------------
# 5. Primary models: EDI, MEQ, CEQ totals (Table 2)
# -------------------------------------------------------------------------
primary_outcomes <- c("EDI Total" = "EDI_Follow_Total",
                      "MEQ Total" = "MEQ_Follow_Total",
                      "CEQ Total" = "CEQ_Follow_Total")

m_primary <- lapply(primary_outcomes, function(y)
  lm(as.formula(paste(y, "~ Cycle +", covar_formula)), data = df))

car::vif(m_primary[["EDI Total"]])   # collinearity check (same for all three)

tab_primary <- bind_rows(lapply(names(m_primary), function(nm)
  clean_lm(m_primary[[nm]], nm)))
print(tab_primary, n = Inf)

# -------------------------------------------------------------------------
# 6. Risk of harm: logistic regression + Firth sensitivity (Table S3)
# -------------------------------------------------------------------------
table(df$Cycle, df$Any_Harm_Follow)
round(prop.table(table(df$Cycle, df$Any_Harm_Follow), margin = 1), 2)

harm_formula <- as.formula(paste("Any_Harm_Follow ~ Cycle +", covar_formula))

m_harm <- glm(harm_formula, data = df, family = binomial)
harm_tab <- broom::tidy(m_harm, conf.int = TRUE) %>%
  dplyr::filter(grepl("^Cycle", term)) %>%
  dplyr::transmute(comparison = sub("^Cycle", "During vs. ", term),
                   OR    = exp(-estimate),
                   lower = exp(-conf.high),
                   upper = exp(-conf.low),
                   p     = p.value)
print(harm_tab, digits = 3)
car::Anova(m_harm, type = "II")

m_harm_firth <- logistf(harm_formula, data = df)
firth_tab <- data.frame(
  OR    = exp(-coef(m_harm_firth)),
  lower = exp(-m_harm_firth$ci.upper),
  upper = exp(-m_harm_firth$ci.lower),
  p     = m_harm_firth$prob
)
firth_tab <- firth_tab[grepl("^Cycle", rownames(firth_tab)), ]
rownames(firth_tab) <- sub("^Cycle", "During vs. ", rownames(firth_tab))
round(firth_tab, 3)

# -------------------------------------------------------------------------
# 7. Exploratory subscale models with FDR correction (Table S1)
#    FDR applied across the 11 subscales, separately for each contrast
# -------------------------------------------------------------------------
subscale_outcomes <- c(
  "MEQ: Mystical"      = "MEQ_Mystical_Follow",
  "MEQ: Positive Mood" = "MEQ_Positive_Follow",
  "MEQ: Transcendence" = "MEQ_Transcendence_Follow",
  "MEQ: Ineffability"  = "MEQ_Ineffability_Follow",
  "CEQ: Fear"          = "CEQ_Fear_Follow",
  "CEQ: Grief"         = "CEQ_Grief_Follow",
  "CEQ: Physical"      = "CEQ_Physical_Follow",
  "CEQ: Insanity"      = "CEQ_Insanity_Follow",
  "CEQ: Isolation"     = "CEQ_Isolation_Follow",
  "CEQ: Death"         = "CEQ_Death_Follow",
  "CEQ: Paranoia"      = "CEQ_Paranoia_Follow"
)

m_sub <- lapply(subscale_outcomes, function(y)
  lm(as.formula(paste(y, "~ Cycle +", covar_formula)), data = df))

tab_sub <- bind_rows(lapply(names(m_sub), function(nm)
  clean_lm(m_sub[[nm]], nm))) %>%
  dplyr::group_by(term) %>%
  dplyr::mutate(p.fdr = p.adjust(p.value, method = "fdr")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(outcome, term)
print(tab_sub, n = Inf)

# -------------------------------------------------------------------------
# 8. Sensitivity analysis: binary cycle status (Table S2)
#    Menstruating during experience vs. all other participants
# -------------------------------------------------------------------------
df <- df %>%
  mutate(Cycle_Bin = factor(
    ifelse(Cycle == CYCLE_LEVELS[1], "During", "All others"),
    levels = c("During", "All others")))

m_bin <- lapply(primary_outcomes, function(y)
  lm(as.formula(paste(y, "~ Cycle_Bin +", covar_formula)), data = df))

tab_bin <- bind_rows(lapply(names(m_bin), function(nm)
  clean_lm(m_bin[[nm]], nm, pattern = "^Cycle_Bin")))
print(tab_bin, n = Inf)

# -------------------------------------------------------------------------
sessionInfo()

