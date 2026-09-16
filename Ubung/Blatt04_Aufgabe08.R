library(tidyverse)
library(visreg)
library(effects)

data <- load("../Daten/wdi.Rdata")

wdi <- wdi %>% filter(year == "2013") %>% # nur Jahr 2013
  filter(!(is.na(CO2emission) | is.na(continent))) %>% # CO2 oder Kontinent müssen vorhanden sein
  mutate(
    CO2emission = log10(CO2emission), # log10 Transformation von CO2 Emission (da Verteilung rechtsschief)
    continent = as.factor(continent)
  ) # Umwandlung von Kontinent zu kategorialer Variable

wdi %>%
  distinct(continent) %>%
  arrange(continent) # alle Kontinente im Datensatz


wdi <- wdi %>% arrange(., continent)

############## Mittelwertsmodell################################################

lm_mw <- lm(formula = CO2emission ~ -1 + continent, data = wdi)
summary(lm_mw)

# Interpretation der Koeffizienten
coefs_mw <- round(coefficients((lm_mw)), 3)
coefs_mw

wdi %>%
  group_by(continent) %>%
  summarize(mean(x = CO2emission, na.rm = TRUE))

# Design Matrix
mm_mw <- model.matrix(lm_mw)
cbind.data.frame(mm_mw, continent = wdi$continent) %>%
  slice_sample(n = 10) %>%
  arrange(continent)

############## Effektkodierung #################################################

lm_effect <- lm(
  formula = CO2emission ~ continent, data = wdi,
  contrasts = list(continent = "contr.sum")
)
summary(lm_effect)

# Interpretation der Koeffizienten
coefs_effect <- round(coefficients((lm_effect)), 3)
coefs_effect

group_means <- wdi %>%
  group_by(continent) %>%
  summarize(mean(x = CO2emission, na.rm = TRUE)) %>%
  tibble::deframe()
group_means
mean(group_means)

group_means - mean(group_means)
sum(group_means - mean(group_means)) # Restriktion: alle Parameter summieren zu 0 auf

# Design Matrix
mm_effekt <- model.matrix(lm_effect)
cbind.data.frame(mm_effekt, continent = wdi$continent) %>%
  slice_sample(n = 10) %>%
  arrange(continent)

# Referenz: Europe
wdi <- wdi %>% mutate(continent2 = fct_relevel(continent, "Oceania")) # Umordnen, sodass Oceania "erste" Kategorie ist
levels(wdi$continent)
levels(wdi$continent2)
lm_effect2 <- lm(formula = CO2emission ~ continent2, data = wdi, contrasts = list(continent2 = "contr.sum"))
summary(lm_effect2)


############## Referenzkodierung ###############################################

lm_ref <- lm(formula = CO2emission ~ continent, data = wdi)
summary(lm_ref)

# Interpretation der Koeffizienten
coefs_ref <- round(coef(lm_ref), 3)
coefs_ref
coefs_mw

mm_ref <- model.matrix(lm_ref)
cbind.data.frame(mm_ref, continent = wdi$continent) %>%
  slice_sample(n = 10) %>%
  arrange(continent)


# Schöne Plots
visreg(lm_mw) # gleich für lm_mw, lm_effect, lm_ref
plot(allEffects(lm_mw)) # gleich für lm_mw, lm_effect, lm_ref


# Konfidenzintervalle für Schätzer bzw. Gruppenmittelwert
confint(lm_ref) # Konfidenzintervalle für Schätzer

predict(lm_mw,
  newdata = data.frame(continent = levels(wdi$continent)),
  interval = "confidence"
) # Konfidenzintervalle für Gruppenmittelwerte (gleich für lm_mw, lm_effect, lm_ref)
