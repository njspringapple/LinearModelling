library(tidyverse)
library(qqplotr)

############################# (a) ##############################################
load("../Daten/wdi.Rdata")
head(wdi)
str(wdi)
summary(wdi)

wdi <- wdi %>%
  dplyr::select(
    country, year, continent, region, continent, year, CO2emission, Area,
    Population, PopulationRural, GDP, Livestock, Employees.M.IND
  ) %>%
  mutate(
    CO2emission = log10(CO2emission), Area = log10(Area),
    Population = log10(Population),
    PopulationRural = log10(PopulationRural + 1), GDP = log10(GDP),
    continent = as.factor(continent),
    year = as.numeric(as.character(year))
  ) %>%
  filter(year %in% 2005:2013) %>%
  na.omit()

wdi_2005 <- wdi %>%
  filter(year == 2005) %>%
  arrange(desc(CO2emission)) %>%
  slice(1:9)
countries <- as.character(wdi_2005$country)
countries

ggplot(
  data = filter(wdi, country %in% countries),
  mapping = aes(x = year, y = CO2emission, col = country)
) +
  geom_point() +
  geom_line() +
  theme_minimal()


############################# (b) ##############################################
model_lm <- lm(
  formula = CO2emission ~ Area + PopulationRural + Population +
    GDP + Livestock + Employees.M.IND + continent,
  data = wdi
)
summary(model_lm)


############################# (e) ##############################################
model_mixed <- lmer(
  formula = CO2emission ~ Area + PopulationRural +
    Population + GDP + Livestock + Employees.M.IND +
    continent + (1 | country),
  data = wdi
)
summary(model_mixed)

random <- ranef(model_mixed)$country
random$country <- rownames(random)
colnames(random)[1] <- "Intercept"
sd_re <- as.data.frame(VarCorr(model_mixed))[1, "sdcor"]

ggplot(data = random, mapping = aes(sample = Intercept)) +
  stat_qq_point() +
  stat_qq_line() +
  theme_minimal()

ggplot(random, aes(x = Intercept)) +
  geom_density(aes(color = "Dichteschätzung"), linewidth = 1) +
  stat_function(aes(color = "Normalverteilung"),
    fun = dnorm,
    args = list(mean = 0, sd = sd_re),
    linewidth = 1
  ) +
  labs(
    x = "Random Intercept",
    y = "Dichte"
  ) +
  scale_color_manual(
    values = c(
      "Dichteschätzung" = "black",
      "Normalverteilung" = "red"
    ),
    name = NULL
  ) +
  theme(
    legend.position = "top"
  ) +
  theme_minimal()
