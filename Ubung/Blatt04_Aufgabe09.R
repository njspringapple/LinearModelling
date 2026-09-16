library(tidyverse)

weight <- readRDS("../Daten/weightGainRats.Rds")


# ein paar Datenanalysen
head(weight)
str(weight)
summary(weight)
table(weight$amount, weight$source)

# Mittelwert für alle Gruppen
weight %>%
  group_by(amount, source) %>%
  summarize(mean(x = weightgain, na.rm = TRUE), .groups = "drop")

# Standardabweichung für alle Gruppen
weight %>%
  group_by(amount, source) %>%
  summarize(sd(x = weightgain, na.rm = TRUE), .groups = "drop")


######### Modellierung #########################################################

model_aov <- aov(formula = weightgain ~ source * amount, data = weight)
summary(model_aov)

# p-Werte der Referenzverteilung
1 - pf(q = c(0.988, 5.812, 3.952), df1 = 1, df2 = 36)

# Nicht durch die beiden Variablen erklärte Reststreuung in der Zielgröße
SSE0 <- sum((weight$weightgain - mean(weight$weightgain))^2)
SSE0 - 221 - 1300 - 884


# Design Matrix
cbind.data.frame(model.matrix(model_aov), weight[, 1:2]) %>%
  slice_sample(n = 10)


# Koeffizienten
coef(model_aov)

# Erwartungswerte der einzelnen Gruppen
fit_table <- weight %>%
  group_by(source, amount) %>%
  summarize(y_hat = mean(weightgain))
fit_table

ggplot(
  data = fit_table,
  mapping = aes(x = amount, y = y_hat, group = source, linetype = source)
) +
  geom_line() +
  ylab("E[weightgain|source, amount]") +
  theme_minimal()


# Schätzung des varianzanalytischen Modells mit lm()
model_lm <- lm(formula = weightgain ~ source * amount, data = weight)
summary(model_lm)
