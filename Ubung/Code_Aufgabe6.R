# Blatt 3

# Aufgabe 6

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)

# a) 
# Daten einlesen
load("Uebung26/Daten/wdi.Rdata")

# Herausfiltern der Daten aus dem Jahr 2013,
# Selektion benötigter Variablen und Umrechnung der Bevölkerung in Tausend,
# Hnzufügen von log10-transformierten Variablen zum Datensatz:
wdi <- wdi %>% select(continent, country, year, Population, Area, CO2emission) %>%
  filter(year == 2013)%>%
  filter(!(is.na(Population) | is.na(Area)))%>%
  mutate(Population = Population / 1000,
         log10.Population = log10(Population),
         log10.CO2emission = log10(CO2emission),
         log10.Area = log10(Area))


# b) 
# Aufteilung der Daten in Trainings- und Testdaten:
set.seed(29042020)
wdi_test <- sample_n(tbl = wdi, size = 100)
wdi_train <- wdi %>% filter(!(country %in% wdi_test$country))
# Schätzung eines linearen Regressionsmodells unter Verwendung der
# log10-transformierten Daten:
lm_log10 <- lm(formula = log10.Population ~ log10.Area, data = wdi_train)
coef_lm_log10 <- round(x = coefficients(lm_log10), digits = 2)
coef_lm_log10


# c) 
# Berechnung der 99% KI: 
KI_beta1 <- lm_log10$coef[2] + c(-1, 1) * summary(lm_log10)$coef[2, 2] *
  qt(p = 0.995, df = nrow(wdi_train) - 2)
KI_beta1



# e) 
# Berechnung eines 95%-Konfidenzintervalls für die erwartete logarithmierte Bevölkerungszahl
# anhand der Beobachtungen des Trainingsdatensatzes:
conf_int <- predict(object = lm_log10, interval = "confidence", level = 0.95)
head(conf_int)

# Graphische Visualisierung der Regressionsgerade mit Konfidenzintervall:
gg_conf <- ggplot(data = wdi_train,
                  mapping = aes(x = log10.Area, y = log10.Population,
                                label = country)) +
  geom_point() + geom_smooth(method = "lm", se = TRUE, level = 0.95) +
  xlab("log10(Fläche) [km^2]") +
  ylab("log10(Bevölkerung) [in Tsd.]") + theme()
gg_conf

# Anteil der Beobachtungen des Trainingsdatensatzes innerhalb des Konfidenzintervalls:
is_covered <- function(y, lwr, upr) {
  mean(lwr <= y & y <= upr)
}
in_ci <- is_covered(wdi_train$log10.Population, conf_int[, "lwr"],
                    conf_int[, "upr"])
in_ci

# Ca. 23.7% der Beobachtungen liegen innerhalb des KI 


# f) 
# Berechnung eines 95%-Prognoseintervalls für die aus der Schätzung herausgelassenen Länder:
pred_int <- predict(object = lm_log10, newdata = wdi_test,
                      interval = "prediction", level = 0.95)
head(pred_int)

# Anteil der Beobachtungen des Testdatensatzes innerhalb des Prognoseintervalls:
in_pi <- is_covered(wdi_test$log10.Population, pred_int[, "lwr"],
pred_int[, "upr"])
in_pi

# 98 der 100 bei der Schätzung herausgelassenen Länder liegen innerhalb des 95%-
# Prognoseintervalls. ⇒ insgesamt gute Coverage

# Graphischer Vergleich zwischen 95%-Konfidenzintervall und 95%-Prognoseintervall:
pred_int <- cbind(wdi_test, pred_int)
gg_prog <- gg_conf +
  geom_line(data = pred_int, mapping = aes(y = lwr), col = "red") +
  geom_line(data = pred_int, mapping = aes(y = upr), col = "red")
gg_prog

# Das Prognoseintervall ist deutlich weiter als das Konfidenzintervall, da die Berechnung
# zusätzlich die Streuung von ϵ_n+1 berücksichtigt.
