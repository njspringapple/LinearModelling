library(tidyverse)
library(ggpubr)

wagen <- read.csv("../Daten/wagen.txt", header = TRUE, sep = "")
head(wagen)
str(wagen)
summary(wagen)


############################# (a)################################################
# Alter
plot1 <- ggplot(data = wagen, mapping = aes(x = alter, y = preis)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  ggtitle("Alter") +
  xlab("Alter") +
  ylab("Verkaufspreis in 1000 Euro") +
  theme_minimal()
# Kilometerstand:
plot2 <- ggplot(data = wagen, mapping = aes(x = kilstand, y = preis)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  ggtitle("Kilometerstand") +
  xlab("Kilometerstand") +
  ylab("Verkaufspreis in 1000 Euro") +
  theme_minimal()
ggarrange(plot1, plot2)

############################# (b)################################################
model_full <- lm(
  formula = preis ~ alter + I(alter^2) + I(alter^3) +
    I(alter^4) + kilstand + I(kilstand^2) + I(kilstand^3) +
    I(kilstand^4),
  data = wagen
)
summary(model_full)

# Grafische Visualisierung
# Alter:
effect_alter <- wagen %>% mutate(kilstand = mean(kilstand))
effect_alter <- effect_alter %>%
  mutate(prediction_full = predict(object = model_full, newdata = effect_alter))
plot_alter <- ggplot(data = effect_alter, mapping = aes(x = alter, y = preis)) +
  geom_point() +
  geom_line(aes(y = prediction_full), col = "blue") +
  theme_minimal() +
  ggtitle("Alter") +
  xlab("Alter") +
  ylab("Verkaufspreis in 1000 Euro")
# Kilometerstand:
effect_km <- wagen %>% mutate(alter = mean(alter))
effect_km <- effect_km %>%
  mutate(prediction_full = predict(object = model_full, newdata = effect_km))
plot_km <- ggplot(data = effect_km, mapping = aes(x = kilstand, y = preis)) +
  geom_point() +
  geom_line(aes(y = prediction_full), col = "blue") +
  theme_minimal() +
  ggtitle("Kilometerstand") +
  xlab("Kilometerstand") +
  ylab("Verkaufspreis in 1000 Euro")
ggarrange(plot_alter, plot_km)

anova(model_full)

SSE <- anova(model_full)$`Sum Sq`[9] # Residual Sum of Squares
SSE
SSM <- sum(anova(model_full)$`Sum Sq`[1:8]) # erklärte Streuung durch Terme im Modell
SSM
SST <- SSM + SSE # Gesamtstreuung
SST
logLik(model_full)

# R^2
1 - SSE / SST
summary(model_full)$r.squared

# adjusted R^2
1 - summary(model_full)$sigma^2 / (SST / (nrow(wagen) - 1))
summary(model_full)$adj.r.squared

# AIC
-2 * logLik(model_full) + 2 * 10
AIC(model_full)

############################# (c)################################################

# definiere Modelle für Polynomgrad 1 - 4 je Variable
model_11 <- lm(formula = preis ~ alter + kilstand, data = wagen)
model_21 <- lm(
  formula = preis ~ alter + I(alter^2) + kilstand,
  data = wagen
)
model_31 <- lm(
  formula = preis ~ alter + I(alter^2) + I(alter^3) + kilstand,
  data = wagen
)
model_41 <- lm(
  formula = preis ~ alter + I(alter^2) + I(alter^3) + I(alter^4) +
    kilstand,
  data = wagen
)
model_12 <- lm(
  formula = preis ~ alter + kilstand + I(kilstand^2),
  data = wagen
)
model_22 <- lm(
  formula = preis ~ alter + I(alter^2) + kilstand + I(kilstand^2),
  data = wagen
)
model_32 <- lm(
  formula = preis ~ alter + I(alter^2) + I(alter^3) + kilstand +
    I(kilstand^2),
  data = wagen
)
model_42 <- lm(
  formula = preis ~ alter + I(alter^2) + I(alter^3) + I(alter^4) +
    kilstand + I(kilstand^2),
  data = wagen
)
model_13 <- lm(
  formula = preis ~ alter + kilstand + I(kilstand^2) +
    I(kilstand^3),
  data = wagen
)
model_23 <- lm(
  formula = preis ~ alter + I(alter^2) + kilstand +
    I(kilstand^2) + I(kilstand^3),
  data = wagen
)
model_33 <- lm(
  formula = preis ~ alter + I(alter^2) + I(alter^3) + kilstand +
    I(kilstand^2) + I(kilstand^3),
  data = wagen
)
model_43 <- lm(
  formula = preis ~ alter + I(alter^2) + I(alter^3) +
    I(alter^4) + kilstand + I(kilstand^2) + I(kilstand^3),
  data = wagen
)
model_14 <- lm(
  formula = preis ~ alter + kilstand + I(kilstand^2) +
    I(kilstand^3) + I(kilstand^4),
  data = wagen
)
model_24 <- lm(
  formula = preis ~ alter + I(alter^2) + kilstand +
    I(kilstand^2) + I(kilstand^3) + I(kilstand^4),
  data = wagen
)
model_34 <- lm(
  formula = preis ~ alter + I(alter^2) + I(alter^3) + kilstand +
    I(kilstand^2) + I(kilstand^3) + I(kilstand^4),
  data = wagen
)
model_44 <- lm(
  formula = preis ~ alter + I(alter^2) + I(alter^3) + I(alter^4) +
    kilstand + I(kilstand^2) + I(kilstand^3) + I(kilstand^4),
  data = wagen
)


# Berechnung von R^2, adj. R^2 und AIC je Modell
selection <- data.frame(matrix(data = 0, nrow = 16, ncol = 5))
colnames(selection) <- c(
  "Grad Alter", "Grad Kilometerstand", "R^2", "R_adj^2",
  "AIC"
)
selection[1, ] <- c(
  1, 1, summary(model_11)$r.squared,
  summary(model_11)$adj.r.squared, AIC(model_11)
)
selection[2, ] <- c(
  2, 1, summary(model_21)$r.squared,
  summary(model_21)$adj.r.squared, AIC(model_21)
)
selection[3, ] <- c(
  3, 1, summary(model_31)$r.squared,
  summary(model_31)$adj.r.squared, AIC(model_31)
)
selection[4, ] <- c(
  4, 1, summary(model_41)$r.squared,
  summary(model_41)$adj.r.squared, AIC(model_41)
)
selection[5, ] <- c(
  1, 2, summary(model_12)$r.squared,
  summary(model_12)$adj.r.squared, AIC(model_12)
)
selection[6, ] <- c(
  2, 2, summary(model_22)$r.squared,
  summary(model_22)$adj.r.squared, AIC(model_22)
)
selection[7, ] <- c(
  3, 2, summary(model_32)$r.squared,
  summary(model_32)$adj.r.squared, AIC(model_32)
)
selection[8, ] <- c(
  4, 2, summary(model_42)$r.squared,
  summary(model_42)$adj.r.squared, AIC(model_42)
)
selection[9, ] <- c(
  1, 3, summary(model_13)$r.squared,
  summary(model_13)$adj.r.squared, AIC(model_13)
)
selection[10, ] <- c(
  2, 3, summary(model_23)$r.squared,
  summary(model_23)$adj.r.squared, AIC(model_23)
)
selection[11, ] <- c(
  3, 3, summary(model_33)$r.squared,
  summary(model_33)$adj.r.squared, AIC(model_33)
)
selection[12, ] <- c(
  4, 3, summary(model_43)$r.squared,
  summary(model_43)$adj.r.squared, AIC(model_43)
)
selection[13, ] <- c(
  1, 4, summary(model_14)$r.squared,
  summary(model_14)$adj.r.squared, AIC(model_14)
)
selection[14, ] <- c(
  2, 4, summary(model_24)$r.squared,
  summary(model_24)$adj.r.squared, AIC(model_24)
)
selection[15, ] <- c(
  3, 4, summary(model_34)$r.squared,
  summary(model_34)$adj.r.squared, AIC(model_34)
)
selection[16, ] <- c(
  4, 4, summary(model_44)$r.squared,
  summary(model_44)$adj.r.squared, AIC(model_44)
)
selection

selection[which(selection$`R^2`
== max(selection$`R^2`)), ]
selection[which(selection$`R_adj^2`
== max(selection$`R_adj^2`)), ]
selection[which(selection$AIC == min(selection$AIC)), ]

# Visualisierung
# Alter
effect_alter <- effect_alter %>%
  mutate(prediction_best = predict(
    object = model_22,
    newdata = effect_alter
  )) %>%
  gather(prediction_full, prediction_best, key = "model", value = "prediction")
plot_alter <- ggplot(data = effect_alter, mapping = aes(x = alter, y = preis)) +
  geom_point() +
  geom_line(aes(y = prediction, col = model)) +
  theme_minimal() +
  ggtitle("Alter") +
  xlab("Alter") +
  ylab("Verkaufspreis in 1000 Euro") +
  scale_color_discrete(
    name = "",
    labels = c("Bestes Modell", "Volles Modell")
  )
# Kilometerstand:
effect_km <- effect_km %>%
  mutate(prediction_best = predict(object = model_22, newdata = effect_km)) %>%
  gather(prediction_full, prediction_best, key = "model", value = "prediction")
plot_km <- ggplot(data = effect_km, mapping = aes(x = kilstand, y = preis)) +
  geom_point() +
  geom_line(aes(y = prediction, col = model)) +
  theme_minimal() +
  ggtitle("Kilometerstand") +
  xlab("Kilometerstand") +
  ylab("Verkaufspreis in 1000 Euro") +
  scale_color_discrete(
    name = "",
    labels = c("Bestes Modell", "Volles Modell")
  )
ggarrange(
  plotlist = list(plot_alter, plot_km), legend = "bottom",
  common.legend = TRUE
)

############################# (d)################################################
model_spline <- lm(
  formula = preis ~ alter + I(alter^2) +
    I((alter > 89.5) * (alter - 90)^2) +
    I((alter > 109.5) * (alter - 110)^2) +
    I((alter > 129.5) * (alter - 130)^2) +
    kilstand + I(kilstand^2) +
    I((kilstand > 99.5) * (kilstand - 100)^2) +
    I((kilstand > 149.5) * (kilstand - 150)^2) +
    I((kilstand > 199.5) * (kilstand - 200)^2),
  data = wagen
)
summary(model_spline)

# Adjustiertes Bestimmtheitsmaß:
summary(model_spline)$adj.r.squared
summary(model_22)$adj.r.squared

# AIC:
AIC(model_spline)
AIC(model_22)

# Alter:
effect_alter <- wagen %>% mutate(kilstand = mean(kilstand))
effect_alter <- effect_alter %>%
  mutate(
    prediction_full = predict(object = model_full, newdata = effect_alter),
    prediction_best = predict(object = model_22, newdata = effect_alter),
    prediction_spline = predict(
      object = model_spline,
      newdata = effect_alter
    )
  ) %>%
  gather(prediction_full, prediction_best, prediction_spline,
    key = "model", value = "prediction"
  )
plot_alter <- ggplot(data = effect_alter, mapping = aes(x = alter, y = preis)) +
  geom_point() +
  geom_line(aes(y = prediction, col = model)) +
  theme_minimal() +
  ggtitle("Alter") +
  xlab("Alter") +
  ylab("Verkaufspreis in 1000 Euro") +
  scale_color_discrete(
    name = "",
    labels = c(
      "Bestes Modell", "Volles Modell",
      "Spline-Modell"
    )
  )
# Kilometerstand:
effect_km <- wagen %>% mutate(alter = mean(alter))
effect_km <- effect_km %>%
  mutate(
    prediction_full = predict(object = model_full, newdata = effect_km),
    prediction_best = predict(object = model_22, newdata = effect_km),
    prediction_spline = predict(
      object = model_spline,
      newdata = effect_km
    )
  ) %>%
  gather(prediction_full, prediction_best, prediction_spline,
    key = "model", value = "prediction"
  )
plot_km <- ggplot(data = effect_km, mapping = aes(x = kilstand, y = preis)) +
  geom_point() +
  geom_line(aes(y = prediction, col = model)) +
  theme_minimal() +
  ggtitle("Kilometerstand") +
  xlab("Kilometerstand") +
  ylab("Verkaufspreis in 1000 Euro") +
  scale_color_discrete(
    name = "",
    labels = c(
      "Bestes Modell", "Volles Modell",
      "Spline-Modell"
    )
  )
# Gemeinsame Visualisierung:
ggarrange(
  plotlist = list(plot_alter, plot_km), legend = "bottom",
  common.legend = TRUE
)
