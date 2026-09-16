library(qqplotr)

############################# (a)###############################################

# Definition von Modellparametern:
beta0 <- 7
beta1 <- 0.01
beta2 <- 1
# Funktion zur Simulation von Datensätzen mit n Beobachtungen anhand des
# gegebenen Regressionsmodells:
simulation <- function(n) {
  eps <- rnorm(n = n, mean = 0, sd = sqrt(1 / 3))
  x1 <- seq(1:n)
  x2 <- x1 %% 2
  y <- beta0 + beta1 * x1 + beta2 * x2 + eps
  data <- data.frame(eps = eps, x1 = x1, x2 = x2, y = y)
  return(data)
}

# Durchführung der Simulation
set.seed(1)
sim10 <- simulation(10)
head(sim10)

sim100 <- simulation(100)
head(sim100)

sim1000 <- simulation(1000)
head(sim1000)

lm10 <- lm(formula = y ~ x1 + x2, data = sim10)
summary(lm10)

lm100 <- lm(formula = y ~ x1 + x2, data = sim100)
summary(lm100)

lm1000 <- lm(formula = y ~ x1 + x2, data = sim1000)
summary(lm1000)

coef(lm10)
coef(lm100)
coef(lm1000)

confint(lm10)
confint(lm100)
confint(lm1000)

abs((summary(lm10)$sigma)^2 - 1 / 3)
abs((summary(lm100)$sigma)^2 - 1 / 3)
abs((summary(lm1000)$sigma)^2 - 1 / 3)


############################# (b)###############################################
par(mfrow = c(2, 2))
plot(x = lm1000, col = rgb(0, 0, 0, alpha = 0.3))

# gefittete Werte und Residuen zu Datensätzen hinzufügen
sim10 <- sim10 %>% mutate(fitted = fitted(lm10), residuals = residuals(lm10))
sim100 <- sim100 %>%
  mutate(fitted = fitted(lm100), residuals = residuals(lm100))
sim1000 <- sim1000 %>%
  mutate(fitted = fitted(lm1000), residuals = residuals(lm1000))


# gefittete vs. beobachtete Werte (Linearität und Homoskedastizität)
plot10 <- ggplot(data = sim10, mapping = aes(x = fitted, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  ggtitle("n = 10")
plot100 <- ggplot(data = sim100, mapping = aes(x = fitted, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  ggtitle("n = 100")
plot1000 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = y)) +
  geom_point(alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  ggtitle("n = 1000")
ggarrange(plot10, plot100, plot1000, ncol = 3)

# gefittete Werte vs. Residuen (Linearität und Homoskedastizität)
plot10 <- ggplot(data = sim10, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  ggtitle("n = 10")
plot100 <- ggplot(data = sim100, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  ggtitle("n = 100")
plot1000 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  ggtitle("n = 1000")
ggarrange(plot10, plot100, plot1000, ncol = 3)

# Normal-Quantil-Plot (Normaverteiltheit)
plot10 <- ggplot(data = sim10, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line() +
  theme_minimal() +
  ggtitle("n = 10")
plot100 <- ggplot(data = sim100, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line() +
  theme_minimal() +
  ggtitle("n = 100")
plot1000 <- ggplot(data = sim1000, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line() +
  theme_minimal() +
  ggtitle("n = 1000")
ggarrange(plot10, plot100, plot1000, ncol = 3)


############################# (c)###############################################

# ohne x1
lm10_ohnex1 <- lm(formula = y ~ x2, data = sim10)
lm100_ohnex1 <- lm(formula = y ~ x2, data = sim100)
lm1000_ohnex1 <- lm(formula = y ~ x2, data = sim1000)


sim10 <- sim10 %>% mutate(
  fitted = fitted(lm10_ohnex1),
  residuals = residuals(lm10_ohnex1)
)
sim100 <- sim100 %>%
  mutate(fitted = fitted(lm100_ohnex1), residuals = residuals(lm100_ohnex1))
sim1000 <- sim1000 %>%
  mutate(fitted = fitted(lm1000_ohnex1), residuals = residuals(lm1000_ohnex1))

# gefittete vs. beobachtete Werte (Linearität und Homoskedastizität)
plot10 <- ggplot(data = sim10, mapping = aes(x = fitted, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  ggtitle("n = 10")
plot100 <- ggplot(data = sim100, mapping = aes(x = fitted, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  ggtitle("n = 100")
plot1000 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = y)) +
  geom_point(alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  ggtitle("n = 1000")
ggarrange(plot10, plot100, plot1000, ncol = 3)

# gefittete Werte vs. Residuen (Linearität und Homoskedastizität)
plot10 <- ggplot(data = sim10, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  ggtitle("n = 10")
plot100 <- ggplot(data = sim100, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  ggtitle("n = 100")
plot1000 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  ggtitle("n = 1000")
ggarrange(plot10, plot100, plot1000, ncol = 3)

# Normal-Quantil-Plot (Normaverteiltheit)
plot10 <- ggplot(data = sim10, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line() +
  theme_minimal() +
  ggtitle("n = 10")
plot100 <- ggplot(data = sim100, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line() +
  theme_minimal() +
  ggtitle("n = 100")
plot1000 <- ggplot(data = sim1000, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line() +
  theme_minimal() +
  ggtitle("n = 1000")
ggarrange(plot10, plot100, plot1000, ncol = 3)

summary(lm1000)
summary(lm1000_ohnex1)


# ohne x2
lm10_ohnex2 <- lm(formula = y ~ x1, data = sim10)
lm100_ohnex2 <- lm(formula = y ~ x1, data = sim100)
lm1000_ohnex2 <- lm(formula = y ~ x1, data = sim1000)

sim10 <- sim10 %>% mutate(
  fitted = fitted(lm10_ohnex2),
  residuals = residuals(lm10_ohnex2)
)
sim100 <- sim100 %>%
  mutate(fitted = fitted(lm100_ohnex2), residuals = residuals(lm100_ohnex2))
sim1000 <- sim1000 %>%
  mutate(fitted = fitted(lm1000_ohnex2), residuals = residuals(lm1000_ohnex2))

# gefittete vs. beobachtete Werte (Linearität und Homoskedastizität)
plot10 <- ggplot(data = sim10, mapping = aes(x = fitted, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  ggtitle("n = 10")
plot100 <- ggplot(data = sim100, mapping = aes(x = fitted, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  ggtitle("n = 100")
plot1000 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = y)) +
  geom_point(alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1) +
  theme_minimal() +
  ggtitle("n = 1000")
ggarrange(plot10, plot100, plot1000, ncol = 3)

# gefittete Werte vs. Residuen (Linearität und Homoskedastizität)
plot10 <- ggplot(data = sim10, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  ggtitle("n = 10")
plot100 <- ggplot(data = sim100, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  ggtitle("n = 100")
plot1000 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  ggtitle("n = 1000")
ggarrange(plot10, plot100, plot1000, ncol = 3)

# Normal-Quantil-Plot (Normaverteiltheit)
plot10 <- ggplot(data = sim10, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line() +
  theme_minimal() +
  ggtitle("n = 10")
plot100 <- ggplot(data = sim100, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line() +
  theme_minimal() +
  ggtitle("n = 100")
plot1000 <- ggplot(data = sim1000, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line() +
  theme_minimal() +
  ggtitle("n = 1000")
ggarrange(plot10, plot100, plot1000, ncol = 3)

summary(lm1000)
summary(lm1000_ohnex2)



############################# (d)###############################################

# (i) Outlier
yneu <- sim1000$y
yneu[10] <- 20 # 10. Beobachtung wird y-Wert von 20

lm1000i <- lm(formula = yneu ~ x1 + x2, data = sim1000)
summary(lm1000i)

# Generierung des benötigten Datensatzes:
sim1000 <- sim1000 %>%
  mutate(fitted = fitted(lm1000i), residuals = residuals(lm1000i))
# Gefittete Wertes vs. beobachtete Werte:
plot1 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = yneu)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  theme_minimal() +
  ggtitle("fitted vs. observed")
# Gefittete Werte vs. Residuen:
plot2 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, col = "red") +
  theme_minimal() +
  ggtitle("fitted vs residuals")
# Normal-Quantil-Plot:
plot3 <- ggplot(data = sim1000, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line(col = "red") +
  theme_minimal() +
  ggtitle("normal quantile plot")
# Gemeinsame Visualisierung:
ggarrange(plot1, plot2, plot3, ncol = 3)

# Cook's distance
head(sort(x = cooks.distance(lm1000i), decreasing = TRUE))
par(mfrow = c(2, 2))
plot(x = lm1000i, col = rgb(0, 0, 0, alpha = 0.3))


# (ii) gleichverteilte Fehler epsilon
# Funktion zur Simulation des neuen Datensatzes mit gleichverteilten
# Fehlern:
simulation_unif <- function(n) {
  eps <- runif(n = n, min = -1, max = 1)
  x1 <- seq(1:n)
  x2 <- x1 %% 2
  y <- beta0 + beta1 * x1 + beta2 * x2 + eps
  data <- data.frame(eps = eps, x1 = x1, x2 = x2, y = y)
  return(data)
}
# Durchführung der Simulationen:
set.seed(1)
sim1000 <- simulation_unif(1000)
head(sim1000)

# Modell
lm1000ii <- lm(formula = y ~ x1 + x2, data = sim1000)
summary(lm1000ii)

# Generierung des benötigten Datensatzes:
sim1000 <- sim1000 %>%
  mutate(fitted = fitted(lm1000ii), residuals = residuals(lm1000ii))
# Gefittete Wertes vs. beobachtete Werte:
plot1 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  theme_minimal() +
  ggtitle("fitted vs. observed")
# Gefittete Werte vs. Residuen:
plot2 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, col = "red") +
  theme_minimal() +
  ggtitle("fitted vs residuals")
# Normal-Quantil-Plot:
plot3 <- ggplot(data = sim1000, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line(col = "red") +
  theme_minimal() +
  ggtitle("normal quantile plot")
# Gemeinsame Visualisierung:
ggarrange(plot1, plot2, plot3, ncol = 3)


# (iii) heteroskedastische Fehler
# Funktion zur Simulation des neuen Datensatzes mit heteoskedastischen
# Fehlern:
simulation_het <- function(n) {
  x1 <- seq(1:n)
  eps <- rnorm(n = n, mean = 0, sd = 2 * x1 / 300) # Fehler hängen von x1 ab
  x2 <- x1 %% 2
  y <- beta0 + beta1 * x1 + beta2 * x2 + eps
  data <- data.frame(eps = eps, x1 = x1, x2 = x2, y = y)
  return(data)
}
# Durchführung der Simulationen:
set.seed(1)
sim1000 <- simulation_het(1000)
head(sim1000)

lm1000iii <- lm(formula = y ~ x1 + x2, data = sim1000)
summary(lm1000iii)

# Generierung des benötigten Datensatzes:
sim1000 <- sim1000 %>%
  mutate(fitted = fitted(lm1000iii), residuals = residuals(lm1000iii))
# Gefittete Wertes vs. beobachtete Werte:
plot1 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  theme_minimal() +
  ggtitle("fitted vs. observed")
# Gefittete Werte vs. Residuen:
plot2 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, col = "red") +
  theme_minimal() +
  ggtitle("fitted vs residuals")
# Normal-Quantil-Plot:
plot3 <- ggplot(data = sim1000, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line(col = "red") +
  theme_minimal() +
  ggtitle("normal quantile plot")
# Gemeinsame Visualisierung:
ggarrange(plot1, plot2, plot3, ncol = 3)



# (iv) nichtlinerarer Zusammenhang
# Funktion zur Simulation des neuen Datensatzes:
simulation_nonlinear <- function(n) {
  eps <- rnorm(n = n, mean = 0, sd = sqrt(1 / 3))
  x1 <- seq(1:n)
  x2 <- x1 %% 2
  y <- beta0 + beta1 * (x1^2 / 100) + beta2 * x2 + eps
  data <- data.frame(eps = eps, x1 = x1, x2 = x2, y = y)
  return(data)
}
# Durchführung der Simulationen:
set.seed(1)
sim1000 <- simulation_nonlinear(1000)
head(sim1000)

lm1000iv <- lm(formula = y ~ x1 + x2, data = sim1000)
summary(lm1000iv)


# Generierung des benötigten Datensatzes:
sim1000 <- sim1000 %>%
  mutate(fitted = fitted(lm1000iv), residuals = residuals(lm1000iv))
# Gefittete Wertes vs. beobachtete Werte:
plot1 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  theme_minimal() +
  ggtitle("fitted vs. observed")
# Gefittete Werte vs. Residuen:
plot2 <- ggplot(data = sim1000, mapping = aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, col = "red") +
  theme_minimal() +
  ggtitle("fitted vs residuals")
# Normal-Quantil-Plot:
plot3 <- ggplot(data = sim1000, mapping = aes(sample = residuals)) +
  stat_qq_point() +
  stat_qq_line(col = "red") +
  theme_minimal() +
  ggtitle("normal quantile plot")
# Gemeinsame Visualisierung:
ggarrange(plot1, plot2, plot3, ncol = 3)
