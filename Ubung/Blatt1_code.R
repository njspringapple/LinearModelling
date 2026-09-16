# Blatt 1

# Laden von benötigten Paketen:
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)

# Definition eines ggplot-Themes:
theme <- theme_classic() +
  theme(text = element_text(size = 12), axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12, hjust = 0),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.text.y = element_text(size = 12),
        strip.placement = "outside", strip.background = element_blank(),
        axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
        axis.title.x = element_text(margin = margin(10, 0, 0, 0)))



# Aufgabe 1
# a)
run <- readRDS("Uebung26/Daten/RunningAgg.Rds")
head(run)
str(run)
summary(run)

# Verteilung der Variable
hist(run$pace)
hist(run$HR)

# Streudiagramm
gg_pace <- ggplot(data = run, mapping = aes(x = pace, y = HR)) + geom_point() +
  xlab("Pace (min/km)") + ylab("HR (bpm)") + theme
gg_pace

lm_pace <- lm(formula = HR ~ pace, data = run)
summary(lm_pace)
beta_pace <- round(x = coef(lm_pace), digits = 2)
beta_pace

gg_pace + geom_smooth(method = "lm", se = FALSE)



# b) 
run$speed <- 60 / run$pace
lm_speed <- lm(formula = HR ~ speed, data = run)
summary(lm_speed)
beta_speed <- round(x = coef(lm_speed), digits = 2)
beta_speed

# Pace und Herzfrequenz:
gg_original <- gg_pace + geom_smooth(method = "lm", se = FALSE) +
  ggtitle("Pace") + theme
# Geschwindigkeit und Herzfrequenz:
gg_transformed <- gg_original + run + aes(x=speed) +
  ggtitle("Geschwindigkeit") + xlab("Geschwindigkeit (km/h)") + theme
# Gemeinsame Visualisierung:
grid.arrange(gg_original, gg_transformed, nrow = 1)


# d)
# Berechnung der transformierten Variable:
run <- run %>% mutate(HRbps = HR / 60, speedMi = speed / 1.61)

# Schätzung eines linearen Regressionsmodells zwischen der Geschwindidkeit und
# der Herzfrequenz auf transformierten Skalen:
lm_trafo <- lm(formula = HRbps ~ speedMi, data = run)
summary(lm_trafo)

# Vergleich des Modelloutputs mit Teilaufgabe (b):
a1 <- 1 / 1.61
b1 <- 1 / 60
(b1 / a1) * coefficients(lm_speed)[2]
b1 * coefficients(lm_speed)[1]
coefficients(lm_trafo)


# e)
mean(run$speed)

run <- run %>% mutate(c_speed = scale(x = speed, center = TRUE, scale = FALSE))

# test if all equal run$c_speed and run$speed - mean(run$speed)
all.equal(c(run$c_speed), run$speed - mean(run$speed), check.attributes = FALSE)

run %>% gather(variable, value, speed, c_speed) %>%
  ggplot(aes(x = variable, y = value)) +
  geom_boxplot() + geom_hline(yintercept=0, col=2, lty=3) +
  xlab("Variable") + ylab("Geschwindigkeit (km/h)") + theme


lm_c_speed <- lm(formula = HR ~ c_speed, data = run)
summary(lm_c_speed)

# g)

run <- run %>% mutate(s_speed = scale(speed), s_HR = scale(HR))
summary(select(run, s_speed, s_HR))

lm_s_speed <- lm(formula = s_HR ~ s_speed, data = run)
summary(lm_s_speed)

cor(x = run$HR, y = run$speed, method = "pearson")
cor(x = run$s_HR, y = run$s_speed, method = "pearson")
coef(lm_s_speed)["s_speed"]
all.equal(cor(run$HR, run$speed), coef(lm_s_speed)["s_speed"],
          check.attributes = FALSE)




# Aufgabe 2
# b)
# Setzen der Zufallszahlen:
set.seed(3456)
# Definition der Modellparameter:
N <- 10000
n <- 100
beta.0 <- -2
beta.1 <- 3.5
sigma.sq <- 100

# Zufällige Ziehung von N Beobachtungen der Einflussgröße aus einer
# Gleichverteilung:
x_unif <-  runif(n = N, min = 0, max = n)
# Alternative: z.B. Ziehung aus einer Exponentialverteilung
# Zufällige Ziehung des Fehlerterms epsilon aus einer Normalverteilung:
epsilon <- rnorm(N, mean = 0, sd = sqrt(sigma.sq))
# Berechnung der Werte der Zielgröße über Regressionsmodell:
y_vals <- beta.0 + beta.1 * x_unif + epsilon
# Abspeichern der Informationen in einem Datensatz:
predictor <- x_unif
response <- y_vals
data_sim <- data.frame(predictor, response)

var_b0_true <- (sigma.sq / n) * (1 + (mean(predictor)^2) / var(predictor))
var_b0_true
var_b1_true <- sigma.sq / (n * var(predictor))
var_b1_true

reps <- 10000
# Matrix der Ergebnisse:
fit <- matrix(ncol = 2, nrow = reps)
# for-Schleife über die Wiederholungen:
for (i in 1:reps){
  sample <- data_sim[sample(1:N, n), ]
  fit[i, ] <- lm(response ~ predictor, data = sample)$coefficients
}

# Erhaltene Varianzschätzungen:
var(fit[, 1])
var(fit[, 2])

par(mfrow = c(1, 2))
# Achsenabschnitt:
hist(x = fit[, 1], cex.main = 1,
     main = bquote(Distribution  ~ of ~ 10000 ~ beta[0] ~ estimates), 
     xlab = bquote(hat(beta)[0]), freq = FALSE)
curve(dnorm(x = x,  mean = -2, sd = sqrt(var_b0_true)), add = TRUE, 
      col = "darkred")
# Steigungsparameter:
hist(x = fit[, 2], cex.main = 1,
     main = bquote(Distribution  ~ of ~ 10000 ~ beta[1] ~ estimates), 
     xlab = bquote(hat(beta)[1]), freq = FALSE)
curve(dnorm(x = x,  mean = 3.5, sd = sqrt(var_b1_true)), add = TRUE, 
      col = "darkred")






# Aufgabe 3
# a) 
data <- data.frame("Groesse" = c(198, 188, 196, 190, 180, 183, 196, 196, 193,
                                 183),
                   "Gewicht" = c(104, 84, 107, 95, 76, 79, 109, 94, 113, 93))
data
summary(data)
# Graphische Visualisierung des Zusammenhangs in einem Streudiagramm:
ggplot(data = data, mapping = aes(x = Groesse, y = Gewicht)) + 
  geom_point() +
  geom_hline(yintercept = mean(data$Gewicht), linetype= "dashed") +
  xlab("Körpergröße [in cm]") + ylab("Gewicht [in kg]") + theme


# b)
SST <- sum((data$Gewicht - mean(data$Gewicht))^2)
SST

# c) 
beta1 <- cov(data$Gewicht, data$Groesse) / var(data$Groesse)
beta1

beta0 <- mean(data$Gewicht) - beta1 * mean(data$Groesse)
beta0

# d)
pred <- beta0 + beta1 * data$Groesse
SSM <- sum((pred - mean(data$Gewicht))^2)
SSM

R2 <- SSM / SST
R2

# e) 
model <- lm(formula = Gewicht ~ Groesse, data = data)
summary(model)

ggplot(data = data, mapping = aes(x = Groesse, y = Gewicht)) + geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  xlab("Körpergröße [in cm]") + ylab("Gewicht [in kg]") + theme

model <- lm(Gewicht ~ Groesse, data=data)
summary(model)
pred <- predict(model)
SSM <- sum((pred-mean(data$Gewicht))^2)
print(SSM)

