# Blatt 3

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(car)

# Aufgabe 7

# a) 
# Einlesen des Datensatzes und Überblick über Daten:
sambia <- read.table(file = "Uebung26/Daten/sambia92.raw", header = TRUE)
head(sambia)
str(sambia)
summary(sambia)

# Beschränkung auf Beobachtungen aus der am häufigsten vorkommenden Region Sambias:
table(sambia$region)
sambia <- sambia %>% filter(region == 2) %>% select(zscore:m_bmi)
head(sambia)

# Für die weitere Aufgabe wird nur noch Region 2 betrachtet.

# Graphische Visualisierung von nach Geschlecht getrennten metrischen Variablen:
# Boxplots:
  ggplot(data = gather(sambia, variable, value, -k_geschl)) +
  geom_boxplot(aes(y = value, x = factor(k_geschl))) +
  facet_wrap(~variable, scale = "free") + theme()

# Dichteschätzungen:
ggplot(data = gather(sambia, variable, value, -k_geschl)) +
  geom_density(aes(x = value, col = factor(k_geschl))) +
  facet_wrap(~variable, scale = "free") + theme()


# Anhand der beiden Visualisierungen sind keine groÿen Unterschiede zwischen den Geschlechtern erkennbar.
# In den Dichteschätzungen werden im Gegensatz zu den Boxplots multimodale Verteilungen erkennbar.

# b) 
model_full <- lm(formula = zscore ~ k_geschl + k_still + k_alter + m_alterg + m_groesse + m_bmi,
                 data = sambia)
summary(model_full)


# c) 
# Berechnung über Package car
A <- matrix(data = 0, nrow = 6, ncol = 7)
diag(A[, 2:7]) <- 1
c <- 0
linearHypothesis(model = model_full, hypothesis.matrix = A, rhs = c,
                 test = "F")

# Manuelle Berechnung:
# Definition der Modellparameter:
X <- model.matrix(model_full)
beta_hat <- as.vector(model_full$coefficients)
p <- length(beta_hat)-1
A <- cbind(rep(0, p), diag(p))
c <- rep(0, p)
a <- nrow(A)
n <- nrow(sambia)
# Berechnung der Teststatistik:
SSH <- as.vector(t(A %*% beta_hat - c) %*%
                   solve(A %*% solve(t(X) %*% X) %*% t(A)) %*%
                   (A %*% beta_hat - c))
SSE <- as.vector(sum(model_full$residuals^2))
TF <- (SSH / a) / (SSE / (nrow(sambia) - length(beta_hat)))
TF
# Testentscheidung:
TF > qf(p = 0.95, df1 = a, df2 = n - p - 1)

# Die Nullhypothese kann zum Signikanzniveau α = 0.05 abgelehnt werden.
# Mindestens eine der betrachteten Kovariablen weist also einen signifikanten Zusammenhang
# mit dem Z-score auf.



# d) 
# package car
A <- c(0, 0, 1, 0, 0, 0, 0)
c <- 0
linearHypothesis(model = model_full, hypothesis.matrix = A, rhs = c,
                 test = "F")

# Vergleich der Quadratsummen über Funktion anova:
model_d <- lm(formula = zscore ~ k_geschl + k_alter + m_alterg + m_groesse +
                  m_bmi,
                data = sambia)
summary(model_d)
anova(model_d, model_full)

# Die Quadratsumme des vollen Modells wird hierbei mit der Quadratsumme des
# gemäß der Nullhypothese restringierten Modells verglichen, d.h. in diesem Fall mit
# der Quadratsumme des Modells ohne die Stilldauer.


# Standardoutput des lm-Objektes (t-Test):
summary(model_full)
summary(model_full)$coef["k_still", ]
summary(model_full)$coef["k_still", "t value"]^2


# e) 
# Package car
A <- rbind(c(0, 1, 0, 0, 0, 0, 0), c(0, 0, 0, 0, 1, 0, 0))
c <- c(0, 0)
linearHypothesis(model = model_full, hypothesis.matrix = A, rhs = c,
                 test = "F")

# Vergleich der Quadratsummen:
model_e <- lm(formula = zscore ~ k_still + k_alter + m_groesse + m_bmi,
                data = sambia)
summary(model_e)
anova(model_e, model_full)

# Die Nullhypothese kann zum Signifikanzniveau α = 0.05 nicht abgelehnt werden.


# f)
# Package car
A <- c(0, 0, 0, 0, 0, 1, -1)
c <- 0
linearHypothesis(model = model_full, hypothesis.matrix = A, rhs = c,
                 test = "F")

# Vergleich der Quadratsummen:
model_f <- lm(formula = zscore ~ k_geschl + k_still  + k_alter + m_alterg +
                  I(m_groesse + m_bmi),
                data = sambia)
summary(model_f)
anova(model_f, model_full)

# Die Nullhypothese kann zum Signifikanzniveau α = 0.05 nicht abgelehnt werden
# Ein gleich großer Zusammenhang der Größe der Mutter und deren BMI ist also möglich.


# g) 
# package car
A <- c(0, 0, 0, 1, 0, 0, 0)
c <- -2
linearHypothesis(model = model_full, hypothesis.matrix = A, rhs = c,
                 test = "F")

# Vergleich der Quadratsummen:
model_g <- lm(formula = zscore ~ k_geschl + k_still + m_alterg + m_groesse +
                  m_bmi + offset(-2 * k_alter),
                data = sambia)
summary(model_g)
anova(model_g, model_full)

# Die Nullhypothese kann zum Signifikanzniveau α = 0.05 nicht abgelehnt werden.
# Der Effekt bezüglich des Kindesalters ist also nicht signifikant vom Wert -2 verschieden.