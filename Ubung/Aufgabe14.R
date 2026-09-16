library(tidyverse)
library(nlme)

############################# (a) ##############################################
data <- read.table(file = "../Daten/zufriedenheit.txt", header = TRUE)
head(data)
str(data)
summary(data)

data <- data %>% mutate(origin = as.factor(origin), gender = as.factor(sex))

plot1 <- ggplot(data = data, mapping = aes(x = origin, y = wellbeing)) +
  geom_boxplot() +
  theme_minimal() +
  ylab("Lebenszufriedenheit") +
  xlab("Herkunft") +
  ggtitle("Herkunft")
plot2 <- ggplot(data = data, mapping = aes(x = as.factor(sex), y = wellbeing)) +
  geom_boxplot() +
  theme_minimal() +
  ylab("Lebenszufriedenheit") +
  xlab("Geschlecht") +
  ggtitle("Geschlecht")
plot3 <- ggplot(data = data, mapping = aes(
  x = as.factor(year),
  y = wellbeing
)) +
  geom_boxplot() +
  theme_minimal() +
  ylab("Lebenszufriedenheit") +
  xlab("Jahr") +
  ggtitle("Jahr") +
  scale_x_discrete(breaks = NULL)
plot4 <- ggplot(data = data, mapping = aes(
  x = as.factor(age),
  y = wellbeing
)) +
  geom_boxplot() +
  theme_minimal() +
  ylab("Lebenszufriedenheit") +
  xlab("Alter") +
  ggtitle("Alter") +
  scale_x_discrete(breaks = NULL)
ggarrange(plot1, plot2, plot3, plot4)


############################# (b) ##############################################
regspline <- function(kovariable, data, K = 3, g = 3) {
  kov <- data[[kovariable]]
  knoten <- seq(min(kov), max(kov), length = K + 2)
  X <- matrix(data = 0, nrow = length(kov), ncol = K + g + 1)
  # Polynome
  for (i in 0:g) {
    X[, i + 1] <- kov^i
  }
  # Trunkierte Potenzen
  for (i in 2:(K + 1)) {
    X[, i + g] <- (kov > knoten[i]) * pmax((kov - knoten[i]), 0)^g
  }
  X <- X[, -1]
  colnames(X) <- c(
    paste0(kovariable, c("", seq_len(g)[-1])),
    paste0(kovariable, "_tp", seq_len(g))
  )
  result <- list(X = X, knoten = knoten)
  return(result)
}

X_reg <- regspline(kovariable = "age", data = data, K = 3, g = 3)$X
X <- cbind(data[, c("id", "wellbeing", "origin", "sex", "year")], X_reg)
head(X)

# Modellierung
model <- lm(
  formula = wellbeing ~ origin + sex + year + age + age2 + age3 +
    age_tp1 + age_tp2 + age_tp3,
  data = X
)
summary(model)

# Visualisierung
data_new <- data[
  rep(1, 100),
  c("id", "wellbeing", "origin", "sex", "year")
]
data_new$age <- seq(from = min(data$age), to = max(data$age), length.out = 100)
reg_new <- regspline(kovariable = "age", data = data_new, K = 3, g = 3)$X
data_new <- cbind(data_new, reg_new[, -1])
data_new$lm <- predict(object = model, newdata = data_new)
data_new <- data_new %>% gather(model, estimate, lm)

plot <- ggplot() +
  geom_point(data = data, mapping = aes(x = age, y = wellbeing)) +
  geom_line(data = data_new, aes(x = age, y = estimate), col = "blue") +
  theme_minimal() +
  xlab("age") +
  ylab("wellbeing")
plot

data <- data %>%
  mutate(fitted = model$fitted.values, residuals = stdres(model))
ggplot(
  data = data %>% filter(id %in% 1:5),
  mapping = aes(x = year, y = residuals, col = as.factor(id))
) +
  geom_line(lwd = 0.3) +
  xlab("Jahr") +
  ylab("Standardisierte Residuen") +
  theme_minimal() +
  theme(legend.position = "none")


############################# (d) ##############################################
model_ar <- gls(
  model = wellbeing ~ origin + sex + year + age + age2 + age3 +
    age_tp1 + age_tp2 + age_tp3,
  data = X, correlation = corAR1(form = ~ 1 | id)
)
summary(model_ar)

phi <- model_ar$modelStruct$corStruct
phi

model_list <- list()
model_list[["model_lm"]] <- summary(model)$coefficients[2:4, 1:2]
model_list[["model_ar"]] <- summary(model_ar)$tTable[2:4, 1:2]
for (i in seq_along(model_list)) {
  colnames(model_list[[i]]) <- c("Estimate", "SE")
}
model_overview <- melt(model_list) %>%
  spread(Var2, value) %>%
  rename(variable = Var1, model = L1) %>%
  mutate(lower = Estimate - 1.96 * SE, upper = Estimate + 1.96 * SE)
model_overview


ggplot(
  data = model_overview,
  mapping = aes(x = variable, y = Estimate, col = model)
) +
  geom_pointrange(
    mapping = aes(ymin = lower, ymax = upper),
    position = position_dodge(w = 0.7)
  ) +
  facet_wrap(~variable) +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    axis.text.y = element_blank()
  )

data_new <- data[
  rep(1, 100),
  c("id", "wellbeing", "origin", "sex", "year")
]
data_new$age <- seq(from = min(data$age), to = max(data$age), length.out = 100)
reg_new <- regspline(kovariable = "age", data = data_new, K = 3, g = 3)$X
data_new <- cbind(data_new, reg_new[, -1])
data_new$lm <- predict(object = model, newdata = data_new)
data_new$ar <- predict(object = model_ar, newdata = data_new)
data_new <- data_new %>% gather(model, estimate, lm:ar)
ggplot(data = data_new, mapping = aes(x = age, y = estimate, col = model)) +
  geom_line(lwd = 0.8) +
  theme_minimal() +
  theme(legend.title = element_blank())
