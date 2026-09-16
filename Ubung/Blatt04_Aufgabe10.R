library(faraway)
library(tidyverse)
library(broom)
library(effects)

data(teengamb)
head(teengamb)
str(teengamb)
summary(teengamb)

teengamb <- teengamb %>% mutate(sex = as.factor(sex))


######################### (a) ##################################################
model_main <- lm(formula = gamble ~ income + sex, data = teengamb)
summary(model_main)

# Viz
fitted_main <- augment(model_main)
ggplot(data = fitted_main, mapping = aes(x = income, y = gamble, col = sex)) +
  geom_point() +
  geom_line(aes(y = .fitted)) +
  theme_minimal()


######################### (b) ##################################################
model_int <- lm(formula = gamble ~ income * sex, data = teengamb)
summary(model_int)

# Viz
fitted_int <- augment(model_int)
ggplot(data = fitted_int, mapping = aes(x = income, y = gamble, col = sex)) +
  geom_point() +
  geom_line(aes(y = .fitted)) +
  theme_minimal()

# Grafischer Vergleich der Residuen
fitted_main <- fitted_main %>% mutate(model = "No Interaction")
fitted_int <- fitted_int %>% mutate(model = "Interaction")
fitted <- rbind(fitted_main, fitted_int)
ggplot(data = fitted, mapping = aes(x = income, y = gamble, col = sex)) +
  geom_point() +
  geom_line(aes(y = .fitted)) +
  theme_minimal() +
  geom_segment(aes(xend = income, yend = .fitted, col = sex)) +
  facet_wrap(~model)


######################### (c) ##################################################

# ohne Interaktion
unname(coef(model_main)[2] * 10 + coef(model_main)[3])

# mit Interaktion
unname(coef(model_int)[2] * 10 + coef(model_int)[3] + coef(model_int)[4] * 12) # income = 2
unname(coef(model_int)[2] * 10 + coef(model_int)[3] + coef(model_int)[4] * 14) # income = 4


######################### (d) ##################################################
model_full <- lm(formula = gamble ~ (income + status + verbal) * sex, data = teengamb)
summary(model_full)

plot(x = allEffects(model_full), rows = 3, cols = 1, ylim = c(-50, 150))


######################### (e) ##################################################
anova(model_int, model_full)
