library(tidyverse)

mod <- lm(mpg ~ wt, data = mtcars)

round(mod$coefficients, 1)

# scatterplot ohne Loesung
scatter <- ggplot(data = mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 6), name = "Gewicht in 1000 Pfund",
                     breaks = seq(0, 6, 1),) +
  scale_y_continuous(limits = c(0, 40), name = "Reichweite in Meilen pro Gallone",
                     breaks = seq(0, 40, 10)) +
  theme_bw()
scatter

# scatterplot mit Loesung
dat <- data.frame(mpg = mod$model$mpg, 
                  wt = mod$model$wt,
                  pred = mod$fitted.values)[c(17, 20, 29),]
scatter_extra <- scatter +
  # Regressionsgerade
  geom_abline(intercept = mod$coefficients[1],
              slope = mod$coefficients[2],
              size = 1) +
  # Beobachtungen x_i, y_i
  geom_point(data = dat, size = 5, shape = 1, col = "blue") +
  # Residuen epsilon_i
  geom_segment(data = dat,
               aes(x = wt, 
                   xend = wt,
                   y = pred, 
                   yend = mpg), col = "orange") +
  # Regressionskoeffizienten
  # Intercept
  geom_segment(aes(x = 0,
                   xend = 0,
                   y = 0,
                   yend = mod$coefficients[1]),
               size = 1,
               col = "darkgreen") +
  # Slope
  geom_segment(x = 2,
               xend = 3,
               y = predict(mod, new = data.frame(wt = c(2))),
               yend = predict(mod, new = data.frame(wt = c(2))),
               col = "darkgreen") +
  geom_segment(x = 3,
               xend = 3,
               y = predict(mod, new = data.frame(wt = c(2))),
               yend = predict(mod, new = data.frame(wt = c(3))),
               size = 1,
               col = "darkgreen") +
  # Beschriftung
  # Regressionsgerade
  annotate("Text", 
           x = c(0.15, 3.25), 
           y = c(18,24), 
           label = c(expression(hat(beta)[0]), expression(hat(beta)[gewicht])),
           col = "darkgreen") +
  # Beobachtungen
  annotate("Text", 
           x = c(5.8, 2.25, 2.75), 
           y = c(14.8, 34.7, 15.8), 
           label = c("Beob. 17", "Beob. 20", "Beob. 29"),
           col = "blue") +
  # Residuen
  annotate("Text", 
           x = c(5.5, 2, 3), 
           y = c(12, 30.5, 18.3), 
           label = c(expression(epsilon[17]),
                     expression(epsilon[20]),
                     expression(epsilon[29])),
           col = "orange",
           size = 5)

scatter_extra

# Wir wollen als nächstes besser verstehen, woher die unsicherheit in unserem
# Modell kommt. Dazu verwenden wir Bootstrap, um die Unsicherheit in den
# Regressionskoeffizienten zu quantifizieren.
# Wir vergleichen das dann mit den uns bekannten Konfidenzintervallen, die wir
# mit geom_smooth() erhalten können.

set.seed(12)
n_boot <- 100

# Bootstrap models
boot_models <- map(1:n_boot, ~ {
  sample_n(mtcars, size = nrow(mtcars), replace = TRUE) %>%
    lm(mpg ~ wt, data = .)
})

# Extract intercept and slope for each model
boot_coefs <- map_dfr(boot_models, ~{
  coefs <- coef(.x)
  tibble(intercept = coefs[1], slope = coefs[2])
}, .id = "boot_id")  # .id gives a unique model number

# Original model
mod <- lm(mpg ~ wt, data = mtcars)

# Plot
scatter_extra + 
  geom_abline(data = boot_coefs, aes(intercept = intercept, slope = slope), 
              alpha = 0.2, color = "gray") +
  geom_abline(intercept = coef(mod)[1], slope = coef(mod)[2], 
              size = 1.2, color = "red") +
  geom_point(data = mtcars, aes(x = wt, y = mpg)) +
  geom_smooth(data = mtcars, aes(x = wt, y = mpg), method = "lm", se = TRUE, 
              size = 1.2, lty = "dashed") +
  theme_bw() +
  labs(title = "Bootstrapped Regressionsgeraden und Ursprungsmodell\nverglichen mit geom_smooth()")