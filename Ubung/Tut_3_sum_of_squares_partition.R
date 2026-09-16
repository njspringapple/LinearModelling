library(tidyverse)
library(ggpubr)

set.seed(1022)
n <- 15
x <- rnorm(n, mean = 5, sd = 1)
y <- 1 + 2 * x + rnorm(n, sd = 1)
z <- 1 + 0.1 * x + rnorm(n, sd = 1)
avg_y <- mean(y)
lm_model <- lm(y ~ x)
lm_model_z <- lm(z ~ x)
summary(lm_model)
summary(lm_model_z)

ggplot(data = tibble(x, y), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = avg_y, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Scatter plot with linear regression line",
       x = "x",
       y = "Y") +
  theme_bw()

# draw a rectangle with the previous line being it's height
gg_sst <- ggplot(data = tibble(x, y), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = avg_y, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_segment(aes(xend = x, yend = avg_y), color = "black", linetype = "dashed", size = 1.2) +
  geom_rect(aes(xmin = pmin(x, x - (avg_y - y)), xmax = pmax(x, x - (avg_y - y)), ymin = pmin(y, avg_y), ymax = pmax(y, avg_y)), fill = "red", alpha = 0.5) +
  labs(title = "(SST) Squared error of y to mean(y)",
       x = "x",
       y = "Y") +
  scale_x_continuous(limits = c(-3, 14)) +
  theme_bw() + 
  coord_equal()
gg_sst

# draw the lines now from the predicted y value to the average y value and adjust the boxes
gg_ssm <- ggplot(data = tibble(x, y), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = avg_y, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_segment(aes(xend = x, y = avg_y, yend = predict(lm_model)), color = "black", linetype = "dashed", size = 1.2) +
  geom_rect(aes(xmin = pmin(x, x - (avg_y - predict(lm_model))), xmax = pmax(x, x - (avg_y - predict(lm_model))), ymin = pmin(predict(lm_model), avg_y), ymax = pmax(predict(lm_model), avg_y)), fill = "blue", alpha = 0.5) +
  labs(title = "(SSM) Squared error of hat{y} to mean(y)",
       x = "x",
       y = "Y") +
  scale_x_continuous(limits = c(-3, 14)) +
  theme_bw() + 
  coord_equal()
gg_ssm

# draw the lines from the predicted y value to the point and adjust the boxes
gg_sse <- ggplot(data = tibble(x, y), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = avg_y, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_segment(aes(xend = x, yend = predict(lm_model)), color = "black", linetype = "dashed", size = 1.2) +
  geom_rect(aes(xmin = pmin(x, x - (predict(lm_model) - y)), xmax = pmax(x, x - (predict(lm_model) - y)), ymin = pmin(y, predict(lm_model)), ymax = pmax(y, predict(lm_model))), fill = "green", alpha = 0.5) +
  labs(title = "(SSE) Squared error of y to hat(y)",
       x = "x",
       y = "Y") +
  scale_x_continuous(limits = c(-3, 14)) +
  theme_bw() + 
  coord_equal()
gg_sse

ggarrange(gg_sst, gg_ssm, gg_sse, 
          ncol = 3, nrow = 1)


# Ohne Intercept
set.seed(1022)
n <- 15
x <- rnorm(n, mean = 5, sd = 1)
y <- 2 * x + rnorm(n, sd = 1)
avg_y <- 0
lm_model <- lm(y ~ - 1 + x)

ggplot(data = tibble(x, y), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = avg_y, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Scatter plot with linear regression line",
       x = "x",
       y = "Y") +
  theme_bw()

# draw a rectangle with the previous line being it's height
ggplot(data = tibble(x, y), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = avg_y, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, formula = "y ~ -1 + x") +
  geom_segment(aes(xend = x, yend = avg_y), color = "black", linetype = "dashed", size = 1.2) +
  geom_rect(aes(xmin = pmin(x, x - (avg_y - y)), xmax = pmax(x, x - (avg_y - y)), ymin = pmin(y, avg_y), ymax = pmax(y, avg_y)), fill = "red", alpha = 0.5) +
  labs(title = "(SST*) Squared error of y to 0",
       x = "x",
       y = "Y") +
  theme_bw() + 
  coord_equal()

# draw the lines now from the predicted y value to the average y value and adjust the boxes
ggplot(data = tibble(x, y), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = avg_y, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, formula = "y ~ -1 + x") +
  geom_segment(aes(xend = x, y = avg_y, yend = predict(lm_model)), color = "black", linetype = "dashed", size = 1.2) +
  geom_rect(aes(xmin = pmin(x, x - (avg_y - predict(lm_model))), xmax = pmax(x, x - (avg_y - predict(lm_model))), ymin = pmin(predict(lm_model), avg_y), ymax = pmax(predict(lm_model), avg_y)), fill = "blue", alpha = 0.5) +
  labs(title = "(SSM*) Squared error of hat{y} to 0",
       x = "x",
       y = "Y") +
  theme_bw() + 
  coord_equal()

# draw the lines from the predicted y value to the point and adjust the boxes
ggplot(data = tibble(x, y), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = avg_y, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, formula = "y ~ 0 + x") +
  geom_segment(aes(xend = x, yend = predict(lm_model)), color = "black", linetype = "dashed", size = 1.2) +
  geom_rect(aes(xmin = pmin(x, x - (predict(lm_model) - y)), xmax = pmax(x, x - (predict(lm_model) - y)), ymin = pmin(y, predict(lm_model)), ymax = pmax(y, predict(lm_model))), fill = "green", alpha = 0.5) +
  labs(title = "(SSE*) Squared error of y to hat(y)",
       x = "x",
       y = "Y") +
  theme_bw() + 
  coord_equal()



# Schlechtes Model verwenden. Dann ist R^2 negativ
set.seed(1022)
n <- 15
x <- rnorm(n, mean = 5, sd = 1)
y <- 1 + 2 * x + rnorm(n, sd = 1)
avg_y <- mean(y)

ggplot(data = tibble(x, y), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = avg_y, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "(SST) Scatter plot with linear regression line",
       x = "x",
       y = "Y") +
  theme_bw()

# draw a rectangle with the previous line being it's height
lm_model <- -x * sin(x) + avg_y
ggplot(data = tibble(x, y), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = avg_y, color = "black", linetype = "dashed") +
  geom_line(aes(y = lm_model)) +
  geom_segment(aes(xend = x, yend = avg_y), color = "black", linetype = "dashed", size = 1.2) +
  geom_rect(aes(xmin = pmin(x, x - (avg_y - y)), xmax = pmax(x, x - (avg_y - y)), ymin = pmin(y, avg_y), ymax = pmax(y, avg_y)), fill = "red", alpha = 0.5) +
  labs(title = "Squared error of y to mean(y)",
       x = "x",
       y = "Y") +
  theme_bw() + 
  coord_equal()

# draw the lines now from the predicted y value to the average y value and adjust the boxes
ggplot(data = tibble(x, y), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = avg_y, color = "black", linetype = "dashed") +
  geom_line(aes(y = lm_model)) +
  geom_segment(aes(xend = x, y = avg_y, yend = lm_model), color = "black", linetype = "dashed", size = 1.2) +
  geom_rect(aes(xmin = pmin(x, x - (avg_y - lm_model)), xmax = pmax(x, x - (avg_y - lm_model)), ymin = pmin(lm_model, avg_y), ymax = pmax(lm_model, avg_y)), fill = "blue", alpha = 0.5) +
  labs(title = "(SSM) Squared error of hat{y} to mean(y)",
       x = "x",
       y = "Y") +
  theme_bw() + 
  coord_equal()

# draw the lines from the predicted y value to the point and adjust the boxes
ggplot(data = tibble(x, y), aes(x = x, y = y)) +
  geom_point() +
  geom_hline(yintercept = avg_y, color = "black", linetype = "dashed") +
  geom_line(aes(y = lm_model)) +
  geom_segment(aes(xend = x, yend = lm_model), color = "black", linetype = "dashed", size = 1.2) +
  geom_rect(aes(xmin = pmin(x, x - (lm_model - y)), xmax = pmax(x, x - (lm_model - y)), ymin = pmin(y, lm_model), ymax = pmax(y, lm_model)), fill = "green", alpha = 0.5) +
  labs(title = "(SSE) Squared error of y to hat(y)",
       x = "x",
       y = "Y") +
  theme_bw() + 
  coord_equal()

