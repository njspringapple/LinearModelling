data(cars)

str(cars)

plot(dist~speed, data=cars)

fit_cars <- lm(dist~speed, data=cars)
abline(fit_cars)
coef(fit_cars)

# 一个离群点就能把直线"拽"过去——平方损失放大了大残差的影响
fit1 <- lm(dist~speed, data=cars)
coef(fit1)

cars2 <- rbind(cars, data.frame(speed = 25, dist = 300))
fit2 <- lm(dist~speed, data=cars2)
coef(fit2)

plot(dist~speed, data=cars2)
abline(fit1,col="blue")
abline(fit2,col="red")







# OLS方差的三个因素实验
# 运行这段代码观察不同因素对估计精度的影响

set.seed(123)

# 真实参数
beta0_true <- 84
beta1_true <- 4

# 模拟函数：重复实验很多次，看beta1估计值的分布
simulate_beta1 <- function(n, sigma, x_range, n_sim = 1000) {
  beta1_estimates <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # 生成x（在指定范围内均匀分布）
    x <- runif(n, min = x_range[1], max = x_range[2])
    # 生成y
    epsilon <- rnorm(n, mean = 0, sd = sigma)
    y <- beta0_true + beta1_true * x + epsilon
    # 估计beta1
    beta1_estimates[i] <- cov(x, y) / var(x)
  }
  
  return(beta1_estimates)
}

# ============================================
# 实验1：sigma的影响
# ============================================
par(mfrow = c(2, 2))

sigma_small <- simulate_beta1(n = 30, sigma = 2, x_range = c(4, 12))
sigma_large <- simulate_beta1(n = 30, sigma = 10, x_range = c(4, 12))

hist(sigma_small, breaks = 30, col = "lightblue", 
     main = "sigma = 2 (噪声小)", xlab = "beta1估计值", xlim = c(0, 8))
abline(v = beta1_true, col = "red", lwd = 2)
abline(v = mean(sigma_small), col = "blue", lwd = 2, lty = 2)
legend("topright", c("真实值", "估计均值"), col = c("red", "blue"), lty = c(1, 2), cex = 0.8)

hist(sigma_large, breaks = 30, col = "lightcoral",
     main = "sigma = 10 (噪声大)", xlab = "beta1估计值", xlim = c(0, 8))
abline(v = beta1_true, col = "red", lwd = 2)
abline(v = mean(sigma_large), col = "blue", lwd = 2, lty = 2)

cat("实验1：sigma的影响\n")
cat("sigma=2 时 V(beta1) =", var(sigma_small), "\n")
cat("sigma=10时 V(beta1) =", var(sigma_large), "\n\n")

# ============================================
# 实验2：样本量n的影响
# ============================================
n_small <- simulate_beta1(n = 10, sigma = 5, x_range = c(4, 12))
n_large <- simulate_beta1(n = 100, sigma = 5, x_range = c(4, 12))

hist(n_small, breaks = 30, col = "lightgreen",
     main = "n = 10 (样本少)", xlab = "beta1估计值", xlim = c(0, 8))
abline(v = beta1_true, col = "red", lwd = 2)

hist(n_large, breaks = 30, col = "lightgoldenrod",
     main = "n = 100 (样本多)", xlab = "beta1估计值", xlim = c(0, 8))
abline(v = beta1_true, col = "red", lwd = 2)

cat("实验2：样本量n的影响\n")
cat("n=10 时 V(beta1) =", var(n_small), "\n")
cat("n=100时 V(beta1) =", var(n_large), "\n\n")

# ============================================
# 实验3：x变异范围的影响
# ============================================
par(mfrow = c(1, 2))

x_narrow <- simulate_beta1(n = 30, sigma = 5, x_range = c(7, 9))   # 窄范围
x_wide <- simulate_beta1(n = 30, sigma = 5, x_range = c(2, 14))    # 宽范围

hist(x_narrow, breaks = 30, col = "plum",
     main = "x范围 [7,9] (窄)", xlab = "beta1估计值", xlim = c(-5, 15))
abline(v = beta1_true, col = "red", lwd = 2)

hist(x_wide, breaks = 30, col = "skyblue",
     main = "x范围 [2,14] (宽)", xlab = "beta1估计值", xlim = c(-5, 15))
abline(v = beta1_true, col = "red", lwd = 2)

cat("实验3：x变异范围的影响\n")
cat("x在[7,9]时 V(beta1) =", var(x_narrow), "\n")
cat("x在[2,14]时 V(beta1) =", var(x_wide), "\n\n")

# ============================================
# 可视化：一次模拟中窄范围 vs 宽范围的区别
# ============================================
par(mfrow = c(1, 2))

set.seed(456)

# 窄范围
x1 <- runif(30, 7, 9)
y1 <- beta0_true + beta1_true * x1 + rnorm(30, 0, 5)
plot(x1, y1, pch = 19, col = "plum", xlim = c(0, 16), ylim = c(80, 160),
     main = "x范围窄: 斜率难判断", xlab = "x", ylab = "y")
abline(lm(y1 ~ x1), col = "purple", lwd = 2)
abline(a = beta0_true, b = beta1_true, col = "red", lwd = 2, lty = 2)
legend("topleft", c("估计线", "真实线"), col = c("purple", "red"), lty = c(1, 2), cex = 0.8)

# 宽范围
x2 <- runif(30, 2, 14)
y2 <- beta0_true + beta1_true * x2 + rnorm(30, 0, 5)
plot(x2, y2, pch = 19, col = "skyblue", xlim = c(0, 16), ylim = c(80, 160),
     main = "x范围宽: 斜率清晰", xlab = "x", ylab = "y")
abline(lm(y2 ~ x2), col = "blue", lwd = 2)
abline(a = beta0_true, b = beta1_true, col = "red", lwd = 2, lty = 2)
legend("topleft", c("估计线", "真实线"), col = c("blue", "red"), lty = c(1, 2), cex = 0.8)









# 展示权重wi的含义：xi离均值越远，对beta1估计贡献越大

set.seed(123)

# ============================================
# 示例数据
# ============================================
x <- c(4, 6, 8, 10, 12)
y <- c(100, 110, 115, 125, 130)
n <- length(x)

# 计算权重
x_bar <- mean(x)
w <- (x - x_bar) / sum((x - x_bar)^2)

# 展示权重
cat("============================================\n")
cat("每个数据点的权重\n")
cat("============================================\n")
data.frame(
  i = 1:n,
  x = x,
  y = y,
  "x - x_bar" = x - x_bar,
  "权重w" = round(w, 4)
)

# 验证：加权和等于beta1
beta1_weighted <- sum(w * y)
beta1_ols <- cov(x, y) / var(x)

cat("\n用加权和计算: beta1 =", beta1_weighted, "\n")
cat("用OLS公式计算: beta1 =", beta1_ols, "\n")

# ============================================
# 可视化1：权重分布
# ============================================
par(mfrow = c(2, 2))

# 图1：权重柱状图
barplot(w, names.arg = paste0("x=", x), col = ifelse(w > 0, "steelblue", "coral"),
        main = "每个点的权重 wi", ylab = "权重", border = NA)
abline(h = 0, lty = 2)
text(1:5 * 1.2 - 0.5, w + sign(w) * 0.01, round(w, 3), cex = 0.9)

# 图2：数据点，点大小表示|权重|
plot(x, y, pch = 19, cex = abs(w) * 50 + 1, col = ifelse(w > 0, "steelblue", "coral"),
     main = "点大小 = |权重|", xlab = "x", ylab = "y", xlim = c(2, 14))
abline(v = x_bar, lty = 2, col = "gray")
text(x_bar, 105, "x均值", pos = 4, col = "gray40")
abline(lm(y ~ x), col = "black", lwd = 2)

# ============================================
# 可视化2：移动一个点对beta1的影响
# ============================================

# 图3：移动边缘点 vs 移动中间点
# 把y[1]（x=4，边缘点）加10
y_move_edge <- y
y_move_edge[1] <- y[1] + 10
beta1_move_edge <- cov(x, y_move_edge) / var(x)

# 把y[3]（x=8，中间点）加10
y_move_center <- y
y_move_center[3] <- y[3] + 10
beta1_move_center <- cov(x, y_move_center) / var(x)

cat("\n============================================\n")
cat("移动不同位置的点对beta1的影响\n")
cat("============================================\n")
cat("原始 beta1 =", round(beta1_ols, 4), "\n")
cat("边缘点(x=4)的y加10后: beta1 =", round(beta1_move_edge, 4), 
    ", 变化 =", round(beta1_move_edge - beta1_ols, 4), "\n")
cat("中间点(x=8)的y加10后: beta1 =", round(beta1_move_center, 4),
    ", 变化 =", round(beta1_move_center - beta1_ols, 4), "\n")

# 图3：移动边缘点
plot(x, y, pch = 19, col = "gray", cex = 2, xlim = c(2, 14), ylim = c(95, 145),
     main = "移动边缘点(x=4): beta1变化大", xlab = "x", ylab = "y")
points(x[1], y_move_edge[1], pch = 19, col = "coral", cex = 2)
arrows(x[1], y[1], x[1], y_move_edge[1] - 2, col = "coral", lwd = 2, length = 0.1)
abline(lm(y ~ x), col = "black", lwd = 2)
abline(lm(y_move_edge ~ x), col = "coral", lwd = 2, lty = 2)
legend("topleft", c("原始线", "移动后"), col = c("black", "coral"), lty = c(1, 2), lwd = 2)

# 图4：移动中间点
plot(x, y, pch = 19, col = "gray", cex = 2, xlim = c(2, 14), ylim = c(95, 145),
     main = "移动中间点(x=8): beta1几乎不变", xlab = "x", ylab = "y")
points(x[3], y_move_center[3], pch = 19, col = "steelblue", cex = 2)
arrows(x[3], y[3], x[3], y_move_center[3] - 2, col = "steelblue", lwd = 2, length = 0.1)
abline(lm(y ~ x), col = "black", lwd = 2)
abline(lm(y_move_center ~ x), col = "steelblue", lwd = 2, lty = 2)
legend("topleft", c("原始线", "移动后"), col = c("black", "steelblue"), lty = c(1, 2), lwd = 2)

# ============================================
# 可视化3：正态性传递
# ============================================
par(mfrow = c(1, 2))

n_sim <- 5000
beta1_estimates <- numeric(n_sim)

for (i in 1:n_sim) {
  epsilon <- rnorm(5, mean = 0, sd = 3)
  y_sim <- 84 + 4 * x + epsilon
  beta1_estimates[i] <- sum(w * y_sim)
}

# 直方图
hist(beta1_estimates, breaks = 50, probability = TRUE, col = "lightblue",
     main = "beta1估计值的分布", xlab = "beta1", border = "white")

# 叠加理论正态曲线
theoretical_var <- 9 / sum((x - x_bar)^2)  # sigma^2 / sum(xi - xbar)^2
curve(dnorm(x, mean = 4, sd = sqrt(theoretical_var)), add = TRUE, col = "red", lwd = 2)
legend("topright", "理论正态分布", col = "red", lwd = 2)

# QQ图
qqnorm(beta1_estimates, main = "QQ图: 验证正态性", pch = 19, col = "steelblue", cex = 0.5)
qqline(beta1_estimates, col = "red", lwd = 2)

cat("\n============================================\n")
cat("正态性验证\n")
cat("============================================\n")
cat("beta1估计值的均值:", round(mean(beta1_estimates), 4), "(理论值: 4)\n")
cat("beta1估计值的方差:", round(var(beta1_estimates), 4), 
    "(理论值:", round(theoretical_var, 4), ")\n")



# 设置画布
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# ===================
# 图1：n=2 时，两点确定一条直线，残差=0
# ===================
x1 <- c(1, 3)
y1 <- c(2, 6)
fit1 <- lm(y1 ~ x1)

plot(x1, y1, pch = 19, cex = 2, col = "blue",
     xlim = c(0, 4), ylim = c(0, 8),
     xlab = "x", ylab = "y",
     main = "n=2: 零自由度\n直线完美穿过两点")
abline(fit1, col = "red", lwd = 2)
text(2, 7, paste("残差平方和 =", round(sum(resid(fit1)^2), 4)), cex = 1.2)
text(2, 6, "自由度 = 2-2 = 0", cex = 1.2, col = "red")

# ===================
# 图2：n=3 时，只有1个自由度
# ===================
x2 <- c(1, 2, 3)
y2 <- c(2, 5, 4)
fit2 <- lm(y2 ~ x2)

plot(x2, y2, pch = 19, cex = 2, col = "blue",
     xlim = c(0, 4), ylim = c(0, 7),
     xlab = "x", ylab = "y",
     main = "n=3: 一个自由度\n知道2个残差，第3个确定")
abline(fit2, col = "red", lwd = 2)
segments(x2, y2, x2, fitted(fit2), col = "green", lwd = 2, lty = 2)
res2 <- round(resid(fit2), 2)
text(x2 + 0.3, (y2 + fitted(fit2))/2, res2, col = "green", cex = 1.1)
text(2, 0.5, paste("残差:", paste(res2, collapse = " + "), "= 0"), cex = 1)
text(2, 6.5, "自由度 = 3-2 = 1", cex = 1.2, col = "red")

# ===================
# 图3：展示两个约束条件
# ===================
set.seed(42)
n <- 10
x3 <- 1:n
y3 <- 2 + 1.5 * x3 + rnorm(n, sd = 2)
fit3 <- lm(y3 ~ x3)
res3 <- resid(fit3)

plot(x3, y3, pch = 19, cex = 1.5, col = "blue",
     xlab = "x", ylab = "y",
     main = "n=10: 两个约束条件")
abline(fit3, col = "red", lwd = 2)
segments(x3, y3, x3, fitted(fit3), col = "green", lwd = 2, lty = 2)

# 显示约束条件
text(5, max(y3), 
     paste("约束1: Σε̂ =", round(sum(res3), 6)), cex = 1)
text(5, max(y3) - 1.5, 
     paste("约束2: Σxε̂ =", round(sum(x3 * res3), 6)), cex = 1)
text(5, min(y3), "自由度 = 10-2 = 8", cex = 1.2, col = "red")

# ===================
# 图4：模拟验证无偏性
# ===================
set.seed(123)
true_sigma2 <- 25  # 真实方差
n <- 30
n_sim <- 1000

estimates_n <- numeric(n_sim)    # 除以 n
estimates_n2 <- numeric(n_sim)   # 除以 n-2

for (i in 1:n_sim) {
  x <- runif(n, 0, 10)
  y <- 2 + 3 * x + rnorm(n, sd = 5)
  fit <- lm(y ~ x)
  res <- resid(fit)
  
  estimates_n[i] <- sum(res^2) / n
  estimates_n2[i] <- sum(res^2) / (n - 2)
}

hist(estimates_n2, breaks = 30, col = "lightblue", border = "white",
     main = "1000次模拟估计σ²",
     xlab = expression(hat(sigma)^2), xlim = c(0, 60))
abline(v = true_sigma2, col = "red", lwd = 3, lty = 1)
abline(v = mean(estimates_n), col = "orange", lwd = 2, lty = 2)
abline(v = mean(estimates_n2), col = "blue", lwd = 2, lty = 2)

legend("topright", 
       legend = c(
         paste("真实σ² =", true_sigma2),
         paste("除以n的均值 =", round(mean(estimates_n), 1), "(有偏)"),
         paste("除以n-2的均值 =", round(mean(estimates_n2), 1), "(无偏)")
       ),
       col = c("red", "orange", "blue"),
       lwd = 2, lty = c(1, 2, 2), cex = 0.9)







set.seed(123)

# 真实参数
beta0 <- 2
beta1 <- 3
sigma <- 5
n <- 20  # 样本量

# 模拟10000次
n_sim <- 10000
chi_sq_values <- numeric(n_sim)

x <- 1:n  # 固定 x

for (i in 1:n_sim) {
  # 生成数据
  epsilon <- rnorm(n, mean = 0, sd = sigma)
  Y <- beta0 + beta1 * x + epsilon
  
  # 拟合回归
  fit <- lm(Y ~ x)
  
  # 计算统计量
  residuals_sq_sum <- sum(residuals(fit)^2)
  chi_sq_values[i] <- residuals_sq_sum / sigma^2
}

# 画图比较
hist(chi_sq_values, breaks = 50, freq = FALSE, 
     main = "模拟 vs 理论 χ²(n-2) 分布",
     xlab = "统计量的值", col = "lightblue")

# 叠加理论分布
curve(dchisq(x, df = n - 2), add = TRUE, col = "red", lwd = 2)

legend("topright", legend = c("模拟结果", "χ²(18) 理论"),
       fill = c("lightblue", NA), border = c("black", NA),
       lty = c(NA, 1), lwd = c(NA, 2), col = c(NA, "red"))



set.seed(456)

# χ² 的定义：k 个独立标准正态的平方和
k <- 5
n_sim <- 10000

# 方法：生成 k 个 N(0,1)，平方后加起来
chi_sq_sim <- numeric(n_sim)

for (i in 1:n_sim) {
  Z <- rnorm(k)          # k 个标准正态
  chi_sq_sim[i] <- sum(Z^2)  # 平方和
}

# 画图
hist(chi_sq_sim, breaks = 50, freq = FALSE,
     main = "Z₁² + Z₂² + ... + Z₅² 服从 χ²(5)",
     xlab = "平方和", col = "lightgreen")

curve(dchisq(x, df = k), add = TRUE, col = "red", lwd = 2)



set.seed(789)

par(mfrow = c(2, 2))

for (df in c(2, 5, 10, 30)) {
  # 模拟
  sim <- replicate(10000, sum(rnorm(df)^2))
  
  hist(sim, breaks = 40, freq = FALSE,
       main = paste("χ² 分布, 自由度 =", df),
       xlab = "", col = "lightblue")
  curve(dchisq(x, df = df), add = TRUE, col = "red", lwd = 2)
}



set.seed(111)

beta0 <- 2; beta1 <- 3; sigma <- 5; n <- 30
x <- runif(n, 0, 10)

n_sim <- 10000

# 两种统计量
stat_n <- numeric(n_sim)    # 假装自由度是 n
stat_n2 <- numeric(n_sim)   # 正确的 n-2

for (i in 1:n_sim) {
  Y <- beta0 + beta1 * x + rnorm(n, sd = sigma)
  fit <- lm(Y ~ x)
  rss <- sum(residuals(fit)^2)
  
  stat_n[i] <- rss / sigma^2
  stat_n2[i] <- rss / sigma^2  # 同一个值，但比较不同理论分布
}

# 比较
par(mfrow = c(1, 2))

hist(stat_n, breaks = 50, freq = FALSE, 
     main = "与 χ²(n) 比较（错误）", col = "lightcoral")
curve(dchisq(x, df = n), add = TRUE, col = "blue", lwd = 2)

hist(stat_n2, breaks = 50, freq = FALSE,
     main = "与 χ²(n-2) 比较（正确）", col = "lightgreen")
curve(dchisq(x, df = n - 2), add = TRUE, col = "red", lwd = 2)









x <- seq(-4, 4, length = 200)

plot(x, dnorm(x), type = "l", lwd = 2, col = "blue",
     ylab = "密度", main = "t 分布 vs 标准正态")

lines(x, dt(x, df = 3), col = "red", lwd = 2)
lines(x, dt(x, df = 10), col = "orange", lwd = 2)
lines(x, dt(x, df = 30), col = "purple", lwd = 2)

legend("topright", 
       legend = c("N(0,1)", "t(3)", "t(10)", "t(30)"),
       col = c("blue", "red", "orange", "purple"), lwd = 2)



set.seed(123)

n <- 5          # 小样本
mu <- 0
sigma <- 1
n_sim <- 10000

# 两种统计量
Z_values <- numeric(n_sim)  # 用真实 σ
T_values <- numeric(n_sim)  # 用估计 s

for (i in 1:n_sim) {
  x <- rnorm(n, mean = mu, sd = sigma)
  xbar <- mean(x)
  s <- sd(x)
  
  Z_values[i] <- (xbar - mu) / (sigma / sqrt(n))  # 知道 σ
  T_values[i] <- (xbar - mu) / (s / sqrt(n))      # 用 s 估计
}

# 画图比较
par(mfrow = c(1, 2))

hist(Z_values, breaks = 50, freq = FALSE, main = "用真实 σ → N(0,1)",
     col = "lightblue", xlim = c(-5, 5))
curve(dnorm(x), add = TRUE, col = "red", lwd = 2)

hist(T_values, breaks = 50, freq = FALSE, main = "用估计 s → t(n-1)",
     col = "lightgreen", xlim = c(-5, 5))
curve(dt(x, df = n - 1), add = TRUE, col = "red", lwd = 2)












# 首先用数据集的数据集建模，估计出 beta0， beta1
fit_cars <- lm(dist ~ speed, data = cars)
new <- data.frame(speed = seq(4, 30, by = 1))

# 用建好的模型做预测

# 均值置信区间
ci <- predict(fit_cars, newdata = new, interval = "confidence")
# 个体预测区间
pi <- predict(fit_cars, newdata = new, interval = "prediction")

plot(dist ~ speed, data = cars, xlim = c(4, 30), ylim = c(-50, 200),
     main = "Konfidenz- und Prognoseintervalle")
abline(fit_cars, col = "black", lwd = 2)

# 均值置信区间（窄）
lines(new$speed, ci[, "lwr"], col = "blue", lty = 2, lwd = 2)
lines(new$speed, ci[, "upr"], col = "blue", lty = 2, lwd = 2)

# 预测区间（宽）
lines(new$speed, pi[, "lwr"], col = "red", lty = 3, lwd = 2)
lines(new$speed, pi[, "upr"], col = "red", lty = 3, lwd = 2)

legend("topleft",
       legend = c("KQ-Gerade", "95% KI Mittelwert", "95% Prognoseintervall"),
       col = c("black", "blue", "red"), lty = c(1, 2, 3), lwd = 2)













set.seed(42)
# 真实模型：Y = x² + ε
x <- c(0, 1, 2)
y <- x^2 + rnorm(3, sd = 0.3)

fit_lin <- lm(y ~ x)
fit_quad <- lm(y ~ x + I(x^2))

# 在 [-1, 3] 范围画图
x_plot <- seq(-1.5, 3, by = 0.01)
y_true <- x_plot^2
y_lin <- predict(fit_lin, newdata = data.frame(x = x_plot))

plot(x, y, xlim = c(-1.5, 3), ylim = c(-3, 10), pch = 19, cex = 1.5,
     main = "Lineare Näherung vs. wahre quadratische Funktion",
     xlab = "x", ylab = "y")
lines(x_plot, y_true, col = "black", lwd = 2)
lines(x_plot, y_lin, col = "red", lwd = 2)
legend("topleft",
       legend = c("wahre Funktion y = x²", "lineare Näherung"),
       col = c("black", "red"), lwd = 2)

# 外推到 x = -1
predict(fit_lin, newdata = data.frame(x = -1))
# 对比真实值 (-1)^2 = 1












set.seed(123)

# 真实参数
beta0_true <- 55
beta1_true <- 3.5

# 生成数据
x <- c(2, 4, 6, 8, 10)
epsilon <- rnorm(5, mean = 0, sd = 5)
Y <- beta0_true + beta1_true * x + epsilon

# 拟合模型
fit <- lm(Y ~ x)
beta0_hat <- coef(fit)[1]
beta1_hat <- coef(fit)[2]
Y_hat <- fitted(fit)
residuals <- Y - Y_hat

# 均值
x_bar <- mean(x)
Y_bar <- mean(Y)

# 画图
par(mar = c(5, 5, 4, 2))
plot(x, Y, 
     xlim = c(0, 12), ylim = c(50, 100),
     pch = 19, cex = 1.5, col = "blue",
     xlab = expression(x[i]),
     ylab = expression(Y[i]),
     main = "线性回归符号示意图")

# 真实直线
abline(a = beta0_true, b = beta1_true, col = "gray50", lty = 2, lwd = 2)

# 拟合直线
abline(fit, col = "red", lwd = 2)

# 残差线段
for(i in 1:5) {
  segments(x[i], Y[i], x[i], Y_hat[i], col = "darkgreen", lty = 3, lwd = 1.5)
}

# 预测值点
points(x, Y_hat, pch = 4, cex = 1.2, col = "red", lwd = 2)

# 均值点
points(x_bar, Y_bar, pch = 18, cex = 2, col = "purple")

# 图例（简化版）
# 图例（带数学符号版）
legend("topleft", 
       legend = c(
         expression(真实直线 ~ beta[0] + beta[1]*x[i]),
         expression(拟合直线 ~ hat(beta)[0] + hat(beta)[1]*x[i]),
         expression(观测点 ~ Y[i]),
         expression(预测点 ~ hat(Y)[i]),
         expression(残差 ~ hat(epsilon)[i]),
         expression(均值点 ~ bar(Y))
       ),
       col = c("gray50", "red", "blue", "red", "darkgreen", "purple"),
       lty = c(2, 1, NA, NA, 3, NA),
       pch = c(NA, NA, 19, 4, NA, 18),
       lwd = c(2, 2, NA, 2, 1.5, NA),
       cex = 0.8,
       bg = "white")

# 标注残差
text(x[3] + 0.6, (Y[3] + Y_hat[3])/2, expression(hat(epsilon)[i]), col = "darkgreen")

# 标注均值点
text(x_bar + 1, Y_bar + 3, expression(bar(Y)), col = "purple")


























data(cars)
head(cars)
fit <- lm(dist~speed,data=cars)
summary(fit)
# 提取基本信息
n <- nrow(cars)
n

# 计算残差的三种方法， 可以手工执行比较

# 方法1：residuals 函数
#e_hat <- residuals(fit)
#head(e_hat)
# 方法2：fitted 函数获取拟合值
#Y <- cars$dist
#Y_hat <- fitted(fit)
#e_hat <- Y - Y_hat
#head(e_hat)
# 方法3：根据公式手工算
Y <- cars$dist
x <- cars$speed
# 截距
beta0_hat <- coef(fit)[1]
# 斜率
beta1_hat <- coef(fit)[2]
Y_hat <- beta0_hat + beta1_hat * x
e_hat <- Y - Y_hat
head(e_hat)

# 计算残差平方和 SSE, 残差的单位是英尺，平方后是平方英尺
SSE <- sum(e_hat^2)
SSE

# 计算方差估计值
sigma2_hat <- SSE / (n-2)
# 这个通过残差估计出来的刹车距离误差开根号应该是实际刹车距离和预测值之间的偏差
sigma2_hat
sigma_hat <- sqrt(sigma2_hat)


# 验算手工计算的对不对
summary(fit)$sigma


# 标准误 用来衡量估计的参数beta的不确定性

# 计算beta_1 的标准误（Standard Error）
# SE(beta1_hat) = sigma_hat / sqrt(Sxx) ，Sxx = Sum((xi - x_bar)^2)
x_bar <- mean(x)
x_bar
Sxx <- sum((x - x_bar)^2)
Sxx
SE_beta1 <- sigma_hat / sqrt(Sxx)
SE_beta1

# 验证下
summary(fit)$coefficients[2, 2]

# 计算beta_0 的标准误 -- 根据公式，不用记忆
SE_beta0 <- sigma_hat * sqrt(1/n + x_bar^2/Sxx)
SE_beta0

# 验证
summary(fit)$coefficients[1, 2]


# 构造 beta1_hat 的 95% 的置信区间
# qt 函数是t分布的分位数函数，因为双侧95%置信区间，左右各留2.5%
t_crit <- qt(0.975, df = n - 2)
t_crit
# 计算置信区间
CI_lower <- beta1_hat - t_crit * SE_beta1
CI_upper <- beta1_hat + t_crit * SE_beta1
CI_lower
CI_upper
# 验证，confint() 是计算回归系数置信区间的函数
confint(fit)

# 计算t统计量
# t = （估计值 - 加设值）/ 标准误
# 下面假设 H0 ： beta1 = 0（速度与刹车距离无关）
t_stat <- beta1_hat / SE_beta1
t_stat
# 验证，数值很大，说明beta1显著不为0，速度与刹车距离确实有关
summary(fit)$coefficients[2, 3]

# 计算 p 值： p = 2 * P(T > |tobs|)
# tobs 是算出来的t统计量
# T 服从自由度 n - 2 的t分布
# t 分布是对称的，中心在 0
# pt(x, df) 返回的是 P(T < x)
# 我们要求的是两边尾巴概率，所以 * 2
p_value <- 2 * pt(-abs(t_stat),df = n - 2)
p_value
# 验证
summary(fit)$coefficients[2, 4]

# 计算R2，R2 模型解释了百分之多少的变异
# SST
Y_bar <- mean(Y)
Y_bar
SST <- sum((Y - Y_bar)^2)
SST
R2 <- 1 - SSE/SST
R2
# 验证
summary(fit)$r.squared


# 预测速度=25mph时的刹车距离
x_new <- 25
Y_pred <- beta0_hat + beta1_hat * x_new
Y_pred
# 计算预测的标准误
SE_pred <- sigma_hat * sqrt(1 + 1/n + (x_new - x_bar)^2/Sxx)
SE_pred
PI_lower <- Y_pred - t_crit * SE_pred
PI_upper <- Y_pred + t_crit * SE_pred
c(PI_lower, PI_upper)
# 验证
new_data <- data.frame(speed = 25)
predict(fit, new_data, interval = "prediction")

# 计算均值的标准误，比预测区间少了那个 1
SE_mean <- sigma_hat * sqrt(1/n + (x_new - x_bar)^2/Sxx)
SE_mean
# 构造均值的95%置信区间
CI_lower <- Y_pred - t_crit * SE_mean
CI_upper <- Y_pred + t_crit * SE_mean
c(CI_lower, CI_upper)
# 验证
predict(fit, new_data, interval = "confidence")





