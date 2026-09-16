# =====================================================
# 目的：演示三种残差的区别，学会找异常值
# 
# 三种残差：
# 1. 普通残差：直接看偏差，但不公平（位置不同方差不同）
# 2. 标准化残差：除以 σ̂√qᵢᵢ，公平比较不同位置的点
# 3. 学生化残差：删掉自己再算 σ̂，防止异常值自我掩护
#
# 最终目的：用学生化残差找异常值，|r*| > 2 或 3
# =====================================================


# =====================
# 1. 生成数据
# 目的：造一个有异常值的数据集
# =====================

set.seed(123)

n <- 20
x <- 1:n
y <- 2 + 0.5 * x + rnorm(n, sd = 1)  # 正常数据

y[18] <- y[18] + 10  # 故意加一个异常值

data <- data.frame(x = x, y = y)
print(data)


# =====================
# 2. 拟合模型
# 目的：得到回归模型，后面从中提取残差
# =====================

model <- lm(y ~ x)
summary(model)


# =====================
# 3. 提取各种残差
# 目的：对比三种残差的计算方式和结果
# =====================

# 普通残差：ε̂ᵢ = yᵢ - ŷᵢ
e_hat <- residuals(model)

# 帽子矩阵对角线 hᵢᵢ（杠杆值）
# hᵢᵢ 大 = 这个点对拟合线影响大
h <- hatvalues(model)

# qᵢᵢ = 1 - hᵢᵢ
# qᵢᵢ 小 = 残差方差小（边缘点）
# qᵢᵢ 大 = 残差方差大（中心点）
q <- 1 - h

# σ̂：用所有点估计的标准差
sigma_hat <- sigma(model)

# 标准化残差（手算）：rᵢ = ε̂ᵢ / (σ̂ × √qᵢᵢ)
# 目的：让不同位置的残差可以公平比较
r_manual <- e_hat / (sigma_hat * sqrt(q))

# 标准化残差（R内置）
r_builtin <- rstandard(model)

# 学生化残差（R内置）
# 目的：删掉第 i 个点再算 σ̂，防止异常值把 σ̂ 拉大
r_student <- rstudent(model)


# =====================
# 4. 整理结果
# 目的：一眼看清每个点的各种数值
# =====================

results <- data.frame(
  i = 1:n,
  x = x,
  y = round(y, 2),
  e_hat = round(e_hat, 2),      # 普通残差
  h_ii = round(h, 3),           # 杠杆值
  q_ii = round(q, 3),           # 1 - 杠杆值
  r_standard = round(r_builtin, 2),  # 标准化残差
  r_student = round(r_student, 2)    # 学生化残差
)

print(results)


# =====================
# 5. 可视化
# 目的：直观看出异常值和各种残差的区别
# =====================

par(mfrow = c(2, 2))

# 图1: 散点图
# 目的：看原始数据，异常值明显偏离拟合线
plot(x, y, pch = 19, main = "1. 散点图：异常值偏离拟合线")
abline(model, col = "blue", lwd = 2)
points(18, y[18], col = "red", pch = 19, cex = 2)
text(18, y[18], "异常值", pos = 2, col = "red")

# 图2: 普通残差
# 目的：能看出异常值，但无法判断"有多异常"
plot(x, e_hat, pch = 19, main = "2. 普通残差：能看出异常")
abline(h = 0, lty = 2)
points(18, e_hat[18], col = "red", pch = 19, cex = 2)
text(18, e_hat[18], "残差大", pos = 2, col = "red")

# 图3: 标准化 vs 学生化
# 目的：学生化残差更大，更能暴露异常值
plot(r_builtin, r_student, pch = 19,
     xlab = "标准化残差", ylab = "学生化残差",
     main = "3. 学生化残差更大")
abline(0, 1, lty = 2, col = "gray")
points(r_builtin[18], r_student[18], col = "red", pch = 19, cex = 2)
text(r_builtin[18], r_student[18], "异常值", pos = 4, col = "red")

# 图4: 学生化残差 + 阈值线
# 目的：用 |r*| > 2 或 3 判断异常值
plot(x, r_student, pch = 19, 
     main = "4. 用学生化残差找异常：|r*| > 2 或 3",
     ylab = "学生化残差 r*")
abline(h = 0)
abline(h = c(-2, 2), lty = 2, col = "orange")  # |r*| > 2：可疑
abline(h = c(-3, 3), lty = 2, col = "red")     # |r*| > 3：异常
points(18, r_student[18], col = "red", pch = 19, cex = 2)
legend("topleft", 
       legend = c("|r*| > 2：可疑", "|r*| > 3：异常"), 
       lty = 2, col = c("orange", "red"), cex = 0.8)


# =====================
# 6. 找异常值
# 目的：实际操作，用代码找出异常值是哪个
# =====================

cat("\n===== 找异常值 =====\n")
cat("|r*| > 2 的点:", which(abs(r_student) > 2), "\n")
cat("|r*| > 3 的点:", which(abs(r_student) > 3), "\n")


# =====================
# 7. 对比：为什么学生化残差更大
# 目的：验证异常值会把 σ̂ 拉大，掩盖自己
# =====================

cat("\n===== 为什么学生化残差更大 =====\n")
cat("第18个点:\n")
cat("  普通残差:", round(e_hat[18], 2), "\n")
cat("  标准化残差:", round(r_builtin[18], 2), "\n")
cat("  学生化残差:", round(r_student[18], 2), "← 更大！\n")

cat("\n用所有点算的 σ̂:", round(sigma_hat, 2), "\n")

# 删掉第18个点重新算
model_no18 <- lm(y ~ x, data = data[-18, ])
sigma_no18 <- sigma(model_no18)
cat("删掉异常值后的 σ̂:", round(sigma_no18, 2), "\n")
cat("异常值把 σ̂ 拉大了", round(sigma_hat / sigma_no18, 1), "倍!\n")

cat("\n结论：标准化残差用了被拉大的 σ̂，所以偏小\n")
cat("      学生化残差删掉自己再算 σ̂，所以更准确\n")


# =====================
# 8. hᵢᵢ 和 qᵢᵢ 的规律
# 目的：理解为什么要除以 √qᵢᵢ
# =====================

par(mfrow = c(1, 2))

# hᵢᵢ：杠杆值
# 两端大（边缘点影响大），中间小
plot(x, h, pch = 19, main = "hᵢᵢ：两端大，中间小",
     ylab = "hᵢᵢ (杠杆值)")

# qᵢᵢ = 1 - hᵢᵢ
# 两端小（边缘点残差方差小），中间大
plot(x, q, pch = 19, main = "qᵢᵢ：两端小，中间大",
     ylab = "qᵢᵢ (残差方差系数)")

cat("\n===== hᵢᵢ 和 qᵢᵢ 的含义 =====\n")
cat("边缘点（x=1,20）：hᵢᵢ大，qᵢᵢ小，残差方差小\n")
cat("中心点（x=10）：hᵢᵢ小，qᵢᵢ大，残差方差大\n")
cat("标准化时除以 √qᵢᵢ，让不同位置可以公平比较\n")


# =====================
# 9. 总结
# =====================

cat("\n===== 总结 =====\n")
cat("1. 普通残差：能看出异常，但不知道'有多异常'\n")
cat("2. 标准化残差：除以 σ̂√qᵢᵢ，公平比较不同位置\n")
cat("3. 学生化残差：删掉自己再算 σ̂，防止自我掩护\n")
cat("4. 找异常值：用 rstudent()，看 |r*| > 2 或 3\n")













# 模型1：只用 x
model1 <- lm(y ~ x)
e1 <- residuals(model1)
q1 <- 1 - hatvalues(model1)
PRESS1 <- sum((e1 / q1)^2)

# 模型2：用 x 和 x²
model2 <- lm(y ~ x + I(x^2))
e2 <- residuals(model2)
q2 <- 1 - hatvalues(model2)
PRESS2 <- sum((e2 / q2)^2)

cat("模型1 PRESS:", PRESS1, "\n")
cat("模型2 PRESS:", PRESS2, "\n")












set.seed(42)

# 生成一组"未知"数据（实际来自指数分布，但假装不知道）
data <- rexp(300, rate = 0.5)

par(mfrow = c(1, 3))

# 猜测1：是不是正态分布？
qqnorm(data, main = "对比正态分布", col = "steelblue", pch = 16)
qqline(data, col = "red", lwd = 2)

# 猜测2：是不是指数分布？
qqplot(qexp(ppoints(length(data)), rate = 0.5), data,
       main = "对比指数分布",
       xlab = "指数分布理论分位数", ylab = "样本分位数",
       col = "darkorange", pch = 16)
qqline(data, distribution = function(p) qexp(p, rate = 0.5), col = "red", lwd = 2)

# 猜测3：是不是均匀分布？
qqplot(qunif(ppoints(length(data)), min = 0, max = max(data)), data,
       main = "对比均匀分布",
       xlab = "均匀分布理论分位数", ylab = "样本分位数",
       col = "forestgreen", pch = 16)
qqline(data, distribution = function(p) qunif(p, min = 0, max = max(data)), col = "red", lwd = 2)












# ============================================================
# 完整示例：误差项非正态的诊断与解决
# 场景：Y 是计数数据（门店每日销量），真实服从泊松分布
# ============================================================

set.seed(42)
n <- 300

# ---- 1. 生成数据（假装我们不知道真实分布）----
X <- runif(n, 0, 3)                # 广告投入
lambda <- exp(0.2 + 0.3 * X)       # 真实条件均值
Y <- rpois(n, lambda)              # 每日销量（计数）

dat <- data.frame(Y = Y, X = X)

# ---- 2. 错误做法：直接用 OLS ----
model_ols <- lm(Y ~ X, data = dat)

# ---- 3. 诊断 ----

# 3a. 偏度和峰度（正态应分别接近 0 和 3）
library(e1071)
cat("偏度:", skewness(residuals(model_ols)), "\n")
cat("峰度:", kurtosis(residuals(model_ols)) + 3, "\n")
# 偏度远大于0 → 右偏；峰度远大于3 → 尖峰厚尾

# 3c. 画图诊断
par(mfrow = c(2, 2))

# 残差直方图
hist(residuals(model_ols), breaks = 30,
     main = "OLS 残差直方图", col = "lightblue", xlab = "残差")

# OLS 残差 QQ 图
qqnorm(residuals(model_ols),
       main = "OLS 残差 QQ图（对比正态）",
       pch = 16, col = "steelblue")
qqline(residuals(model_ols), col = "red", lwd = 2)
# 右上方明显上翘 → 右偏，误差不是正态

# ---- 4. 解决：用泊松 GLM ----
model_glm <- glm(Y ~ X, family = poisson, data = dat)

# GLM 偏差残差直方图
hist(residuals(model_glm, type = "deviance"), breaks = 30,
     main = "GLM 偏差残差直方图", col = "lightyellow", xlab = "偏差残差")

# GLM 偏差残差 QQ 图
qqnorm(residuals(model_glm, type = "deviance"),
       main = "GLM 偏差残差 QQ图（对比正态）",
       pch = 16, col = "darkorange")
qqline(residuals(model_glm, type = "deviance"), col = "red", lwd = 2)
# 点基本贴线 → 模型选对了

# ---- 5. 对比系数估计 ----
cat("\n===== OLS 结果 =====\n")
print(coef(summary(model_ols)))

cat("\n===== GLM 结果 =====\n")
print(coef(summary(model_glm)))
# GLM 估计的截距 ≈ 1，斜率 ≈ 0.5，接近真实值

# ---- 6. 对比预测区间（影响最大的地方）----
par(mfrow = c(1, 2))

newdata <- data.frame(X = seq(0, 3, length.out = 100))

# OLS 预测区间（基于正态假设，可能出负数！）
ols_pred <- predict(model_ols, newdata, interval = "prediction", level = 0.95)
plot(dat$X, dat$Y, pch = 16, col = "grey60",
     main = "OLS 预测区间", xlab = "X", ylab = "Y")
lines(newdata$X, ols_pred[, "fit"], col = "blue", lwd = 2)
lines(newdata$X, ols_pred[, "lwr"], col = "red", lty = 2, lwd = 2)
lines(newdata$X, ols_pred[, "upr"], col = "red", lty = 2, lwd = 2)
abline(h = 0, col = "black", lty = 3)
# 下界出现负数 → 计数数据不可能为负，说明模型不合适

# GLM 预测区间（用模拟方式生成）
glm_mu <- predict(model_glm, newdata, type = "response")
glm_lwr <- qpois(0.025, lambda = glm_mu)
glm_upr <- qpois(0.975, lambda = glm_mu)
plot(dat$X, dat$Y, pch = 16, col = "grey60",
     main = "GLM 预测区间", xlab = "X", ylab = "Y")
lines(newdata$X, glm_mu, col = "darkorange", lwd = 2)
lines(newdata$X, glm_lwr, col = "red", lty = 2, lwd = 2)
lines(newdata$X, glm_upr, col = "red", lty = 2, lwd = 2)
# 下界始终 ≥ 0，区间不对称，符合计数数据的特征








set.seed(42)
n <- 300

par(mfrow = c(1, 2))

# ============================================================
# 左图：适合 OLS 的数据（同方差）
# Y = 2 + 3X + 正态误差（方差恒定）
# ============================================================
X1 <- runif(n, 0, 5)
Y1 <- 2 + 3 * X1 + rnorm(n, 0, 2)   # 方差始终是 2²=4，不随 X 变化

model_good <- lm(Y1 ~ X1)

plot(fitted(model_good), residuals(model_good),
     pch = 16, col = "steelblue",
     xlab = "拟合值", ylab = "残差",
     main = "✅ 适合 OLS（等宽带状）")
abline(h = 0, col = "red", lwd = 2)

# ============================================================
# 右图：不适合 OLS 的数据（异方差）
# Y = 2 + 3X + 正态误差（方差随 X 增大）
# ============================================================
X2 <- runif(n, 0, 5)
Y2 <- 2 + 3 * X2 + rnorm(n, 0, 0.5 * X2)  # 标准差 = 0.5X，越大越散

model_bad <- lm(Y2 ~ X2)

plot(fitted(model_bad), residuals(model_bad),
     pch = 16, col = "darkorange",
     xlab = "拟合值", ylab = "残差",
     main = "❌ 不适合 OLS（喇叭形）")
abline(h = 0, col = "red", lwd = 2)














set.seed(42)
n <- 300

# 生成异方差数据（和刚才右图一样）
X <- runif(n, 0, 5)
Y <- 2 + 3 * X + rnorm(n, 0, 0.5 * X)   # 标准差 = 0.5X

model_bad <- lm(Y ~ X)

par(mfrow = c(2, 2))

# ---- 原始 OLS（有病）----
plot(fitted(model_bad), residuals(model_bad),
     pch = 16, col = "darkorange",
     xlab = "拟合值", ylab = "残差",
     main = "❌ 原始 OLS（喇叭形）")
abline(h = 0, col = "red", lwd = 2)

# ---- 治疗1：log 变换 ----
model_log <- lm(log(Y) ~ X, data = data.frame(X = X[Y > 0], Y = Y[Y > 0]))

plot(fitted(model_log), residuals(model_log),
     pch = 16, col = "steelblue",
     xlab = "拟合值", ylab = "残差",
     main = "治疗1：log(Y)")
abline(h = 0, col = "red", lwd = 2)

# ---- 治疗2：WLS（权重 = 1/X²）----
# 因为我们知道 Var ∝ X²，所以权重取 1/X²
wts <- 1 / X^2
model_wls <- lm(Y ~ X, weights = wts)

plot(fitted(model_wls), residuals(model_wls) * sqrt(wts),
     pch = 16, col = "forestgreen",
     xlab = "拟合值", ylab = "加权残差",
     main = "治疗2：WLS（权重=1/X²）")
abline(h = 0, col = "red", lwd = 2)

# ---- 治疗3：稳健标准误（模型不变，修正推断）----
library(sandwich)
library(lmtest)

# 残差图不变（模型本身没改），但系数检验变了
plot(fitted(model_bad), residuals(model_bad),
     pch = 16, col = "purple",
     xlab = "拟合值", ylab = "残差",
     main = "治疗3：稳健标准误\n（图不变，推断修正）")
abline(h = 0, col = "red", lwd = 2)

# ============================================================
# 对比系数估计和标准误
# ============================================================
cat("\n===== 原始 OLS =====\n")
print(coef(summary(model_bad)))

cat("\n===== WLS =====\n")
print(coef(summary(model_wls)))

cat("\n===== OLS + 稳健标准误 (HC3) =====\n")
print(coeftest(model_bad, vcov = vcovHC(model_bad, type = "HC3")))

# ============================================================
# BP 检验对比
# ============================================================
cat("\n===== BP 检验 =====\n")
cat("原始 OLS:  p =", bptest(model_bad)$p.value, "\n")
cat("log 变换:  p =", bptest(model_log)$p.value, "\n")
cat("WLS:       p =", bptest(model_wls)$p.value, "\n")


















# ============================================================
# 误差相关 vs 误差独立：逐个诊断方法对比
# ============================================================

set.seed(42)
n <- 200
X <- (1:n) / 50

# ---- 生成两组数据 ----

# 数据1：误差相关（AR(1), ρ=0.8）
eps_bad <- numeric(n)
eps_bad[1] <- rnorm(1)
for (i in 2:n) {
  eps_bad[i] <- 0.8 * eps_bad[i-1] + rnorm(1)
}
Y_bad <- 2 + 3 * X + eps_bad

# 数据2：误差独立
eps_good <- rnorm(n)
Y_good <- 2 + 3 * X + eps_good

model_bad  <- lm(Y_bad ~ X)
model_good <- lm(Y_good ~ X)

res_bad  <- residuals(model_bad)
res_good <- residuals(model_good)

# ============================================================
# 诊断方法1：残差 vs 时间
# ============================================================
par(mfrow = c(1, 2))

plot(1:n, res_bad, type = "l", col = "steelblue", lwd = 1.5,
     xlab = "时间", ylab = "残差",
     main = "❌ 有自相关：残差 vs 时间")
abline(h = 0, col = "red", lwd = 2)

plot(1:n, res_good, type = "l", col = "forestgreen", lwd = 1.5,
     xlab = "时间", ylab = "残差",
     main = "✅ 无自相关：残差 vs 时间")
abline(h = 0, col = "red", lwd = 2)

# ============================================================
# 诊断方法2：ε_t vs ε_{t-1}（滞后散点图）
# ============================================================
par(mfrow = c(1, 2))

plot(res_bad[-n], res_bad[-1],
     pch = 16, col = "darkorange",
     xlab = expression(hat(epsilon)[t-1]),
     ylab = expression(hat(epsilon)[t]),
     main = "❌ 有自相关：滞后散点图")
abline(lm(res_bad[-1] ~ res_bad[-n]), col = "red", lwd = 2)

plot(res_good[-n], res_good[-1],
     pch = 16, col = "darkorange",
     xlab = expression(hat(epsilon)[t-1]),
     ylab = expression(hat(epsilon)[t]),
     main = "✅ 无自相关：滞后散点图")
abline(lm(res_good[-1] ~ res_good[-n]), col = "red", lwd = 2)

# ============================================================
# 诊断方法3：ACF 图
# ============================================================
par(mfrow = c(1, 2))

acf(res_bad,  main = "❌ 有自相关：ACF")
acf(res_good, main = "✅ 无自相关：ACF")

# ============================================================
# Durbin-Watson 检验
# ============================================================
library(lmtest)
cat("\n===== Durbin-Watson 检验 =====\n")
cat("有自相关: DW =", dwtest(model_bad)$statistic,
    " p =", dwtest(model_bad)$p.value, "\n")
cat("无自相关: DW =", dwtest(model_good)$statistic,
    " p =", dwtest(model_good)$p.value, "\n")









# ============================================================
# Durbin-Watson 判定：两组数据分开画
# ============================================================

set.seed(42)
n <- 200
X <- (1:n) / 50

# ---- 数据1：有自相关（ρ=0.8）----
eps_bad <- numeric(n)
eps_bad[1] <- rnorm(1)
for (i in 2:n) {
  eps_bad[i] <- 0.8 * eps_bad[i-1] + rnorm(1)
}
Y_bad <- 2 + 3 * X + eps_bad
model_bad <- lm(Y_bad ~ X)

# ---- 数据2：无自相关 ----
Y_good <- 2 + 3 * X + rnorm(n)
model_good <- lm(Y_good ~ X)

# ---- 计算 DW 值 ----
res_bad <- residuals(model_bad)
res_good <- residuals(model_good)
dw_bad  <- sum(diff(res_bad)^2) / sum(res_bad^2)
dw_good <- sum(diff(res_good)^2) / sum(res_good^2)

# ============================================================
# 画图：左右分开
# ============================================================
par(mfrow = c(1, 2))

# ---- 左图：有自相关 ----
plot(NULL, xlim = c(0, 4), ylim = c(0, 1),
     xlab = "DW 统计量", ylab = "",
     main = "❌ 有自相关：DW 判定", yaxt = "n")

rect(0, 0, 0.8, 1, col = rgb(1, 0, 0, 0.2), border = NA)
rect(0.8, 0, 1.2, 1, col = rgb(1, 1, 0, 0.2), border = NA)
rect(1.2, 0, 2.8, 1, col = rgb(0, 1, 0, 0.2), border = NA)
rect(2.8, 0, 3.2, 1, col = rgb(1, 1, 0, 0.2), border = NA)
rect(3.2, 0, 4, 1, col = rgb(1, 0, 0, 0.2), border = NA)

text(0.4, 0.3, "正自相关\n拒绝H0", cex = 0.7)
text(1.0, 0.3, "?", cex = 1)
text(2.0, 0.3, "无自相关\n不拒绝H0", cex = 0.7)
text(3.0, 0.3, "?", cex = 1)
text(3.6, 0.3, "负自相关\n拒绝H0", cex = 0.7)

abline(v = dw_bad, col = "red", lwd = 3, lty = 2)
text(dw_bad + 0.3, 0.8, paste0("DW=", round(dw_bad, 2)),
     col = "red", cex = 1.2, font = 2)
arrows(dw_bad + 0.2, 0.7, dw_bad + 0.02, 0.6,
       col = "red", lwd = 2, length = 0.1)

# ---- 右图：无自相关 ----
plot(NULL, xlim = c(0, 4), ylim = c(0, 1),
     xlab = "DW 统计量", ylab = "",
     main = "✅ 无自相关：DW 判定", yaxt = "n")

rect(0, 0, 0.8, 1, col = rgb(1, 0, 0, 0.2), border = NA)
rect(0.8, 0, 1.2, 1, col = rgb(1, 1, 0, 0.2), border = NA)
rect(1.2, 0, 2.8, 1, col = rgb(0, 1, 0, 0.2), border = NA)
rect(2.8, 0, 3.2, 1, col = rgb(1, 1, 0, 0.2), border = NA)
rect(3.2, 0, 4, 1, col = rgb(1, 0, 0, 0.2), border = NA)

text(0.4, 0.3, "正自相关\n拒绝H0", cex = 0.7)
text(1.0, 0.3, "?", cex = 1)
text(2.0, 0.3, "无自相关\n不拒绝H0", cex = 0.7)
text(3.0, 0.3, "?", cex = 1)
text(3.6, 0.3, "负自相关\n拒绝H0", cex = 0.7)

abline(v = dw_good, col = "darkgreen", lwd = 3, lty = 2)
text(dw_good + 0.3, 0.8, paste0("DW=", round(dw_good, 2)),
     col = "darkgreen", cex = 1.2, font = 2)
arrows(dw_good + 0.2, 0.7, dw_good + 0.02, 0.6,
       col = "darkgreen", lwd = 2, length = 0.1)





















# ============================================================
# 杠杆值 (Leverage) 可视化解释
# ============================================================

library(ggplot2)
library(gridExtra)
library(MASS)

set.seed(42)

# ============================================================
# 图1: 杠杆值的几何意义 —— X空间中的距离
# ============================================================

# 生成正常数据
n <- 30
x_normal <- rnorm(n - 1, mean = 5, sd = 1.5)
y_normal <- 2 + 1.5 * x_normal + rnorm(n - 1, sd = 1)

# 添加一个高杠杆点（X空间远离中心）
x_all <- c(x_normal, 15)        # 高杠杆点 x=15，远离均值5
y_all <- c(y_normal, 2 + 1.5 * 15 + rnorm(1, sd = 1))  # 符合真实模式

# 计算杠杆值
X_mat <- cbind(1, x_all)
H <- X_mat %*% solve(t(X_mat) %*% X_mat) %*% t(X_mat)
h_ii <- diag(H)

# 判定标准
p_prime <- 2  # 截距 + 1个自变量
threshold <- 2 * p_prime / length(x_all)

df1 <- data.frame(
  x = x_all, 
  y = y_all, 
  leverage = h_ii,
  high_lev = ifelse(h_ii > threshold, "高杠杆点", "正常点")
)

p1 <- ggplot(df1, aes(x = x, y = y, color = high_lev, size = leverage)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = mean(x_all), linetype = "dashed", color = "gray50") +
  annotate("text", x = mean(x_all) + 0.3, y = max(y_all), 
           label = paste0("x̄ = ", round(mean(x_all), 1)), hjust = 0, size = 3.5) +
  scale_color_manual(values = c("高杠杆点" = "red", "正常点" = "steelblue")) +
  scale_size_continuous(range = c(2, 8)) +
  labs(title = "图1: 杠杆值 = X空间中距数据中心的距离",
       subtitle = paste0("判定标准: h_ii > 2p'/n = ", round(threshold, 3)),
       x = "X", y = "Y", color = "类型", size = "h_ii") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# 图2: h_ii 与距离中心的关系
# ============================================================

dist_from_center <- (x_all - mean(x_all))^2

df2 <- data.frame(
  dist = dist_from_center,
  leverage = h_ii,
  high_lev = ifelse(h_ii > threshold, "高杠杆点", "正常点")
)

p2 <- ggplot(df2, aes(x = dist, y = leverage, color = high_lev)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = p_prime / length(x_all), linetype = "dotted", color = "gray50") +
  annotate("text", x = max(dist_from_center) * 0.6, y = threshold + 0.02,
           label = paste0("阈值 2p'/n = ", round(threshold, 3)), color = "red", size = 3.5) +
  annotate("text", x = max(dist_from_center) * 0.6, y = p_prime/length(x_all) + 0.015,
           label = paste0("平均值 p'/n = ", round(p_prime/length(x_all), 3)), 
           color = "gray40", size = 3.5) +
  scale_color_manual(values = c("高杠杆点" = "red", "正常点" = "steelblue")) +
  labs(title = "图2: h_ii 随 (x_i - x̄)² 单调递增",
       subtitle = "简单回归中: h_ii = 1/n + (x_i - x̄)² / Σ(x_j - x̄)²",
       x = expression((x[i] - bar(x))^2), y = expression(h[ii]),
       color = "类型") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# 图3: 高杠杆点对回归线的"拉力"
# ============================================================

# 场景A: 无高杠杆点
x_base <- rnorm(25, mean = 5, sd = 1.5)
y_base <- 2 + 1.5 * x_base + rnorm(25, sd = 1)

# 场景B: 高杠杆点符合模式（good leverage）
x_good <- c(x_base, 14)
y_good <- c(y_base, 2 + 1.5 * 14)  # 在回归线上

# 场景C: 高杠杆点偏离模式（bad leverage = influential point）
x_bad <- c(x_base, 14)
y_bad <- c(y_base, 5)  # 严重偏离

# 拟合模型
fit_base <- lm(y_base ~ x_base)
fit_good <- lm(y_good ~ x_good)
fit_bad <- lm(y_bad ~ x_bad)

df3 <- data.frame(
  x = c(x_base, x_good, x_bad),
  y = c(y_base, y_good, y_bad),
  scenario = rep(c("A: 无高杠杆点", "B: 好的高杠杆点(符合模式)", "C: 坏的高杠杆点(偏离模式)"),
                 c(length(x_base), length(x_good), length(x_bad))),
  point_type = c(rep("normal", 25), 
                 rep("normal", 25), "high_lev",
                 rep("normal", 25), "high_lev")
)

p3 <- ggplot(df3, aes(x = x, y = y)) +
  geom_point(aes(color = point_type, size = point_type), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen", linewidth = 1) +
  facet_wrap(~scenario, ncol = 3) +
  scale_color_manual(values = c("normal" = "steelblue", "high_lev" = "red"), guide = "none") +
  scale_size_manual(values = c("normal" = 2, "high_lev" = 5), guide = "none") +
  labs(title = "图3: 高杠杆点对回归线的影响",
       subtitle = "注意C中回归线被高杠杆异常值严重拉偏！",
       x = "X", y = "Y") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 10))

# ============================================================
# 图4: h_ii 作为自身权重的理解
# Ŷ_i = Σ h_ij * Y_j，h_ii 是 Y_i 对 Ŷ_i 的贡献权重
# ============================================================

# 用场景C展示
X_bad_mat <- cbind(1, x_bad)
H_bad <- X_bad_mat %*% solve(t(X_bad_mat) %*% X_bad_mat) %*% t(X_bad_mat)
h_bad_diag <- diag(H_bad)

# 最后一个点（高杠杆点）的权重行
weights_last <- H_bad[nrow(H_bad), ]

df4 <- data.frame(
  index = 1:length(x_bad),
  weight = weights_last,
  type = c(rep("其他观测", 25), "高杠杆点自身")
)

p4 <- ggplot(df4, aes(x = index, y = weight, fill = type)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("其他观测" = "steelblue", "高杠杆点自身" = "red")) +
  annotate("text", x = 26, y = weights_last[26] + 0.03,
           label = paste0("h_ii = ", round(weights_last[26], 3)), 
           color = "red", fontface = "bold", size = 4) +
  labs(title = expression(paste("图4: ", hat(Y)[i], " = ", Sigma, " h"[ij], " Y"[j], 
                                "  中各 Y 的权重")),
       subtitle = "高杠杆点的拟合值几乎完全由自身决定 (h_ii ≈ 1 时回归线被迫穿过该点)",
       x = "观测编号 j", y = expression(h[ij]), fill = "") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# 组合输出
# ============================================================

cat("=" , rep("=", 60), "\n")
cat("  杠杆值 (Leverage) 诊断总结\n")
cat("=" , rep("=", 60), "\n\n")

cat("模型参数:\n")
cat(sprintf("  n = %d, p' = %d\n", length(x_all), p_prime))
cat(sprintf("  平均杠杆值 p'/n = %.4f\n", p_prime/length(x_all)))
cat(sprintf("  高杠杆阈值 2p'/n = %.4f\n\n", threshold))

cat("各观测点杠杆值:\n")
cat(sprintf("  正常点范围: [%.4f, %.4f]\n", min(h_ii[-length(h_ii)]), max(h_ii[-length(h_ii)])))
cat(sprintf("  高杠杆点 h_ii = %.4f  ⚠️ (> %.4f)\n\n", h_ii[length(h_ii)], threshold))

cat("场景C中高杠杆异常点:\n")
cat(sprintf("  h_ii = %.4f (该点拟合值中 %.1f%% 来自自身)\n", 
            h_bad_diag[26], h_bad_diag[26] * 100))
cat(sprintf("  原始回归: Y = %.2f + %.2f * X\n", coef(fit_base)[1], coef(fit_base)[2]))
cat(sprintf("  加入后:   Y = %.2f + %.2f * X  ← 严重偏移!\n", coef(fit_bad)[1], coef(fit_bad)[2]))

# 绘图
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
             top = "杠杆值 (Leverage) h_ii 可视化解释")



















# ============================================================
# Cook 距离 (Cook's Distance) 可视化解释
# ============================================================

library(ggplot2)
library(gridExtra)

set.seed(123)

# ============================================================
# 生成数据：包含不同类型的异常点
# ============================================================

n <- 40
x_base <- rnorm(n - 3, mean = 5, sd = 1.5)
y_base <- 2 + 1.5 * x_base + rnorm(n - 3, sd = 1.2)

# 添加三种特殊点：
# 点A: 高杠杆 + 大残差（最危险！）
# 点B: 高杠杆 + 小残差（好的杠杆点）
# 点C: 低杠杆 + 大残差（普通异常值）

x_all <- c(x_base, 13, 13.5, 5)
y_all <- c(y_base, 
           5,                    # A: 高杠杆 + 偏离模式 (真实值应约 21.5)
           2 + 1.5 * 13.5,      # B: 高杠杆 + 符合模式
           15)                   # C: 低杠杆 + 大残差

point_labels <- c(rep("普通点", n - 3), 
                  "A: 高杠杆+大残差", 
                  "B: 高杠杆+小残差", 
                  "C: 低杠杆+大残差")

# 拟合模型
fit <- lm(y_all ~ x_all)

# 计算诊断量
X_mat <- cbind(1, x_all)
H <- X_mat %*% solve(t(X_mat) %*% X_mat) %*% t(X_mat)
h_ii <- diag(H)
p_prime <- 2
sigma2 <- sum(residuals(fit)^2) / (length(x_all) - p_prime)

# 内学生化残差
e_i <- residuals(fit)
r_i <- e_i / (sqrt(sigma2) * sqrt(1 - h_ii))

# Cook 距离（使用分解公式）
D_i <- (r_i^2 / p_prime) * (h_ii / (1 - h_ii))

# 阈值
threshold_cook <- 4 / length(x_all)

# ============================================================
# 图1: 散点图 + 标注三种特殊点
# ============================================================

df1 <- data.frame(
  x = x_all, y = y_all, 
  label = point_labels,
  cook = D_i, leverage = h_ii, resid = abs(r_i)
)

# 用点大小表示Cook距离
p1 <- ggplot(df1, aes(x = x, y = y)) +
  geom_smooth(method = "lm", se = FALSE, color = "darkgreen", 
              linewidth = 0.8, linetype = "dashed") +
  geom_point(aes(color = label, size = cook), alpha = 0.8) +
  scale_color_manual(values = c("普通点" = "gray50", 
                                "A: 高杠杆+大残差" = "red",
                                "B: 高杠杆+小残差" = "blue",
                                "C: 低杠杆+大残差" = "orange")) +
  scale_size_continuous(range = c(2, 10)) +
  annotate("text", x = 13, y = 5 - 1.5, label = "A", color = "red", 
           fontface = "bold", size = 5) +
  annotate("text", x = 13.5, y = 2 + 1.5*13.5 + 1.5, label = "B", color = "blue", 
           fontface = "bold", size = 5) +
  annotate("text", x = 5, y = 15 + 1, label = "C", color = "orange", 
           fontface = "bold", size = 5) +
  labs(title = "图1: Cook距离 = 点大小（越大影响越大）",
       subtitle = "A最危险(高杠杆+大残差), B安全(高杠杆但符合模式), C有限影响(大残差但低杠杆)",
       x = "X", y = "Y", color = "点类型", size = expression(D[i])) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

# ============================================================
# 图2: Cook 距离的分解 D_i = (r²/p') × (h/(1-h))
# ============================================================

df2 <- data.frame(
  leverage_component = h_ii / (1 - h_ii),
  residual_component = r_i^2 / p_prime,
  cook = D_i,
  label = point_labels
)

# 绘制等高线（等Cook距离线）
lev_seq <- seq(0, max(df2$leverage_component) * 1.1, length.out = 200)
res_seq <- seq(0, max(df2$residual_component) * 1.1, length.out = 200)
grid <- expand.grid(lev = lev_seq, res = res_seq)
grid$cook <- grid$lev * grid$res

p2 <- ggplot(df2, aes(x = leverage_component, y = residual_component)) +
  # 等Cook距离线
  geom_contour(data = grid, aes(x = lev, y = res, z = cook),
               breaks = c(threshold_cook, 0.5, 1), 
               color = "gray70", linetype = "dashed") +
  geom_point(aes(color = label, size = cook), alpha = 0.8) +
  scale_color_manual(values = c("普通点" = "gray50", 
                                "A: 高杠杆+大残差" = "red",
                                "B: 高杠杆+小残差" = "blue",
                                "C: 低杠杆+大残差" = "orange")) +
  scale_size_continuous(range = c(2, 8)) +
  annotate("text", x = df2$leverage_component[n-2] + 0.1, 
           y = df2$residual_component[n-2] + 0.5,
           label = "A", color = "red", fontface = "bold", size = 5) +
  annotate("text", x = df2$leverage_component[n-1] + 0.1, 
           y = df2$residual_component[n-1] + 0.5,
           label = "B", color = "blue", fontface = "bold", size = 5) +
  annotate("text", x = df2$leverage_component[n] + 0.05, 
           y = df2$residual_component[n] + 0.5,
           label = "C", color = "orange", fontface = "bold", size = 5) +
  labs(title = expression(paste("图2: Cook距离分解  ", D[i], " = ", 
                                frac(r[i]^2, "p'"), " × ", frac(h[ii], "1-"*h[ii]))),
       subtitle = "虚线为等Cook距离线; 只有两个分量都大时Cook距离才大",
       x = expression(paste("杠杆分量: ", h[ii]/(1 - h[ii]))),
       y = expression(paste("残差分量: ", r[i]^2/p*"'")),
       color = "点类型") +
  guides(size = "none") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# 图3: 删除第i个观测后回归线的变化
# ============================================================

# 计算删除不同点后的回归线
fit_no_A <- lm(y_all[-( n-2)] ~ x_all[-(n-2)])
fit_no_B <- lm(y_all[-(n-1)] ~ x_all[-(n-1)])
fit_no_C <- lm(y_all[-n] ~ x_all[-n])

x_pred <- seq(min(x_all) - 1, max(x_all) + 1, length.out = 100)

df3 <- data.frame(
  x = rep(x_pred, 4),
  y = c(coef(fit)[1] + coef(fit)[2] * x_pred,
        coef(fit_no_A)[1] + coef(fit_no_A)[2] * x_pred,
        coef(fit_no_B)[1] + coef(fit_no_B)[2] * x_pred,
        coef(fit_no_C)[1] + coef(fit_no_C)[2] * x_pred),
  model = rep(c("全部数据", "删除A后", "删除B后", "删除C后"), each = 100)
)

p3 <- ggplot() +
  geom_point(data = df1, aes(x = x, y = y, color = label), alpha = 0.6, size = 2) +
  geom_line(data = df3, aes(x = x, y = y, linetype = model, color = model), linewidth = 0.9) +
  scale_color_manual(values = c("普通点" = "gray60",
                                "A: 高杠杆+大残差" = "red",
                                "B: 高杠杆+小残差" = "blue",
                                "C: 低杠杆+大残差" = "orange",
                                "全部数据" = "black",
                                "删除A后" = "red",
                                "删除B后" = "blue",
                                "删除C后" = "orange")) +
  scale_linetype_manual(values = c("全部数据" = "solid", 
                                   "删除A后" = "dashed",
                                   "删除B后" = "dotdash",
                                   "删除C后" = "dotted")) +
  labs(title = "图3: 删除单个观测后回归线的变化",
       subtitle = sprintf("删A后斜率变化 Δβ₁=%.2f; 删B后 Δβ₁=%.2f; 删C后 Δβ₁=%.2f",
                          coef(fit_no_A)[2] - coef(fit)[2],
                          coef(fit_no_B)[2] - coef(fit)[2],
                          coef(fit_no_C)[2] - coef(fit)[2]),
       x = "X", y = "Y", linetype = "回归线", color = "") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

# ============================================================
# 图4: Cook 距离柱状图（诊断图）
# ============================================================

df4 <- data.frame(
  index = 1:length(x_all),
  cook = D_i,
  label = point_labels,
  high = ifelse(D_i > threshold_cook, "超阈值", "正常")
)

p4 <- ggplot(df4, aes(x = index, y = cook, fill = label)) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_hline(yintercept = threshold_cook, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "darkred", linewidth = 0.6) +
  annotate("text", x = 3, y = threshold_cook + max(D_i) * 0.03,
           label = paste0("4/n = ", round(threshold_cook, 3)), color = "red", size = 3.5) +
  annotate("text", x = 3, y = 1 + max(D_i) * 0.03,
           label = "D = 1", color = "darkred", size = 3.5) +
  scale_fill_manual(values = c("普通点" = "gray60",
                               "A: 高杠杆+大残差" = "red",
                               "B: 高杠杆+小残差" = "blue",
                               "C: 低杠杆+大残差" = "orange")) +
  labs(title = "图4: Cook距离诊断图",
       subtitle = "寻找明显突出的观测——A远超其他点",
       x = "观测编号", y = expression(D[i]), fill = "点类型") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# 图5: 拟合值变化向量 ||Ŷ_(i) - Ŷ||
# ============================================================

yhat_full <- fitted(fit)

# 逐一删除计算拟合值变化
delta_yhat <- sapply(1:length(x_all), function(i) {
  fit_i <- lm(y_all[-i] ~ x_all[-i])
  yhat_i <- cbind(1, x_all) %*% coef(fit_i)  # 对所有点预测
  sum((yhat_i - yhat_full)^2)
})
# 标准化
delta_yhat_std <- delta_yhat / (p_prime * sigma2)

df5 <- data.frame(
  index = 1:length(x_all),
  delta = delta_yhat_std,
  label = point_labels
)

p5 <- ggplot(df5, aes(x = index, y = delta, fill = label)) +
  geom_col(alpha = 0.8, width = 0.7) +
  scale_fill_manual(values = c("普通点" = "gray60",
                               "A: 高杠杆+大残差" = "red",
                               "B: 高杠杆+小残差" = "blue",
                               "C: 低杠杆+大残差" = "orange")) +
  labs(title = expression(paste("图5: 等价形式  ", D[i], " = ", 
                                frac(paste("(", hat(Y)["(i)"] - hat(Y), ")'(", 
                                           hat(Y)["(i)"] - hat(Y), ")"), 
                                     paste("p'", hat(sigma)^2)))),
       subtitle = "Cook距离 = 删除第i点后，全部拟合值的标准化变化量",
       x = "观测编号", 
       y = expression(paste("拟合值变化: ", group("||", hat(Y)["(i)"] - hat(Y), "||")^2 / (p*"'"*hat(sigma)^2))),
       fill = "点类型") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# 数值总结
# ============================================================

cat("\n", paste(rep("=", 65), collapse = ""), "\n")
cat("  Cook 距离诊断总结\n")
cat(paste(rep("=", 65), collapse = ""), "\n\n")

cat("三个特殊点的诊断量:\n")
cat(sprintf("  %-25s  h_ii    |r_i|   D_i      判断\n", "点"))
cat(paste(rep("-", 70), collapse = ""), "\n")

special_idx <- c(n-2, n-1, n)
special_names <- c("A: 高杠杆+大残差", "B: 高杠杆+小残差", "C: 低杠杆+大残差")

for (k in 1:3) {
  idx <- special_idx[k]
  flag <- ifelse(D_i[idx] > threshold_cook, " ⚠️ 影响点!", "  正常")
  cat(sprintf("  %-25s  %.4f  %.3f   %.4f  %s\n", 
              special_names[k], h_ii[idx], abs(r_i[idx]), D_i[idx], flag))
}

cat(paste(rep("-", 70), collapse = ""), "\n")
cat(sprintf("  普通点 D_i 范围: [%.5f, %.5f]\n", 
            min(D_i[1:(n-3)]), max(D_i[1:(n-3)])))
cat(sprintf("  阈值 4/n = %.4f\n\n", threshold_cook))

cat("关键结论:\n")
cat("  • 点A: 高杠杆 × 大残差 → Cook距离最大 → 最具影响力 ⚠️\n")
cat("  • 点B: 高杠杆 × 小残差 → Cook距离小 → 虽有杠杆但无害\n")
cat("  • 点C: 低杠杆 × 大残差 → Cook距离中等 → 影响有限\n")
cat("  → D_i 大 ⟺ 残差大 AND 杠杆高（二者缺一不可）\n\n")

cat("回归系数变化:\n")
cat(sprintf("  全部数据:   Y = %.3f + %.3f X\n", coef(fit)[1], coef(fit)[2]))
cat(sprintf("  删除A后:    Y = %.3f + %.3f X  (Δβ₁ = %+.3f)\n", 
            coef(fit_no_A)[1], coef(fit_no_A)[2], coef(fit_no_A)[2]-coef(fit)[2]))
cat(sprintf("  删除B后:    Y = %.3f + %.3f X  (Δβ₁ = %+.3f)\n", 
            coef(fit_no_B)[1], coef(fit_no_B)[2], coef(fit_no_B)[2]-coef(fit)[2]))
cat(sprintf("  删除C后:    Y = %.3f + %.3f X  (Δβ₁ = %+.3f)\n", 
            coef(fit_no_C)[1], coef(fit_no_C)[2], coef(fit_no_C)[2]-coef(fit)[2]))

# ============================================================
# 输出图形
# ============================================================

# 上面三张
grid.arrange(p1, p2, p4, ncol = 2, nrow = 2, 
             layout_matrix = rbind(c(1, 2), c(3, 3)),
             top = "Cook 距离 (Cook's Distance) 可视化解释 — 第1页")

# 下面两张
grid.arrange(p3, p5, ncol = 1, nrow = 2,
             top = "Cook 距离 (Cook's Distance) 可视化解释 — 第2页")



















# ============================================================
# 模型设定错误 & Partial Leverage Plot（偏杠杆图）可视化
# ============================================================

library(ggplot2)
library(gridExtra)

set.seed(42)

# ============================================================
# 场景设计：y = β0 + β1*x1 + β2*x2 + β3*x2² + ε
# 真实模型中 x2 有非线性效应，但分析者只用了线性项
# ============================================================

n <- 120

# 两个预测变量（有一定相关性）
x1 <- rnorm(n, mean = 3, sd = 1.5)
x2 <- 0.4 * x1 + rnorm(n, mean = 0, sd = 1.2)  # x2 与 x1 有相关

# 真实模型：y 与 x2 有二次关系
beta0 <- 2
beta1 <- 1.5
beta2 <- -0.8
beta3 <- 0.6   # 二次项系数

y <- beta0 + beta1 * x1 + beta2 * x2 + beta3 * x2^2 + rnorm(n, sd = 1.0)

# ============================================================
# 模型对比
# ============================================================

# 错误模型（遗漏二次项）
fit_wrong <- lm(y ~ x1 + x2)

# 正确模型
fit_correct <- lm(y ~ x1 + x2 + I(x2^2))

# ============================================================
# 图1: 错误模型的残差 vs 拟合值（诊断）
# ============================================================

df_resid <- data.frame(
  fitted = fitted(fit_wrong),
  residuals = residuals(fit_wrong)
)

p1 <- ggplot(df_resid, aes(x = fitted, y = residuals)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.6, color = "steelblue", size = 2) +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 1) +
  labs(title = "图1: 错误模型的残差图 (y ~ x1 + x2)",
       subtitle = "残差呈现明显弯曲模式 → 提示遗漏了非线性项",
       x = expression(hat(y)[i]), 
       y = expression(hat(epsilon)[i])) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# 图2: 正确模型的残差 vs 拟合值
# ============================================================

df_resid2 <- data.frame(
  fitted = fitted(fit_correct),
  residuals = residuals(fit_correct)
)

p2 <- ggplot(df_resid2, aes(x = fitted, y = residuals)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.6, color = "forestgreen", size = 2) +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 1) +
  labs(title = "图2: 正确模型的残差图 (y ~ x1 + x2 + x2²)",
       subtitle = "残差随机分布，无系统模式 → 模型设定正确",
       x = expression(hat(y)[i]), 
       y = expression(hat(epsilon)[i])) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# 图3-5: Partial Leverage Plot 的构造过程
# ============================================================

# --- 步骤演示：构造 x2 的 Partial Leverage Plot ---

# Step 1: 用其他变量拟合 y（不含 x2）
fit_y_on_others <- lm(y ~ x1)
y_star <- residuals(fit_y_on_others)  # Q_(k) * y

# Step 2: 用其他变量拟合 x2（不含 x2 自身）
fit_x2_on_others <- lm(x2 ~ x1)
x2_star <- residuals(fit_x2_on_others)  # Q_(k) * x_k

# --- 图3: 原始 y vs x2 的散点图（被 x1 混淆）---

df_raw <- data.frame(x2 = x2, y = y)

p3 <- ggplot(df_raw, aes(x = x2, y = y)) +
  geom_point(alpha = 0.5, color = "gray40", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.8) +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 0.8, 
              linetype = "dashed") +
  labs(title = "图3: 原始散点图 y vs x2",
       subtitle = "被 x1 混淆，难以看清 x2 的真实非线性效应",
       x = expression(x[2]), y = "y") +
  annotate("text", x = min(x2) + 0.5, y = max(y) - 1, 
           label = "蓝=线性拟合, 红=LOESS", hjust = 0, size = 3.5) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# --- 图4: 构造过程示意 ---

df_step <- data.frame(
  y_star = y_star,
  x2_star = x2_star,
  x1_val = x1
)

# 验证：x2_star 对 y_star 的回归斜率 = 完整模型中 x2 的系数
fit_partial <- lm(y_star ~ x2_star)

p4 <- ggplot(df_step, aes(x = x2_star, y = y_star)) +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_hline(yintercept = 0, color = "gray80") +
  geom_point(alpha = 0.6, color = "darkorchid", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.8) +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 1.2,
              linetype = "dashed") +
  labs(title = expression(paste("图4: Partial Leverage Plot — ", x[2], "的偏杠杆图")),
       subtitle = "去除x1影响后, y*与x2*的关系 — 明显看到二次弯曲!",
       x = expression(paste(x[2]^"*", " = ", Q["(k)"] * x[2], "  (去除x1影响后的x2)")),
       y = expression(paste(y^"*", " = ", Q["(k)"] * y, "  (去除x1影响后的y)"))) +
  annotate("text", x = min(x2_star) + 0.3, y = max(y_star) - 0.5,
           label = "红色LOESS曲线揭示\n非线性关系!", 
           color = "red", size = 3.5, hjust = 0, fontface = "bold") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# 图5: 对比 — x1 的 Partial Leverage Plot（线性关系）
# ============================================================

# x1 的偏杠杆图
fit_y_no_x1 <- lm(y ~ x2)
y_star_x1 <- residuals(fit_y_no_x1)

fit_x1_no_x1 <- lm(x1 ~ x2)
x1_star <- residuals(fit_x1_no_x1)

df_partial_x1 <- data.frame(x1_star = x1_star, y_star = y_star_x1)

p5 <- ggplot(df_partial_x1, aes(x = x1_star, y = y_star)) +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_hline(yintercept = 0, color = "gray80") +
  geom_point(alpha = 0.6, color = "steelblue", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.8) +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 1.2,
              linetype = "dashed") +
  labs(title = expression(paste("图5: Partial Leverage Plot — ", x[1], "的偏杠杆图")),
       subtitle = "去除x2影响后, y*与x1*的关系 — 线性关系良好, 无需变换",
       x = expression(paste(x[1]^"*", " = ", Q["(k)"] * x[1], "  (去除x2影响后的x1)")),
       y = expression(paste(y^"*", " = ", Q["(k)"] * y, "  (去除x2影响后的y)"))) +
  annotate("text", x = min(x1_star) + 0.3, y = max(y_star_x1) - 0.5,
           label = "LOESS几乎与线性重合\n→ x1的线性假设正确", 
           color = "blue", size = 3.5, hjust = 0) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# 图6: 修正后 — 加入 x2² 后的 Partial Leverage Plot
# ============================================================

# 在正确模型下，x2的偏杠杆图（控制 x1 和 x2²）
fit_y_correct_no_x2 <- lm(y ~ x1 + I(x2^2))
y_star_corrected <- residuals(fit_y_correct_no_x2)

fit_x2_correct_no_x2 <- lm(x2 ~ x1 + I(x2^2))
x2_star_corrected <- residuals(fit_x2_correct_no_x2)

df_corrected <- data.frame(x2_star = x2_star_corrected, y_star = y_star_corrected)

p6 <- ggplot(df_corrected, aes(x = x2_star, y = y_star)) +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_hline(yintercept = 0, color = "gray80") +
  geom_point(alpha = 0.6, color = "forestgreen", size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.8) +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 1.2,
              linetype = "dashed") +
  labs(title = expression(paste("图6: 修正后的偏杠杆图 — 控制", x[1], "和", x[2]^2, "后")),
       subtitle = "加入二次项后，x2的线性剩余关系良好 → 模型设定正确",
       x = expression(paste(x[2]^"*", "  (去除x1和x2²影响后)")),
       y = expression(paste(y^"*", "  (去除x1和x2²影响后)"))) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# ============================================================
# 输出
# ============================================================

grid.arrange(p1, p2, ncol = 2,
             top = "模型设定错误诊断 — 残差图对比")

grid.arrange(p3, p4, ncol = 2,
             top = "Partial Leverage Plot 构造过程")

grid.arrange(p5, p6, ncol = 2,
             top = "Partial Leverage Plot — 正确设定 vs 需要变换")

# ============================================================
# 数值总结
# ============================================================

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("  模型设定错误 & Partial Leverage Plot 诊断总结\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("真实数据生成过程 (DGP):\n")
cat(sprintf("  y = %.1f + %.1f·x1 + (%.1f)·x2 + %.1f·x2² + ε\n\n", 
            beta0, beta1, beta2, beta3))

cat("错误模型 (遗漏 x2²):\n")
cat(sprintf("  y = %.3f + %.3f·x1 + %.3f·x2\n", 
            coef(fit_wrong)[1], coef(fit_wrong)[2], coef(fit_wrong)[3]))
cat(sprintf("  R² = %.4f,  RSE = %.4f\n", 
            summary(fit_wrong)$r.squared, summary(fit_wrong)$sigma))
cat(sprintf("  β1 偏差: %.3f (真值 %.1f, 估计 %.3f)\n",
            coef(fit_wrong)[2] - beta1, beta1, coef(fit_wrong)[2]))
cat(sprintf("  β2 偏差: %.3f (真值 %.1f, 估计 %.3f)  ← 严重偏误!\n\n",
            coef(fit_wrong)[3] - beta2, beta2, coef(fit_wrong)[3]))

cat("正确模型:\n")
cat(sprintf("  y = %.3f + %.3f·x1 + %.3f·x2 + %.3f·x2²\n", 
            coef(fit_correct)[1], coef(fit_correct)[2], 
            coef(fit_correct)[3], coef(fit_correct)[4]))
cat(sprintf("  R² = %.4f,  RSE = %.4f\n", 
            summary(fit_correct)$r.squared, summary(fit_correct)$sigma))
cat(sprintf("  β1 偏差: %.3f (真值 %.1f, 估计 %.3f)\n",
            coef(fit_correct)[2] - beta1, beta1, coef(fit_correct)[2]))
cat(sprintf("  β2 偏差: %.3f (真值 %.1f, 估计 %.3f)\n",
            coef(fit_correct)[3] - beta2, beta2, coef(fit_correct)[3]))
cat(sprintf("  β3 偏差: %.3f (真值 %.1f, 估计 %.3f)\n\n",
            coef(fit_correct)[4] - beta3, beta3, coef(fit_correct)[4]))

# F-检验
cat("F-检验 (错误模型 vs 正确模型):\n")
anova_test <- anova(fit_wrong, fit_correct)
cat(sprintf("  F = %.2f,  p-value = %.2e  → 强烈拒绝遗漏二次项的模型\n\n",
            anova_test$F[2], anova_test$`Pr(>F)`[2]))

cat("Partial Leverage Plot 关键性质:\n")
cat(sprintf("  x2偏杠杆图中线性拟合斜率 = %.4f\n", coef(fit_partial)[2]))
cat(sprintf("  完整模型中 x2 系数       = %.4f  (二者相等!)\n", coef(fit_wrong)[3]))
cat("  → 偏杠杆图中的回归斜率 = 完整模型中对应变量的回归系数\n\n")

cat(paste(rep("-", 70), collapse = ""), "\n")
cat("Partial Leverage Plot 构造步骤:\n")
cat("  ① 用其他变量(x1)拟合 y → 残差 y* = Q_(k)·y\n")
cat("  ② 用其他变量(x1)拟合 xk → 残差 xk* = Q_(k)·xk\n")
cat("  ③ 画 y* vs xk* 散点图\n")
cat("  ④ 若散点呈非线性 → xk 需要变换或加入高阶项\n")
cat("  ⑤ 若散点呈线性 → xk 的线性假设成立\n")
cat(paste(rep("-", 70), collapse = ""), "\n")















# ============================================================
# Partial Leverage Plot（偏杠杆图）— 按定义步骤构造
# ============================================================

library(ggplot2)
library(gridExtra)

set.seed(42)
n <- 120

# 数据生成：y = 2 + 1.5*x1 - 0.8*x2 + 0.6*x2² + ε
x1 <- rnorm(n, mean = 3, sd = 1.5)
x2 <- 0.4 * x1 + rnorm(n, mean = 0, sd = 1.2)
y <- 2 + 1.5 * x1 - 0.8 * x2 + 0.6 * x2^2 + rnorm(n, sd = 1.0)

# ============================================================
# x2 的偏杠杆图（只用线性模型，控制 x1）
# ============================================================

# 步骤1: Q_(k) 为不包含 x2 的模型的残差矩阵
# 步骤2: y* = Q_(k) * y （从 y 中去除 x1 的影响）
y_star <- residuals(lm(y ~ x1))

# 步骤3: x2* = Q_(k) * x2 （从 x2 中去除 x1 的影响）
x2_star <- residuals(lm(x2 ~ x1))

# 步骤4: 画 y* 对 x2* 的散点图
df1 <- data.frame(x_star = x2_star, y_star = y_star)

p1 <- ggplot(df1, aes(x = x_star, y = y_star)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.8) +
  labs(title = expression(paste("Partial Leverage Plot: ", x[2])),
       subtitle = "控制 x1 后，y* 与 x2* 的关系 — 明显非线性，需要变换",
       x = expression(x[2]^"*"),
       y = expression(y^"*")) +
  theme_minimal(base_size = 12)

# ============================================================
# x1 的偏杠杆图（控制 x2）
# ============================================================

y_star_x1 <- residuals(lm(y ~ x2))
x1_star <- residuals(lm(x1 ~ x2))

df2 <- data.frame(x_star = x1_star, y_star = y_star_x1)

p2 <- ggplot(df2, aes(x = x_star, y = y_star)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.8) +
  labs(title = expression(paste("Partial Leverage Plot: ", x[1])),
       subtitle = "控制 x2 后，y* 与 x1* 的关系 — 线性良好，无需变换",
       x = expression(x[1]^"*"),
       y = expression(y^"*")) +
  theme_minimal(base_size = 12)

grid.arrange(p1, p2, ncol = 2)

















# ============================================================
# 共线性诊断：条件数 & VIF
# ============================================================

library(ggplot2)
library(gridExtra)

set.seed(42)
n <- 100

# ============================================================
# 构造不同程度共线性的数据
# ============================================================

x1 <- rnorm(n)

# 三种情况
x2_low <- rnorm(n)                          # 与 x1 几乎无关
x2_mid <- 0.9 * x1 + rnorm(n, sd = 0.5)    # 中等共线性
x2_high <- 0.99 * x1 + rnorm(n, sd = 0.1)  # 严重共线性

y_low <- 1 + 2 * x1 + 3 * x2_low + rnorm(n)
y_mid <- 1 + 2 * x1 + 3 * x2_mid + rnorm(n)
y_high <- 1 + 2 * x1 + 3 * x2_high + rnorm(n)

# ============================================================
# 计算条件数
# ============================================================

cond_number <- function(X) {
  XtX <- t(X) %*% X
  ev <- eigen(XtX)$values
  sqrt(max(ev) / min(ev))
}

X_low <- cbind(1, x1, x2_low)
X_mid <- cbind(1, x1, x2_mid)
X_high <- cbind(1, x1, x2_high)

K_low <- cond_number(X_low)
K_mid <- cond_number(X_mid)
K_high <- cond_number(X_high)

# ============================================================
# 计算 VIF
# ============================================================

vif_x1 <- function(x1, x2) {
  R2 <- summary(lm(x1 ~ x2))$r.squared
  1 / (1 - R2)
}

VIF_low <- vif_x1(x1, x2_low)
VIF_mid <- vif_x1(x1, x2_mid)
VIF_high <- vif_x1(x1, x2_high)

# ============================================================
# 展示后果：β 的置信区间宽度
# ============================================================

fit_low <- lm(y_low ~ x1 + x2_low)
fit_mid <- lm(y_mid ~ x1 + x2_mid)
fit_high <- lm(y_high ~ x1 + x2_high)

ci_low <- confint(fit_low)
ci_mid <- confint(fit_mid)
ci_high <- confint(fit_high)

# 画置信区间对比图
df_ci <- data.frame(
  scenario = rep(c("低共线性", "中等共线性", "严重共线性"), each = 2),
  variable = rep(c("β1 (x1)", "β2 (x2)"), 3),
  estimate = c(coef(fit_low)[2:3], coef(fit_mid)[2:3], coef(fit_high)[2:3]),
  lower = c(ci_low[2:3,1], ci_mid[2:3,1], ci_high[2:3,1]),
  upper = c(ci_low[2:3,2], ci_mid[2:3,2], ci_high[2:3,2])
)
df_ci$scenario <- factor(df_ci$scenario, levels = c("低共线性", "中等共线性", "严重共线性"))

p1 <- ggplot(df_ci, aes(x = variable, y = estimate, color = scenario)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), 
                  position = position_dodge(width = 0.5), size = 0.8) +
  geom_hline(yintercept = c(2, 3), linetype = "dashed", color = "gray50") +
  labs(title = "共线性对参数估计的影响",
       subtitle = "虚线 = 真值 (β1=2, β2=3)，线段 = 95% 置信区间",
       x = "", y = "参数估计值", color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p1)

# ============================================================
# 数值总结
# ============================================================

cat("\n条件数 K(X):\n")
cat(sprintf("  低共线性:   K = %.1f\n", K_low))
cat(sprintf("  中等共线性: K = %.1f\n", K_mid))
cat(sprintf("  严重共线性: K = %.1f\n", K_high))

cat("\nVIF (x1):\n")
cat(sprintf("  低共线性:   VIF = %.2f\n", VIF_low))
cat(sprintf("  中等共线性: VIF = %.2f\n", VIF_mid))
cat(sprintf("  严重共线性: VIF = %.2f\n", VIF_high))

cat("\nβ1 的标准误:\n")
cat(sprintf("  低共线性:   SE = %.3f\n", summary(fit_low)$coefficients[2,2]))
cat(sprintf("  中等共线性: SE = %.3f\n", summary(fit_mid)$coefficients[2,2]))
cat(sprintf("  严重共线性: SE = %.3f\n", summary(fit_high)$coefficients[2,2]))






















# ============================================================
# 测量误差 & 向零偏误（Attenuation Bias）演示
# ============================================================

library(ggplot2)
library(gridExtra)

set.seed(42)
n <- 200

# ============================================================
# 真实模型：y = 1 + 3*x + ε
# ============================================================

beta0 <- 1
beta1 <- 3

x_true <- rnorm(n, mean = 2, sd = 1)         # 真实的 x
y <- beta0 + beta1 * x_true + rnorm(n, sd = 1)  # y 无测量误差

# x 的测量误差：x_obs = x_true + 测量噪声
x_small  <- x_true + rnorm(n, sd = 0.3)   # 小误差
x_mid    <- x_true + rnorm(n, sd = 1.0)   # 中等误差
x_large  <- x_true + rnorm(n, sd = 2.0)   # 大误差

# ============================================================
# 回归对比
# ============================================================

fit_true  <- lm(y ~ x_true)
fit_small <- lm(y ~ x_small)
fit_mid   <- lm(y ~ x_mid)
fit_large <- lm(y ~ x_large)

# ============================================================
# 图：散点图 + 回归线对比
# ============================================================

df <- data.frame(
  x = c(x_true, x_small, x_mid, x_large),
  y = rep(y, 4),
  group = rep(c("无误差 (真实x)", "小误差 (σ=0.3)", 
                "中等误差 (σ=1.0)", "大误差 (σ=2.0)"), each = n)
)
df$group <- factor(df$group, levels = c("无误差 (真实x)", "小误差 (σ=0.3)", 
                                        "中等误差 (σ=1.0)", "大误差 (σ=2.0)"))

p1 <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  geom_abline(intercept = beta0, slope = beta1, 
              linetype = "dashed", color = "blue", linewidth = 0.8) +
  facet_wrap(~group, ncol = 2, scales = "free_x") +
  labs(title = "测量误差导致向零偏误 (Attenuation Bias)",
       subtitle = "蓝色虚线 = 真实关系 (β=3)，红色实线 = 估计的回归线",
       x = "观测到的 x", y = "y") +
  theme_minimal(base_size = 11)

print(p1)

# ============================================================
# 理论公式验证
# ============================================================

# 向零偏误公式：E[β̂] = β * σ²_x / (σ²_x + σ²_u)
# 其中 σ²_x = 真实 x 的方差, σ²_u = 测量误差方差

sigma2_x <- var(x_true)
theoretical <- function(sigma_u) beta1 * sigma2_x / (sigma2_x + sigma_u^2)

cat("\n向零偏误演示:\n")
cat(sprintf("  真实 β1 = %.1f\n\n", beta1))
cat(sprintf("  %-20s  估计β1    理论E[β̂]   偏误%%\n", "测量误差"))
cat(sprintf("  %-20s  %.3f    %.3f      %.1f%%\n", "无误差", 
            coef(fit_true)[2], beta1, 0))
cat(sprintf("  %-20s  %.3f    %.3f      %.1f%%\n", "小 (σ_u=0.3)",
            coef(fit_small)[2], theoretical(0.3),
            (1 - coef(fit_small)[2]/beta1)*100))
cat(sprintf("  %-20s  %.3f    %.3f      %.1f%%\n", "中等 (σ_u=1.0)",
            coef(fit_mid)[2], theoretical(1.0),
            (1 - coef(fit_mid)[2]/beta1)*100))
cat(sprintf("  %-20s  %.3f    %.3f      %.1f%%\n", "大 (σ_u=2.0)",
            coef(fit_large)[2], theoretical(2.0),
            (1 - coef(fit_large)[2]/beta1)*100))

cat("\n公式: E[β̂] = β · σ²_x / (σ²_x + σ²_u)\n")
cat("测量误差越大 → σ²_u 越大 → β̂ 越接近 0\n")











library(nlme)

# ============================================
# 模拟函数
# ============================================

simulate_once <- function() {
  rho <- 0.7
  n <- 50
  beta0_true <- 2
  beta1_true <- 0.5
  
  epsilon <- arima.sim(model = list(ar = rho), n = n)
  y <- beta0_true + beta1_true * (1:n) + epsilon
  data <- data.frame(x = 1:n, y = y)
  
  fit_ml <- gls(y ~ x, correlation = corAR1(), data = data, method = "ML")
  fit_reml <- gls(y ~ x, correlation = corAR1(), data = data, method = "REML")
  
  # 提取 ρ
  rho_ml <- coef(fit_ml$modelStruct$corStruct, unconstrained = FALSE)
  rho_reml <- coef(fit_reml$modelStruct$corStruct, unconstrained = FALSE)
  
  # 提取 β
  beta0_ml <- coef(fit_ml)[1]
  beta1_ml <- coef(fit_ml)[2]
  beta0_reml <- coef(fit_reml)[1]
  beta1_reml <- coef(fit_reml)[2]
  
  c(
    rho_ml = rho_ml, 
    rho_reml = rho_reml,
    beta0_ml = beta0_ml,
    beta1_ml = beta1_ml,
    beta0_reml = beta0_reml,
    beta1_reml = beta1_reml
  )
}

set.seed(456)
results <- t(replicate(500, simulate_once()))  # 模拟 500 次
# 看平均值
colMeans(results)
par(mfrow = c(2, 2))  # 2x2 四张图
# 图1：ρ 的比较
boxplot(
  results[, "rho_ml.Phi"], results[, "rho_reml.Phi"],
  names = c("ML", "REML"),
  main = "ρ 的估计",
  ylab = expression(hat(rho)),
  col = c("lightblue", "lightgreen")
)
abline(h = 0.7, col = "red", lty = 2, lwd = 2)
legend("bottomright", "真实值 0.7", col = "red", lty = 2)

# 图2：β0 的比较
boxplot(
  results[, "beta0_ml.(Intercept)"], results[, "beta0_reml.(Intercept)"],
  names = c("ML", "REML"),
  main = "β₀ (截距) 的估计",
  ylab = expression(hat(beta)[0]),
  col = c("lightblue", "lightgreen")
)
abline(h = 2, col = "red", lty = 2, lwd = 2)
legend("bottomright", "真实值 2", col = "red", lty = 2)

# 图3：β1 的比较
boxplot(
  results[, "beta1_ml.x"], results[, "beta1_reml.x"],
  names = c("ML", "REML"),
  main = "β₁ (斜率) 的估计",
  ylab = expression(hat(beta)[1]),
  col = c("lightblue", "lightgreen")
)
abline(h = 0.5, col = "red", lty = 2, lwd = 2)
legend("bottomright", "真实值 0.5", col = "red", lty = 2)

# 图4：汇总表格
plot.new()
text(0.5, 0.9, "估计结果汇总", cex = 1.5, font = 2)
text(0.5, 0.7, paste("ρ:  ML =", round(mean(results[, 1]), 3), 
                     " REML =", round(mean(results[, 2]), 3),
                     " 真实 = 0.7"), cex = 1.1)
text(0.5, 0.5, paste("β₀: ML =", round(mean(results[, 3]), 3), 
                     " REML =", round(mean(results[, 5]), 3),
                     " 真实 = 2"), cex = 1.1)
text(0.5, 0.3, paste("β₁: ML =", round(mean(results[, 4]), 3), 
                     " REML =", round(mean(results[, 6]), 3),
                     " 真实 = 0.5"), cex = 1.1)

par(mfrow = c(1, 1))  # 恢复默认

# ============================================
# 计算偏差
# ============================================

cat("\n===== 偏差比较 =====\n")
cat("ρ:  ML 偏差 =", round(mean(results[, 1]) - 0.7, 4), 
    " REML 偏差 =", round(mean(results[, 2]) - 0.7, 4), "\n")
cat("β₀: ML 偏差 =", round(mean(results[, 3]) - 2, 4), 
    " REML 偏差 =", round(mean(results[, 5]) - 2, 4), "\n")
cat("β₁: ML 偏差 =", round(mean(results[, 4]) - 0.5, 4), 
    " REML 偏差 =", round(mean(results[, 6]) - 0.5, 4), "\n")

