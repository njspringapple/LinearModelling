set.seed(34)

# 造数据, 数据生成过程
n <- 100
size <- runif(n, 40, 150)
year <- 1950 + round(runif(n, 0, 70))
rooms <- round(runif(n, 1, 6))
# 真实模型：price = 50 + 3*size + 0.5*year + 2*rooms + 噪声
price <- 50 + 3*size + 0.5*year + 2*rooms + rnorm(n, 0, 15)

# price 完全随机，和 size/year/rooms 无关
# price <- rnorm(n, mean = 500, sd = 100)

dat <- data.frame(price,size,year,rooms)

# 对 year 进行中心化（减去均值）
#dat$year_c <- dat$year - mean(dat$year)

# 标准化所有自变量（减去均值，除以标准差）
dat$size <- scale(dat$size)
dat$year <- scale(dat$year)
dat$rooms <- scale(dat$rooms)

head(dat)

# ---------------------------------------------------------------------
# Step -1 手工构造X，Y，计算 β̂ = (X'X)^(-1) X'Y
# ---------------------------------------------------------------------
X <- cbind(1,dat$size, dat$year, dat$rooms)
Y <- dat$price
colnames(X) <- c("Intercept", "size", "year", "rooms")

n <- nrow(X)
p_ <- ncol(X)
p <- p_ - 1

XtX <- t(X) %*% X
XtX_inv <- solve(XtX)
beta_hat <- XtX_inv %*% t(X) %*% Y
beta_hat

# 和lm函数对比
fit <- lm(price ~ size+year+rooms,data=dat)
coef(fit)

# ---------------------------------------------------------------------
# Step -2 Hat矩阵，H = X(X'X)^(-1)X' 和 I-H
# 性质: H' = H, H^2 = H(幂等)，tr(H) = p‘， tr(I-H) = n - p'
# ---------------------------------------------------------------------
H <- X %*% XtX_inv %*% t(X)
I_n <- diag(n)
M <- I_n - H   # 残差投影矩阵，有时候记作 Q

cat("tr(H)     = ", sum(diag(H)),  " 应该 = ", p_, "\n")
cat("tr(I-H)   = ", sum(diag(M)),  " 应该 = ", n - p_, "\n")
cat("H是否幂等 :", all.equal(H %*% H, H), "\n")
cat("H是否对称 :", all.equal(H, t(H)), "\n")

Y_hat <- H %*% Y           # 拟合值是H矩阵（投影矩阵）变换原始Y（也成为PY）
e <- M %*% Y               # 残差值是（I-H）矩阵变换原始Y

# ---------------------------------------------------------------------
# Step -3  平方和分解 SST = SSM + SSE
# ---------------------------------------------------------------------
Ybar <- mean(Y)
SST <- sum((Y - Ybar)^2)      # 总误差
SSM <- sum((Y_hat - Ybar)^2)  # 总模型解释部分误差 
SSE <- sum(e^2)               # 总残差

cat("SST = ", SST, "\n")
cat("SSM = ", SSM, "\n")
cat("SSE = ", SSE, "\n")
cat("SSM + SSE = ", SSM + SSE, " (应 = SST)\n")

R2 <- 1 - SSE/SST
cat("R2 = ", R2, " lm里面的：", summary(fit)$r.squared, "\n")

# ---------------------------------------------------------------------
# Step 4  σ² 的无偏估计 σ² = SSE / (n - p')
# ---------------------------------------------------------------------
sigma2_hat <- SSE / (n - p_)
sigma_hat <- sqrt(sigma2_hat)
cat(" σ² = ", sigma_hat, " lm里面的：", summary(fit)$sigma,"\n")


# ---------------------------------------------------------------------
# Step 5  β̂ ~ N(β, σ²(X'X)^(-1))
# 标准误 se(β̂_k) = sqrt(σ̂² * [(X'X)^(-1)]_kk)
# ---------------------------------------------------------------------
Var_beta <- sigma2_hat * XtX_inv
se_beta <- sqrt(diag(Var_beta))
se_beta
print.default(summary(fit))
summary(fit)$coefficients[, "Std. Error"]   # 对比 ✓

# ---------------------------------------------------------------------
# Step 6 单个系数 t 检验 H0: β_k = 0
# T = (β̂_k - 0) / se(β̂_k) ~ t(n - p')
# ---------------------------------------------------------------------
t_stat <- beta_hat / se_beta
t_stat
#  pt - t 分布的累积分布函数（CDF）
#  t_stat 取绝对值，因为双侧检验只关心偏离0的程度，不关心方向
# lower.tail = TRUE 左尾概率，因为对称，所以乘以 2
p_value <- 2 * pt(abs(t_stat), df = n - p_, lower.tail = FALSE)
p_value
result_t <- data.frame(
  Estimate = round(beta_hat,4),
  Std.Err  = round(se_beta,4),
  t.value  = round(t_stat,4),
  p.value  = round(p_value,4)
)
result_t
summary(fit)$coefficients

# 以 rooms 系数的 t 统计量为例
t_val <- abs(t_stat[4])  # rooms 的 t 值
df <- n - p_

# 生成 t 分布曲线，t 分布的大部分概率密度都集中在 -4 到 4 之间
x <- seq(-5, 5, length.out = 300)
# dt 是 t 分布的密度函数
y <- dt(x, df = df)

# 画图
plot(x, y, type = "l", lwd = 2,
     main = paste0("双侧 t 检验 (df = ", df, ", |t| = ", round(t_val, 2), ")"),
     xlab = "t", ylab = "密度",
     col = "black")

# 填充左尾 (t < -|t|)
x_left <- seq(-5, -t_val, length.out = 100)
y_left <- dt(x_left, df = df)
polygon(c(x_left, rev(x_left)), c(y_left, rep(0, length(y_left))),
        col = rgb(1, 0, 0, 0.4), border = NA)

# 填充右尾 (t > |t|)
x_right <- seq(t_val, 5, length.out = 100)
y_right <- dt(x_right, df = df)
polygon(c(x_right, rev(x_right)), c(y_right, rep(0, length(y_right))),
        col = rgb(1, 0, 0, 0.4), border = NA)

# 标注 t 值位置
abline(v = c(-t_val, t_val), lty = 2, col = "red")

# 添加图例
p_val <- 2 * pt(t_val, df = df, lower.tail = FALSE)
legend("topright",
       legend = c(
         paste0("右尾: P(T > ", round(t_val, 2), ") = ", round(p_val/2, 4)),
         paste0("左尾: P(T < -", round(t_val, 2), ") = ", round(p_val/2, 4)),
         paste0("p-value = ", round(p_val, 4))
       ),
       fill = c(rgb(1, 0, 0, 0.4), rgb(1, 0, 0, 0.4), NA),
       border = c("black", "black", NA),
       bty = "n")

# p值大于0.05，红色区域较大，说明t值不够强，不能拒绝H0: β_rooms = 0（房间数对房价没有影响）

# ---------------------------------------------------------------------
# Step 7 整体F检验：H0：所有参数均 = 0，即无效
# F = (SSM/p) / (SSE/(n-p')) ~ F(p, n-p')
# ---------------------------------------------------------------------
MSM <- SSM / p
MSE <- SSE / (n - p_)
F_stat <- MSM / MSE
p_F    <- pf(F_stat, df1 = p, df2 = n - p_, lower.tail = FALSE)

cat("F 值 =", F_stat, "   p 值 =", p_F, "\n")
summary(fit)$fstatistic          # 对比


par(mfrow = c(1, 2))  # 一行两个图

# ========== 图1：完整范围，显示 F 统计量位置 ==========
x <- seq(0, F_stat + 50, length.out = 300)
y <- df(x, df1 = df1, df2 = df2)

plot(x, y, type = "l", lwd = 2,
     main = paste0("完整视图 (F = ", round(F_stat, 2), ")"),
     xlab = "F", ylab = "密度", col = "black")
abline(v = F_stat, lty = 2, col = "red", lwd = 2)
text(F_stat, max(y) * 0.5, labels = "F 值在这里\n(极端！)", col = "red", pos = 2)

# ========== 图2：放大 F 分布主体，展示拒绝域 ==========
F_crit <- qf(0.95, df1 = df1, df2 = df2)  # α=0.05 的临界值

x2 <- seq(0, 8, length.out = 300)
y2 <- df(x2, df1 = df1, df2 = df2)

plot(x2, y2, type = "l", lwd = 2,
     main = paste0("放大视图 (临界值 F_0.05 = ", round(F_crit, 2), ")"),
     xlab = "F", ylab = "密度", col = "black")

# 填充拒绝域 (F > F_crit)
x_rej <- seq(F_crit, 8, length.out = 100)
y_rej <- df(x_rej, df1 = df1, df2 = df2)
polygon(c(x_rej, 8, F_crit), c(y_rej, 0, 0),
        col = rgb(1, 0, 0, 0.4), border = NA)

abline(v = F_crit, lty = 2, col = "red", lwd = 2)
text(F_crit, df(F_crit, df1, df2) + 0.02, 
     labels = paste0("F_crit = ", round(F_crit, 2)), 
     col = "red", pos = 4)

# 动态生成结论
if (F_stat > F_crit) {
  comparison <- ">>"
  conclusion <- "结论：拒绝 H0，模型显著"
} else {
  comparison <- "<"
  conclusion <- "结论：不能拒绝 H0，模型不显著"
}

legend("topright",
       legend = c(
         paste0("红色区域 = 拒绝域 (α=0.05)"),
         paste0("你的 F = ", round(F_stat, 2), " ", comparison, " ", round(F_crit, 2)),
         conclusion
       ),
       bty = "n", cex = 0.8)

par(mfrow = c(1, 1))  # 恢复默认




# ---------------------------------------------------------------------
# 检验线性约束 的三种方法 ： Wald法 和 SSE 比较 和 似然比
# ---------------------------------------------------------------------




# ---------------------------------------------------------------------
# STEP 8:线性假设(Wald)  Aβ = c
#   例:检验 size 与 rooms 效应是否相等  β_size - β_rooms = 0
# ---------------------------------------------------------------------
# 对应 X 矩阵  (Intercept, size, year, rooms)
A <- matrix(c(0, 1, 0, -1), nrow = 1)    # 对应 (Intercept, size, year, rooms)
cc <- 0
a  <- nrow(A)      # 约束条件个数

Abeta_c   <- A %*% beta_hat - cc
mid       <- solve(A %*% XtX_inv %*% t(A))
SSH       <- t(Abeta_c) %*% mid %*% Abeta_c
MSH       <- SSH / a
F_wald    <- as.numeric( MSH / MSE )
p_wald    <- pf(F_wald, a, n - p_, lower.tail = FALSE)
cat("Wald F =", F_wald, "  p =", p_wald, "\n")

# 用 car::linearHypothesis 验证
# install.packages("car")
library(car)
linearHypothesis(fit, "size - rooms = 0")

# F 远远大于1， p 很小，拒绝假设，也就是size和rooms对 price 的影响系数不相等

# ---------------------------------------------------------------------
# STEP 9:条件约束(两模型 SSE 比较法)
#   SSH = SSE(H0) - SSE
#   F   = (SSH/q) / (SSE/(n-p')) ~ F(q, n-p')
# ---------------------------------------------------------------------
fit_restricted <- lm(price ~ I(size + rooms) + year, data = dat)  # 强制 β_size = β_rooms
SSE_H0 <- sum(resid(fit_restricted)^2)
SSH    <- SSE_H0 - SSE
MSH    <- SSH / a
F_res  <- MSH / MSE
cat("SSE(H0)=", SSE_H0, "  SSH=", SSH, "  F=", F_res, "\n")
# 和 STEP 8 的 F_wald 完全相同 ✓
anova(fit_restricted, fit)

# ---------------------------------------------------------------------
# STEP 10:似然比检验(LRT)
#   τ_LR = (SSE(H0)/SSE)^(-n/2)
#   -2 log τ 近似 ~ χ²(q)
# ---------------------------------------------------------------------
logLik_full <- logLik(fit)
logLik_res  <- logLik(fit_restricted)
LR_stat <- as.numeric(-2 * (logLik_res - logLik_full))
p_LR    <- pchisq(LR_stat, df = q, lower.tail = FALSE)
cat("LR 统计量 =", LR_stat, "  p =", p_LR, "\n")

# 用 lmtest::lrtest 验证
# install.packages("lmtest")
library(lmtest)
lrtest(fit_restricted, fit)

# ---------------------------------------------------------------------
# STEP 11:顺序平方和(Type I) vs 偏平方和(Type III)
# ---------------------------------------------------------------------
cat("\n--- Type I(顺序)---\n")
anova(lm(price ~ size + year + rooms, data=dat))   # 先 size 再 year 再 rooms
cat("\n--- 换顺序 ---\n")
anova(lm(price ~ year + size + rooms, data=dat))   # 顺序不同,第一项差很多!


cat("\n--- Type III(偏平方和)---\n")
Anova(fit, type = 3)                               # 每个变量都是在控制其他变量之后

# ---------------------------------------------------------------------
# STEP 12:偏杠杆图(Partial Leverage Plot / Added-Variable Plot)
#   "净化"size:把 size 和 price 都对 (year+rooms) 回归取残差,再画散点
# ---------------------------------------------------------------------
Y_star    <- resid(lm(price ~ year + rooms, data = dat))
size_star <- resid(lm(size  ~ year + rooms, data = dat))

plot(size_star, Y_star,
     main = "偏杠杆图:size 对 price(已控制 year, rooms)",
     xlab = "size* (残差)", ylab = "price* (残差)", pch = 19, col = "steelblue")
abline(lm(Y_star ~ size_star), col = "red", lwd = 2)

# 斜率应该 = β̂_size
cat("avplot 斜率 =", coef(lm(Y_star ~ size_star))[2],
    "  β̂_size =", beta_hat["size",], "\n")

# car 包内置函数
avPlots(fit)

# ---------------------------------------------------------------------
# STEP 13:预测区间 vs 置信区间
#   置信区间(均值):σ̂² · x'(X'X)^(-1)x
#   预测区间(单值):σ̂² · (1 + x'(X'X)^(-1)x)   ← 多了个 1
# ---------------------------------------------------------------------
# 获取标准化参数
size_center <- attr(dat$size, "scaled:center")
size_scale <- attr(dat$size, "scaled:scale")
year_center <- attr(dat$year, "scaled:center")
year_scale <- attr(dat$year, "scaled:scale")
rooms_center <- attr(dat$rooms, "scaled:center")
rooms_scale <- attr(dat$rooms, "scaled:scale")

# 标准化新数据
size_new <- (100 - size_center) / size_scale
year_new <- (2000 - year_center) / year_scale
rooms_new <- (4 - rooms_center) / rooms_scale

# 手动计算
x_new <- c(1, size_new, year_new, rooms_new)
y_new <- as.numeric(x_new %*% beta_hat)

var_conf <- sigma2_hat * t(x_new) %*% XtX_inv %*% x_new
var_pred <- sigma2_hat * (1 + t(x_new) %*% XtX_inv %*% x_new)

tcrit <- qt(0.975, df = n - p_)
ci <- y_new + c(-1,1) * tcrit * sqrt(var_conf)
pi <- y_new + c(-1,1) * tcrit * sqrt(var_pred)
cat("预测点估计:", y_new, "\n")
cat("置信区间(均值):", ci, "\n")
cat("预测区间(单值):", pi, "\n")

# 用 predict 验证（转成向量解决类型问题）
dat2 <- data.frame(
  price = dat$price,
  size = as.numeric(dat$size),
  year = as.numeric(dat$year),
  rooms = as.numeric(dat$rooms)
)
fit2 <- lm(price ~ size + year + rooms, data = dat2)

newdat <- data.frame(size = size_new, year = year_new, rooms = rooms_new)
predict(fit2, newdat, interval = "confidence")
predict(fit2, newdat, interval = "prediction")

# ---------------------------------------------------------------------
# STEP 15:渐近性(P23)—— n 大时用正态,n 小时用 t
# ---------------------------------------------------------------------
# 比较 t(n-p') 分位数 和 N(0,1) 分位数
cat("t_{0.975}(n-p') =", qt(0.975, n - p_), "\n")
cat("t_{0.975}(10*n-p') =", qt(s0.975, 10*n - p_), "\n")
cat("t_{0.975}(100*n-p') =", qt(0.975, 100*n - p_), "\n")
cat("t_{0.975}(1000*n-p') =", qt(0.975, 1000*n - p_), "\n")
cat("t_{0.975}(10000*n-p') =", qt(0.975, 10000*n - p_), "\n")
cat("z_{0.975}       =", qnorm(0.975),       "  <-- n 大时两者接近\n")












