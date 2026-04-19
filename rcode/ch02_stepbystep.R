data(cars)
head(cars)
fit <- lm(dist~speed,data=cars)
summary(fit)
# 提取基本信息
n <- nrow(cars)
n

# 简单线性回归模型：
#   Y_i = beta0 + beta1 * x_i + epsilon_i
# 假设：
#   (1) 线性性：E(Y|x) = beta0 + beta1 * x
#   (2) 独立性：epsilon_i 相互独立
#   (3) 同方差：Var(epsilon_i) = sigma^2 （与 x 无关）
#   (4) 正态性：epsilon_i ~ N(0, sigma^2)  —— 做 t 检验/CI 时需要



# 先观察数据：是否存在线性关系？
plot(cars$speed, cars$dist,
     xlab = "速度 (mph)", ylab = "刹车距离 (ft)",
     main = "刹车距离 vs 速度 —— 看起来是线性的吗？")
cor(cars$speed, cars$dist)   # 相关系数 ≈ 0.807

# summary讲解
summary(fit)




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
# 为什么是 n-2 而不是 n？
# 因为我们用数据估计了 2 个参数（beta0, beta1），
# 损失了 2 个自由度，剩下 n-2 个"独立"的残差信息用来估计 sigma^2。
# 这样 sigma2_hat 才是无偏估计：E(sigma2_hat) = sigma^2
sigma2_hat
sigma_hat <- sqrt(sigma2_hat)


# 验算手工计算的对不对
summary(fit)$sigma


# 标准误 用来衡量估计的参数beta的不确定性

# 计算beta_1 的标准误（Standard Error）
# SE(beta1_hat) = sigma_hat / sqrt(Sxx) ，Sxx = Sum((xi - x_bar)^2)
# beta1的标准误 和 beta1估计值没有关系，取决于数据噪声程度 和 x的分散程度
# 直觉：
#   - sigma_hat 越大（数据越嘈杂） → SE 越大（估计越不确定）
#   - Sxx 越大（x 越分散） → SE 越小（杠杆越长，斜率越稳）
# 启示：实验设计时，让 x 分散一些，能得到更精确的斜率估计。
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
# 关键区别：
#   置信区间 CI：估计"平均刹车距离 E(Y|x=25)"的不确定性 —— 只来自参数估计
#   预测区间 PI：预测"某一辆车的刹车距离 Y_new"的不确定性 
#                —— 参数估计不确定性 + 个体随机误差 sigma^2（就是那个多出来的 "1"）
# 所以 PI 永远比 CI 宽，而且即使 n → ∞，PI 宽度也不会收敛到 0。
SE_mean <- sigma_hat * sqrt(1/n + (x_new - x_bar)^2/Sxx)
SE_mean
# 构造均值的95%置信区间
CI_lower <- Y_pred - t_crit * SE_mean
CI_upper <- Y_pred + t_crit * SE_mean
c(CI_lower, CI_upper)
# 验证
predict(fit, new_data, interval = "confidence")



# 可视化
# 生成预测数据
speed_seq <- seq(4, 25, length.out = 100)
new_df <- data.frame(speed = speed_seq)

# 预测区间和置信区间
pred_PI <- predict(fit, new_df, interval = "prediction")
pred_CI <- predict(fit, new_df, interval = "confidence")

# 画图
plot(cars$speed, cars$dist,
     xlab = "速度 (mph)",
     ylab = "刹车距离 (ft)",
     main = "简单线性回归：置信区间 vs 预测区间")

# 回归线
abline(fit, col = "red", lwd = 2)

# 置信区间（窄）
lines(speed_seq, pred_CI[, 2], col = "blue", lty = 2, lwd = 2)
lines(speed_seq, pred_CI[, 3], col = "blue", lty = 2, lwd = 2)

# 预测区间（宽）
lines(speed_seq, pred_PI[, 2], col = "darkgreen", lty = 3, lwd = 2)
lines(speed_seq, pred_PI[, 3], col = "darkgreen", lty = 3, lwd = 2)

# 图例
legend("topleft",
       legend = c("回归线", "95%置信区间", "95%预测区间"),
       col = c("red", "blue", "darkgreen"),
       lty = c(1, 2, 3),
       lwd = 2)


# 模型诊断：检查假设是否成立
par(mfrow = c(2, 2))
plot(fit)
par(mfrow = c(1, 1))