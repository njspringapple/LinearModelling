# 模拟更合理的数据 (n=50)
set.seed(123)
n <- 50
temp  <- runif(n, 50, 90)
hours <- runif(n, 1, 8)

# 用真实参数生成数据
eta_true <- -10 + 0.08 * temp + 0.6 * hours
pi_true  <- 1 / (1 + exp(-eta_true))
failure  <- rbinom(n, 1, pi_true)

data <- data.frame(temp, hours, failure)

# 拟合模型
model <- glm(failure ~ temp + hours, data = data, family = binomial)
summary(model)

# 参数
coef(model)

# 预测概率
predict(model, type = "response")

# Odds Ratio
exp(coef(model))





# ============================================
# 数据：10个人，1个自变量（年龄），Y=是否患病
# ============================================
x_raw <- c(25, 30, 35, 40, 45, 50, 55, 60, 65, 70)
Y     <- c( 0,  0,  0,  0,  1,  0,  1,  1,  1,  1)

# 设计矩阵（加截距列）
X <- cbind(1, x_raw)  # n x 2 矩阵
n <- nrow(X)
p <- ncol(X)

cat("===== 设计矩阵 X =====\n")
print(X)



# sigmoid 函数
sigmoid <- function(z) 1 / (1 + exp(-z))

# 初始值
beta <- c(0, 0)  # 从全0开始猜

for (iter in 1:20) {
  # 第1步：用当前 beta 算 pi
  pi_vec <- sigmoid(X %*% beta)
  
  # 第2步：构造 W 对角矩阵
  W <- diag(as.vector(pi_vec * (1 - pi_vec)))
  
  # 第3步：一阶导 Score
  Score <- t(X) %*% (Y - pi_vec)
  
  # 第4步：二阶导 Hessian（Fisher 信息矩阵）
  Fisher <- t(X) %*% W %*% X
  
  # 第5步：Newton-Raphson 更新
  beta_new <- beta + solve(Fisher) %*% Score
  
  # 检查收敛
  diff <- max(abs(beta_new - beta))
  cat(sprintf("迭代%2d: beta0=%8.4f, beta1=%8.5f, 变化=%e\n", 
              iter, beta_new[1], beta_new[2], diff))
  
  beta <- beta_new
  
  if (diff < 1e-8) {
    cat(sprintf(">>> 第 %d 次迭代收敛！\n", iter))
    break
  }
}

# 手工结果
cat("\n===== 手工计算结果 =====\n")
cat(sprintf("beta0 (截距) = %.6f\n", beta[1]))
cat(sprintf("beta1 (年龄) = %.6f\n", beta[2]))

# 协方差矩阵 = Fisher逆
pi_final <- sigmoid(X %*% beta)
W_final  <- diag(as.vector(pi_final * (1 - pi_final)))
Fisher_final <- t(X) %*% W_final %*% X
Cov_beta <- solve(Fisher_final)

cat("\n协方差矩阵 (Fisher逆):\n")
print(Cov_beta)

# 标准误
SE <- sqrt(diag(Cov_beta))
cat(sprintf("\nSE(beta0) = %.6f\n", SE[1]))
cat(sprintf("SE(beta1) = %.6f\n", SE[2]))

# 95% Wald 置信区间
z_crit <- qnorm(0.975)  # 1.96

CI_beta0 <- beta[1] + c(-1, 1) * z_crit * SE[1]
CI_beta1 <- beta[2] + c(-1, 1) * z_crit * SE[2]

cat(sprintf("\nbeta0 的 95%% CI: (%.4f, %.4f)\n", CI_beta0[1], CI_beta0[2]))
cat(sprintf("beta1 的 95%% CI: (%.5f, %.5f)\n", CI_beta1[1], CI_beta1[2]))

# 优势比 OR 及其置信区间
OR <- exp(beta[2])
CI_OR <- exp(CI_beta1)
cat(sprintf("\nOR(年龄) = %.4f\n", OR))
cat(sprintf("OR 的 95%% CI: (%.4f, %.4f)\n", CI_OR[1], CI_OR[2]))

# ============================================
# 方法二：R 内置 glm 函数
# ============================================

cat("\n\n===== glm 函数结果 =====\n")
fit <- glm(Y ~ x_raw, family = binomial(link = "logit"))
summary(fit)

# 置信区间
cat("\nbeta 的 95% CI (glm):\n")
print(confint.default(fit))  # Wald 置信区间

cat("\nOR 及 95% CI:\n")
cat(sprintf("OR = %.4f\n", exp(coef(fit)[2])))
cat(sprintf("OR 的 95%% CI: (%.4f, %.4f)\n", 
            exp(confint.default(fit)[2,1]), 
            exp(confint.default(fit)[2,2])))

# ============================================
# 对比
# ============================================
cat("\n\n===== 对比两种方法 =====\n")
cat(sprintf("          手工计算     glm函数\n"))
cat(sprintf("beta0:   %9.6f    %9.6f\n", beta[1], coef(fit)[1]))
cat(sprintf("beta1:   %9.6f    %9.6f\n", beta[2], coef(fit)[2]))
cat(sprintf("SE(b1):  %9.6f    %9.6f\n", SE[2], summary(fit)$coef[2,2]))








# 大模型
fit_full <- glm(Y ~ x_raw, family = binomial)

# 小模型（只有截距）
fit_null <- glm(Y ~ 1, family = binomial)

# 似然比检验
LQ <- -2 * (logLik(fit_null) - logLik(fit_full))
p_value <- 1 - pchisq(LQ, df = 1)

















# ============================================
# 1. 生成模拟数据
# ============================================
set.seed(42)
n <- 500

# 两个特征, 例如年龄，体重
x1 <- rnorm(n)
x2 <- rnorm(n)

sigmoid <- function(z) 1 / (1 + exp(-z))

# 真实关系：logit(P) = -1 + 2*x1 + 0.5*x2
#prob <- plogis(-1 + 2*x1 + 0.5*x2)  # plogis = logistic函数
prob <- sigmoid(-1 + 2*x1 + 0.5*x2)  # 
y <- rbinom(n, 1, prob)

data <- data.frame(y, x1, x2)

# 查看正负样本比例
table(data$y)

# ============================================
# 2. 拆分训练集和测试集
# ============================================
train_idx <- sample(1:n, 350)
train <- data[train_idx, ]
test  <- data[-train_idx, ]

# ============================================
# 3. 拟合逻辑回归
# ============================================
model <- glm(y ~ x1 + x2, data = train, family = binomial)
summary(model)

# ============================================
# 4. 在测试集上预测（预测分数）
# ============================================
pred_score <- predict(model, newdata = test, type = "response")

# 看几个预测分数
head(data.frame(真实值 = test$y, 预测分数 = round(pred_score, 3)))

# ============================================
# 5. 画 ROC 曲线 + 计算 AUC
# ============================================
# 安装包（如果没有的话）
# install.packages("pROC")
library(pROC)

roc_obj <- roc(test$y, pred_score)

# 画ROC曲线
plot(roc_obj, 
     col = "blue", 
     lwd = 2,
     main = "ROC 曲线",
     xlab = "FPR (1 - 特异度)",
     ylab = "TPR (敏感度)")
abline(a = 0, b = 1, lty = 2, col = "gray")  # 对角线（随机猜）

# 输出AUC
auc_value <- auc(roc_obj)
cat("AUC =", round(auc_value, 4), "\n")

# 在图上标注AUC
legend("bottomright", 
       legend = paste("AUC =", round(auc_value, 4)),
       col = "blue", lwd = 2)


