# 场景：10个班级，每班20个学生，研究辅导对成绩的影响
library(lme4)
set.seed(123)

# 参数设定
n_class <- 10       # 10个班级
n_student <- 20     # 每班20人
beta_0 <- 60        # 总体截距（平均成绩）
beta_1 <- 5         # 辅导的效果（固定效应）
sigma_gamma <- 8    # 班级间标准差（随机效应）
sigma <- 10         # 个体误差标准差

# 生成数据
data <- data.frame()
for(i in 1:n_class){
  gamma_i <- rnorm(1, 0, sigma_gamma)  # 每个班的随机截距
  tutoring <- rbinom(n_student, 1, 0.5) # 是否接受辅导
  epsilon <- rnorm(n_student, 0, sigma) # 个体误差
  
  score <- beta_0 + beta_1 * tutoring + gamma_i + epsilon
  
  data <- rbind(data, data.frame(
    class = i,
    student = 1:n_student,
    tutoring = tutoring,
    score = score
  ))
}

data$class <- as.factor(data$class)



# 混合模型：固定效应 + 随机截距
model_mixed <- lmer(score ~ tutoring + (1|class), data = data)
summary(model_mixed)


# ICC = σ²_γ / (σ²_γ + σ²)
46.55 / (46.55 + 93.41)  # ≈ 0.33

# 普通线性回归（忽略班级）
model_lm <- lm(score ~ tutoring, data = data)
summary(model_lm)


# 把班级当固定效应
data$class <- as.factor(data$class)
model_fixed <- lm(score ~ tutoring + class, data = data)
summary(model_fixed)

library(ggplot2)

# 提取随机效应
ranef(model_mixed)$class

# 画图
ggplot(data, aes(x = factor(tutoring), y = score, color = class)) +
  geom_boxplot() +
  labs(x = "是否辅导", y = "成绩", title = "各班级成绩分布")