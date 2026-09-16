# Step1 简单线性回归

# 数据：年龄和降压效果
age <- c(30, 40, 50, 60, 70)
effect <- c(10, 15, 18, 22, 25)

# 拟合模型
model1 <- lm(effect ~ age)
summary(model1)

# 设计矩阵长什么样
model.matrix(model1)

# 画图
plot(age, effect, pch = 19, main = "简单线性回归")
abline(model1, col = "red")

# Step2 多元线性回归


# 数据：年龄、体重和降压效果
age <- c(30, 40, 50, 60, 70, 35, 45, 55)
weight <- c(60, 70, 80, 75, 85, 65, 72, 78)
effect <- c(10, 15, 20, 18, 25, 12, 16, 19)

# 拟合模型
model2 <- lm(effect ~ age + weight)
summary(model2)

# 设计矩阵
model.matrix(model2)

# Step3 单因素ANOVA（分类变量）

# 数据：三种药，每种6人
drug <- factor(rep(c("A", "B", "C"), each = 6))
effect <- c(10, 12, 11, 13, 10, 12,   # A组
            20, 22, 21, 23, 20, 22,   # B组
            15, 17, 16, 18, 15, 17)   # C组

# 方法1：用aov
model3a <- aov(effect ~ drug)
summary(model3a)

# 方法2：用lm（完全等价）
model3b <- lm(effect ~ drug)
summary(model3b)

# 设计矩阵（参考编码，A是参考组）
model.matrix(model3b)

# 看每组均值
tapply(effect, drug, mean)


# Step4 双因素ANOVA（无交互）

# 数据：药物（A/B）× 性别（男/女）
drug <- factor(rep(c("A", "B"), each = 6))
gender <- factor(rep(c("男", "女"), times = 6))
effect <- c(10, 12, 11, 13, 14, 15,   # 药A
            20, 22, 21, 23, 24, 25)   # 药B

# 无交互模型
model4 <- lm(effect ~ drug + gender)
summary(model4)

# 设计矩阵
model.matrix(model4)

# Step5 双因素ANOVA（有交互）

# 数据：药物效果因性别而异
drug <- factor(rep(c("A", "B"), each = 6))
gender <- factor(rep(c("男", "女"), times = 6))
effect <- c(10, 20, 11, 21, 12, 22,   # 药A：男女差10
            15, 35, 16, 36, 17, 37)   # 药B：男女差20（交互！）

# 有交互模型
model5 <- lm(effect ~ drug * gender)  # drug * gender = drug + gender + drug:gender
summary(model5)

# 设计矩阵
model.matrix(model5)

# 交互图
interaction.plot(drug, gender, effect, col = c("blue", "red"), lwd = 2)

# Step6 协方差分析（离散+连续，无交互）
# 数据：三组，每组有年龄和降压效果
group <- factor(rep(c("A", "B", "C"), each = 5))
age <- c(30, 40, 50, 60, 70,    # A组
         35, 45, 55, 65, 75,    # B组
         32, 42, 52, 62, 72)    # C组
effect <- c(12, 16, 20, 24, 28,  # A组
            17, 21, 25, 29, 33,  # B组
            22, 26, 30, 34, 38)  # C组

# 无交互：平行线
model6 <- lm(effect ~ group + age)
summary(model6)

# 设计矩阵
model.matrix(model6)

# 画图
library(ggplot2)
df <- data.frame(group, age, effect)
ggplot(df, aes(x = age, y = effect, color = group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("无交互：平行线")


# Step7 协方差分析（离散+连续，有交互）
# 数据：斜率不同
group <- factor(rep(c("A", "B", "C"), each = 5))
age <- c(30, 40, 50, 60, 70,
         35, 45, 55, 65, 75,
         32, 42, 52, 62, 72)
effect <- c(10, 14, 18, 22, 26,  # A组：斜率0.4
            15, 21, 27, 33, 39,  # B组：斜率0.6
            20, 30, 40, 50, 60)  # C组：斜率1.0

# 有交互：不同斜率
model7 <- lm(effect ~ group * age)
summary(model7)

# 设计矩阵
model.matrix(model7)

# 画图
df <- data.frame(group, age, effect)
ggplot(df, aes(x = age, y = effect, color = group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("有交互：斜率不同")

# Step8 平方和类型
library(car)

# 双因素数据
set.seed(456)
A <- factor(rep(c("a1", "a2"), each = 10))
B <- factor(rep(c("b1", "b2"), times = 10))
y <- c(10, 12, 11, 13, 14, 20, 22, 21, 23, 24,
       15, 17, 16, 18, 19, 30, 32, 31, 33, 34) + rnorm(20, 0, 1)

# 不平衡设计（删几个观测）
df <- data.frame(A, B, y)[-c(1, 5, 15), ]

model <- lm(y ~ A + B, data = df)

# Type I
cat("=== Type I: 先A后B ===\n")
anova(model)

# 换顺序
model_rev <- lm(y ~ B + A, data = df)
cat("\n=== Type I: 先B后A ===\n")
anova(model_rev)  # 不一样！

# Type III
cat("\n=== Type III ===\n")
Anova(model, type = 3)  # 顺序无关

# Step9 前后对比（临床试验）
# 数据
group <- factor(c(rep("真药", 10), rep("安慰剂", 10)))
x1 <- c(150, 145, 160, 155, 140, 148, 152, 158, 143, 147,  # 真药组治疗前
        148, 152, 145, 150, 155, 142, 157, 149, 153, 146)  # 安慰剂组治疗前
x2 <- c(130, 128, 140, 135, 125, 132, 138, 142, 128, 130,  # 真药组治疗后
        145, 150, 143, 148, 152, 140, 154, 147, 150, 144)  # 安慰剂组治疗后

# 方法1：差值t检验
d <- x2 - x1
t.test(d ~ group)

# 方法2：ANCOVA
model9 <- lm(x2 ~ x1 + group)
summary(model9)

# 看group的系数：真药比安慰剂多降多少

# Step10 查看设计矩阵的技巧

# 先创建factor变量
group <- factor(c("A", "B", "C", "A", "B", "C"))

# 查看默认编码（treatment coding）
cat("=== 默认编码（参考编码）===\n")
contrasts(group)

# 改成sum coding（效应编码）
contrasts(group) <- contr.sum(3)
cat("\n=== Sum coding（效应编码）===\n")
contrasts(group)

# 改成helmert coding
contrasts(group) <- contr.helmert(3)
cat("\n=== Helmert coding ===\n")
contrasts(group)








