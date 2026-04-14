data(cars)
# cars 是list，如果 speed, dist 两列又是 list，所以需要用 [["xxx"]] 来取出列的元素值（向量）
plot(cars[["speed"]],cars[["dist"]],xlab="speed",ylab="stopping distance",main="Average trend plus random scatter")
# fit 类型是list类型，可以用 typeof(fit) 或者 is.list(fit) 来查看
# R 语言默认回归算法 是 OLS（最小二乘法），dist ~ speed 是 公式语法，左侧是要预测的Y，波浪线是”被建模“的意思，右侧是用于预测的 X
#  dist~speed 意思是通过已有数据，要预测不同速度下的刹车举例distance
fit <- lm(dist~speed, data=cars)
# lm 类型，直接打印用的lm类的打印方式组织输出
print(fit)
# 强制用list类型打印
print.default(fit)
abline(fit,col="red",lwd=4)
# 提取器函数，用于从回归模型中提取系数：Intercept ， speed （截距和斜率）
coef(fit)



data(mtcars)
# mpg 油耗， wt 车重， 根据车重预测油耗
fit <- lm(mpg~wt,data=mtcars)
summary(fit)
# 计算回归系数的置信区间 (Confidence Interval)
confint(fit)
# 用训练好的模型，对新数据做预测，fit 就是训练好的模型
# interval= "prediction" 不仅要预测值，还要给个预测区间
# 预测车重是2500磅 3500磅的车重量
# 输出 lwr，upr 是下限和上限，fit 是均值
predict(fit,newdata = data.frame(wt=c(2.5,3.5)), interval= "prediction")
# 绘制回归线和置信带
plot(mtcars$wt, mtcars$mpg,
     xlab="Weight (1000 lbs)",
     ylab="Miles per Gallon",
     main="mpg ~ wt 回归分析")

abline(fit, col="red", lwd=2)

# 添加置信区间带
new_wt <- data.frame(wt=seq(1.5, 5.5, length=100))
# pred 返回的是矩阵，不是数据框
pred <- predict(fit, new_wt, interval="confidence")
lines(new_wt$wt, pred[,2], col="blue", lty=2)
lines(new_wt$wt, pred[,3], col="blue", lty=2)

legend("topright", 
       legend=c("回归线", "95% 置信带"),
       col=c("red", "blue"), lty=c(1, 2))




data(iris)

# RSE ： Residual standard error 残差，预测值与实际值的平均误差 越小越好
# R² (Multiple R-squared) 决定系数 模型解释了 Y 的多少百分比的变化（解释力多强） 0-1，越大越好

# 连续响应 + 连续解释变量
# 根据花的花瓣长度预测花萼长度
fit_lm <- lm(Sepal.Length~Petal.Length, data=iris)
# 预测值 = 截距 + 斜率 × 花瓣长度
# (Intercept)   4.30660    0.07839   54.94   <2e-16 ***
# 截距行，也就是 花瓣 长度为0时候，花萼长度预测为4.31 cm
# Petal.Length  0.40892    0.01889   21.65   <2e-16 ***
# 连续系统的系数是斜率
# 花瓣长度每增加1cm，花萼长度平均增加0.41cm
# Residual standard error: 0.4071 on 148 degrees of freedom
# Multiple R-squared:   0.76,	Adjusted R-squared:  0.7583 
summary(fit_lm)

# 连续响应 + 分类解释变量
# 根据花的物种预测花萼长度
fit_group <- lm(Sepal.Length~Species, data=iris)
# R 自动将分类变量转换为虚拟变量 (dummy variables)，并选择第一个类别作为基准组 (baseline)
# (Intercept)         5.0060     0.0728  68.762  < 2e-16 ***
# setosa 物种的平均花萼长度
# Speciesversicolor   0.9300     0.1030   9.033 8.77e-16 ***
# versicolor 比 setosa 长 0.93
# p < 0.001 看 Pr(>|t|) 这一列，说明差异 极显著， p 越小越显著，*号越多越显著
# Multiple R-squared:  0.6187
# R-squared 物种能解释 61.9% 的花萼长度变化
# Residual standard error: 0.5148
# RSE = 0.5148 预测误差约 ±0.51 cm
# 分类变量的系数不是传统斜率，而是"相对于基准组的差值"，所以只有 k-1 行（k=类别数）。
summary(fit_group)

# 二元响应：是否setosa
# 根据花瓣长度，预测这朵花是不是 setosa
iris[["is_setosa"]] <- ifelse(iris[["Species"]] == "setosa" , 1, 0)
# family = binomial() 它告诉 R：Y 是二项分布，使用逻辑回归
fit_logit <- glm(is_setosa~Petal.Length, data=iris, family = binomial())
# (Intercept)     91.67   47334.35   0.002    0.998
# Petal.Length   -37.22   18357.58  -0.002    0.998
# p值接近1，不可信
# Residual deviance: 7.3324e-09  on 148  degrees of freedom
# 残差接近 0，过拟合
summary(fit_logit)

plot(iris$Petal.Length, iris$is_setosa,
     xlab="Petal.Length", ylab="is_setosa")


# sigmod 逻辑回归函数，用来把线性预测值转换成概率**，也就是Y是二元变量
z <- seq(-10, 10, length=100)
p <- 1 / (1 + exp(-z))

plot(z, p, type="l", lwd=2,
     xlab="线性预测值 z", 
     ylab="概率 P(Y=1)",
     main="逻辑函数：把任意 z 映射到 [0,1]")
abline(h=c(0,1), col="gray", lty=2)
abline(v=0, col="red", lty=2)
points(0, 0.5, col="red", pch=16)

















