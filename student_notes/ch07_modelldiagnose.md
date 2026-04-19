# 第 7 章学生笔记：模型诊断 / Kapitel 7: Modelldiagnose

材料 / Grundlage:

- 讲义 / Folien: `chapter_pdfs/LinearModellingPPT/07_Modelldiagnose.pdf`
- 主教材 / Hauptlehrbuch: `Regression.pdf`, Chapter 3.4.4, Chapter 4.1.3, Chapter 4.1.4, Chapter 4.2
- 辅助教材 / Hilfslehrbuch: Faraway, `Linear Models with R`, Chapter 4, Chapter 5, Chapter 6
- 扩展教材 / Erweiterung: Faraway, `Extending The Linear Model With R`, Chapter 8.4 as GLM diagnostics preview

## 覆盖审核 / Abdeckungspruefung

| 来源 / Quelle | 内容点 / Inhaltspunkt | 笔记位置 / Wo |
|---|---|---|
| 讲义 7 | 回归问题与诊断总览 | 7.1 |
| 讲义 7 | 残差类型：普通、标准化、学生化、交叉验证、递归残差 | 7.2 |
| 讲义 7 | 非正态误差 | 7.3 |
| 讲义 7 | 异方差 | 7.4 |
| 讲义 7 | 误差相关、Durbin-Watson | 7.5 |
| 讲义 7 | 离群点、高杠杆点、Cook's distance | 7.6 |
| 讲义 7 | 模型设定错误、partial leverage plot | 7.7 |
| 讲义 7 | 共线性、VIF、条件数 | 7.8 |
| 讲义 7 | 测量误差 | 7.9 |
| 主教材 | 线性模型诊断、异方差、自相关、正则化背景 | 全章 |
| Faraway | R 诊断图、influence measures、预测变量/误差问题 | R 练习 |

## 7.1 为什么要诊断 / Warum Modelldiagnose?

**遇到的问题 / Problem.**  
线性模型输出表看起来很正式，但它依赖线性、同方差、独立、正态性和正确模型设定等假设。如果假设严重不成立，系数解释、标准误、$p$ 值和预测区间都会失真。

**可能的解决路径 / Moegliche Wege.**

- 只看 `summary()`：快但危险。
- 系统检查残差、杠杆、影响点和共线性：能定位问题。
- 发现问题后选择变换、WLS/GLS、稳健方法、非线性项或 mixed model。

**本课采用的方法及好处 / Methode und Vorteil.**  
诊断不是为了“证明模型完美”，而是为了回答：

- 哪条假设最可疑？
- 这个问题会影响解释、推断还是预测？
- 有哪些合理修正路径？

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)
par(mfrow = c(2, 2))
plot(fit)
par(mfrow = c(1, 1))
```

## 7.2 残差类型 / Verschiedene Typen von Residuen

**遇到的问题 / Problem.**  
普通残差 $\widehat{\varepsilon}_i$ 的方差不完全相同，高杠杆点的残差天然更小。因此直接比较普通残差可能误导。

**本课采用的方法及好处 / Methode und Vorteil.**  
普通残差：

$$
\widehat{\varepsilon}_i=Y_i-\widehat{Y}_i
$$

Hat value：

$$
h_{ii}=P_{ii}
$$

标准化残差：

$$
r_i=\frac{\widehat{\varepsilon}_i}{\widehat{\sigma}\sqrt{1-h_{ii}}}
$$

学生化残差进一步使用删除第 $i$ 个观测后的方差估计。交叉验证残差可写为：

$$
\widehat{\varepsilon}_{(i)}=\frac{\widehat{\varepsilon}_i}{1-h_{ii}}
$$

好处是：不同残差类型服务不同诊断目的。

**具体解决 / Konkrete Loesung.**

- 普通残差：看拟合偏差。
- 标准化/学生化残差：找异常响应。
- 交叉验证残差：看单点预测误差。
- 递归残差：常用于时间序列或结构变化检查。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)

head(residuals(fit))
head(rstandard(fit))
head(rstudent(fit))
head(residuals(fit) / (1 - hatvalues(fit)))
```

## 7.3 非正态误差 / Nicht-normalverteilte Stoerterme

**遇到的问题 / Problem.**  
正态性主要影响小样本精确推断。如果误差有重尾或明显偏态，$t$ 检验和预测区间可能不可靠。

**可能的解决路径 / Moegliche Wege.**

- QQ 图检查。
- 变换响应变量。
- 稳健回归。
- 大样本下依赖渐近近似，但仍要谨慎。

**本课采用的方法及好处 / Methode und Vorteil.**  
QQ 图将残差分位数与正态分位数比较。若点大致在直线上，正态性较合理；若尾部严重偏离，可能有重尾或异常点。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)
qqnorm(rstandard(fit))
qqline(rstandard(fit), col = "red")
shapiro.test(residuals(fit))
```

## 7.4 异方差 / Heterogene Varianzen

**遇到的问题 / Problem.**  
如果误差方差随 fitted values 或某个协变量改变，OLS 系数可能仍无偏，但标准误、检验和区间会不可靠。

**本课采用的方法及好处 / Methode und Vorteil.**  
同方差假设：

$$
V(\varepsilon_i)=\sigma^2
$$

异方差：

$$
V(\varepsilon_i)=\sigma_i^2
$$

残差图中常见信号是扇形扩散或收缩。修正路径包括变换、加权最小二乘和稳健标准误。

**R 练习 / R-Uebung.**

```r
data(cars)
fit <- lm(dist ~ speed, data = cars)
plot(fitted(fit), residuals(fit),
     xlab = "fitted values",
     ylab = "residuals")
abline(h = 0, col = "gray")

plot(fitted(fit), sqrt(abs(rstandard(fit))),
     xlab = "fitted values",
     ylab = "sqrt absolute standardized residuals")
```

## 7.5 相关误差与 Durbin-Watson / Autokorrelation

**遇到的问题 / Problem.**  
时间序列、空间数据或重复测量中，误差可能相关。若忽略相关性，标准误通常会错误。

**本课采用的方法及好处 / Methode und Vorteil.**  
独立假设：

$$
\operatorname{Cov}(\varepsilon_i,\varepsilon_j)=0,\quad i\neq j
$$

一阶自相关常写为：

$$
\varepsilon_t=\rho\varepsilon_{t-1}+u_t
$$

Durbin-Watson 统计量基于相邻残差差：

$$
DW=
\frac{\sum_{t=2}^n(\widehat{\varepsilon}_t-\widehat{\varepsilon}_{t-1})^2}
{\sum_{t=1}^n\widehat{\varepsilon}_t^2}
$$

值接近 $2$ 表示一阶自相关较弱；明显低于 $2$ 常提示正自相关。

**R 练习 / R-Uebung.**

```r
data(Nile)
t <- seq_along(Nile)
fit <- lm(as.numeric(Nile) ~ t)

plot(t, residuals(fit), type = "b")
acf(residuals(fit))

res <- residuals(fit)
DW <- sum(diff(res)^2) / sum(res^2)
DW
```

## 7.6 离群点、高杠杆点和 Cook's distance

**遇到的问题 / Problem.**  
有些点的 $Y$ 很异常，有些点的 $x$ 位置很极端，还有些点会强烈改变拟合结果。它们不是同一类问题。

**本课采用的方法及好处 / Methode und Vorteil.**

- 离群点 / Ausreisser: 响应方向异常，残差大。
- 高杠杆点 / hoher Leverage: $x$ 空间位置极端，$h_{ii}$ 大。
- 强影响点 / einflussreicher Punkt: 删除后模型显著改变。

Cook's distance 结合残差和杠杆衡量影响：

$$
D_i\ \text{large} \Rightarrow \text{observation may be influential}
$$

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)

head(hatvalues(fit))
head(cooks.distance(fit))

plot(hatvalues(fit), rstandard(fit),
     xlab = "hat values",
     ylab = "standardized residuals")
abline(h = c(-2, 2), col = "gray", lty = 2)
```

## 7.7 模型设定错误与 partial leverage plot

**遇到的问题 / Problem.**  
如果真实关系非线性或遗漏了重要交互，残差图会显示系统模式。系数显著不能弥补模型设定错误。

**本课采用的方法及好处 / Methode und Vorteil.**  
partial leverage plot 通过调整其他协变量，显示某个变量的条件关系。若调整后仍有曲线形状，说明可能需要变换、多项式或 spline。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)

res_y <- residuals(lm(mpg ~ hp + cyl, data = mtcars))
res_wt <- residuals(lm(wt ~ hp + cyl, data = mtcars))

plot(res_wt, res_y,
     xlab = "wt adjusted for hp and cyl",
     ylab = "mpg adjusted for hp and cyl")
abline(lm(res_y ~ res_wt), col = "red", lwd = 2)
```

## 7.8 共线性 / Kollinearitaet

**遇到的问题 / Problem.**  
解释变量之间高度相关时，单个系数的标准误会变大，系数可能不稳定。模型整体预测可能还不错，但单个参数解释变难。

**本课采用的方法及好处 / Methode und Vorteil.**  
VIF 思想：

$$
VIF_j=\frac{1}{1-R_j^2}
$$

其中 $R_j^2$ 来自把 $x_j$ 对其他解释变量回归。$VIF_j$ 越大，说明 $x_j$ 被其他变量解释得越好，独立信息越少。

条件数 / condition number 也用于检查设计矩阵接近奇异的程度。

**R 练习 / R-Uebung.**

```r
data(mtcars)
X <- model.matrix(lm(mpg ~ wt + hp + cyl, data = mtcars))[, -1]

vif_manual <- sapply(seq_len(ncol(X)), function(j) {
  xj <- X[, j]
  others <- X[, -j, drop = FALSE]
  r2 <- summary(lm(xj ~ others))[["r.squared"]]
  1 / (1 - r2)
})
names(vif_manual) <- colnames(X)
vif_manual

kappa(model.matrix(lm(mpg ~ wt + hp + cyl, data = mtcars)))
```

## 7.9 测量误差 / Messfehler

**遇到的问题 / Problem.**  
如果解释变量 $x$ 本身有测量误差，普通线性模型把 $x$ 当作精确观测的假设被破坏。斜率可能被系统性低估或产生其他偏差。

**可能的解决路径 / Moegliche Wege.**

- 改进测量设计。
- 使用重复测量估计测量误差。
- 使用 errors-in-variables model。
- 做敏感性分析。

**本课采用的方法及好处 / Methode und Vorteil.**  
讲义在本章主要提醒风险：诊断不只看 $Y$ 方向的残差，也要思考解释变量是否可靠。

## 7.10 诊断后的处理路线 / Was tun nach der Diagnose?

| 诊断发现 | 可能处理 |
|---|---|
| 非线性 | 变换、多项式、spline |
| 异方差 | 变换、WLS、稳健标准误 |
| 自相关 | GLS、时间序列误差结构、mixed model |
| 离群点 | 核查数据、稳健回归、敏感性分析 |
| 高杠杆强影响 | 报告影响、比较删除前后结果 |
| 共线性 | 合并变量、删减变量、主成分、ridge/LASSO |
| 测量误差 | 改进测量、误差模型 |

## 7.11 选择题 / Multiple-Choice-Fragen

1. 高杠杆点主要指什么？  
   A. 在 $x$ 空间位置极端的点  
   B. $Y$ 一定最大  
   C. 残差一定为零  
   D. 样本量太小  
   答案 / Antwort: A

2. 异方差主要影响什么？  
   A. 标准误和检验可靠性  
   B. 变量名  
   C. 样本编号  
   D. 截距是否存在  
   答案 / Antwort: A

3. VIF 大说明什么？  
   A. 某解释变量可被其他解释变量很好解释  
   B. 误差一定正态  
   C. 没有共线性  
   D. $R^2$ 必定为零  
   答案 / Antwort: A

4. Cook's distance 用来衡量什么？  
   A. 单个观测对拟合结果的影响  
   B. 样本均值  
   C. 正态分布分位数  
   D. 变量单位  
   答案 / Antwort: A

5. Durbin-Watson 主要用于检查什么？  
   A. 一阶自相关  
   B. 分类变量编码  
   C. dummy trap  
   D. 多重检验  
   答案 / Antwort: A

## 7.12 本章总结 / Zusammenfassung

中文总结：

- 模型诊断检查线性、同方差、独立、正态性、异常点、共线性和测量误差等问题。
- 残差有不同类型，标准化/学生化残差更适合比较异常程度。
- 异方差主要破坏标准误和检验；自相关常见于时间或重复测量数据。
- 高杠杆点不一定影响大，Cook's distance 同时考虑杠杆和残差。
- 共线性让单个系数不稳定，VIF 和条件数是常用诊断。
- 诊断后要根据问题选择变换、WLS/GLS、稳健方法、mixed model 或正则化等。

Deutsche Zusammenfassung:

- Modelldiagnose prueft zentrale Annahmen und auffaellige Beobachtungen.
- Standardisierte und studentisierte Residuen helfen beim Erkennen von Ausreissern.
- Heteroskedastizitaet und Autokorrelation betreffen vor allem Standardfehler und Inferenz.
- Leverage und Cook's Distanz messen unterschiedliche Aspekte einflussreicher Punkte.
- Kollinearitaet erschwert die Interpretation einzelner Koeffizienten.

## 7.13 专业德语单词汇总 / Fachwortschatz

| Deutsch | 中文 | 说明 |
|---|---|---|
| Modelldiagnose | 模型诊断 | 检查假设和影响点 |
| Residuum | 残差 | $\widehat{\varepsilon}_i$ |
| standardisiertes Residuum | 标准化残差 | 调整方差 |
| studentisiertes Residuum | 学生化残差 | 删除点方差 |
| Kreuzvalidierungs-Residuum | 交叉验证残差 | leave-one-out |
| Heteroskedastizitaet | 异方差 | 方差不恒定 |
| Autokorrelation | 自相关 | 误差相关 |
| Durbin-Watson-Test | Durbin-Watson 检验 | 一阶自相关 |
| Ausreisser | 离群点 | 响应异常 |
| Leverage | 杠杆值 | $h_{ii}$ |
| Cook's Distanz | Cook 距离 | 影响度 |
| Kollinearitaet | 共线性 | 解释变量相关 |
| Varianzinflationsfaktor | 方差膨胀因子 | VIF |
| Messfehler | 测量误差 | 变量测量不准 |
