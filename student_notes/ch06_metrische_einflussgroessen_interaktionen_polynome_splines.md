# 第 6 章学生笔记：连续解释变量、交互、多项式与样条 / Kapitel 6: Metrische Einflussgroessen

材料 / Grundlage:

- 讲义 / Folien: `chapter_pdfs/LinearModellingPPT/06_Metrische_Einflugroen_Interaktionen_Polynomiale_Regression_Trigonometrische_Polynome_Regression.pdf`
- 主教材 / Hauptlehrbuch: `Regression.pdf`, Chapter 3.1.3, Chapter 8.1, Chapter 9.3
- 辅助教材 / Hilfslehrbuch: Faraway, `Linear Models with R`, Chapter 7, Chapter 13
- 扩展教材 / Erweiterung: Faraway, `Extending The Linear Model With R`, Chapter 14, Chapter 15

## 覆盖审核 / Abdeckungspruefung

| 来源 / Quelle | 内容点 / Inhaltspunkt | 笔记位置 / Wo |
|---|---|---|
| 讲义 6 | 连续变量交互 | 6.1 |
| 讲义 6 | 连续解释变量的处理策略 | 6.2 |
| 讲义 6 | 变量变换、多项式、trigonometric polynomials | 6.3, 6.4 |
| 讲义 6 | basis functions / 基函数一般框架 | 6.5 |
| 讲义 6 | trend model | 6.5 |
| 讲义 6 | 多重检验、同时置信区间 | 6.6 |
| 主教材 | 非线性协变量效应、splines、interactions | 全章 |
| Faraway | transformation、ANCOVA 和 R 中 `poly()`、`I()` | R 练习 |

## 6.1 连续变量交互 / Interaktionen bei metrischen Variablen

**遇到的问题 / Problem.**  
两个连续解释变量可能不是简单相加。例如汽车重量对油耗的影响可能随马力不同而变化。

**可能的解决路径 / Moegliche Wege.**

- 只用加性模型：假设每个变量作用不依赖其他变量。
- 加入乘积项：允许一个变量的斜率随另一个变量改变。

**本课采用的方法及好处 / Methode und Vorteil.**  
交互模型：

$$
E(Y\mid x,z)=\beta_0+\beta_1x+\beta_2z+\beta_3xz
$$

$x$ 的边际斜率为：

$$
\frac{\partial E(Y\mid x,z)}{\partial x}=\beta_1+\beta_3z
$$

好处是：模型仍对参数线性，但可以表达“作用随条件改变”。

**具体解决 / Konkrete Loesung.**  
解释交互时，不要孤立解释 $\beta_1$ 或 $\beta_2$。要说明在给定 $z$ 水平下，$x$ 的作用是多少。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit_add <- lm(mpg ~ wt + hp, data = mtcars)
fit_int <- lm(mpg ~ wt * hp, data = mtcars)

anova(fit_add, fit_int)
coef(fit_int)
```

## 6.2 连续解释变量的处理策略 / Behandlung metrischer Einflussgroessen

**遇到的问题 / Problem.**  
连续变量的关系未必是直线。若直接放入 $x$，可能错设函数形式。

**可能的解决路径 / Moegliche Wege.**

- 线性项 $x$：最简单、最可解释。
- 变量变换 $g(x)$：如 $\log(x)$、$1/x$。
- 多项式 $x,x^2,x^3$：表达平滑曲率。
- 分段函数和样条：局部灵活。
- 非参数或 additive model：更灵活，后续章节展开。

**本课采用的方法及好处 / Methode und Vorteil.**  
仍保持线性模型框架：

$$
E(Y\mid x)=\beta_0+\beta_1b_1(x)+\dots+\beta_kb_k(x)
$$

其中 $b_j(x)$ 是基函数。好处是：模型对参数 $\beta$ 仍是线性的，所以可继续用最小二乘。

**具体解决 / Konkrete Loesung.**  
先画散点图和残差图。若弯曲明显，再考虑变换、多项式或 spline，而不是只机械看 $p$ 值。

**R 练习 / R-Uebung.**

```r
data(women)
fit_linear <- lm(weight ~ height, data = women)
fit_quad <- lm(weight ~ height + I(height^2), data = women)

anova(fit_linear, fit_quad)
plot(women[["height"]], women[["weight"]])
abline(fit_linear, col = "red", lwd = 2)
h <- seq(min(women[["height"]]), max(women[["height"]]), length.out = 100)
lines(h, predict(fit_quad, newdata = data.frame(height = h)),
      col = "blue", lwd = 2)
```

## 6.3 变量变换 / Transformationen

**遇到的问题 / Problem.**  
关系可能在原始尺度上非线性，但在变换尺度上近似线性。

**本课采用的方法及好处 / Methode und Vorteil.**  
常见形式：

$$
E(Y\mid z)=\beta_0+\beta_1\log(z)
$$

或：

$$
E(Y\mid z)=\beta_0+\beta_1\frac{1}{z}
$$

若变换响应变量：

$$
\log(Y)=\beta_0+\beta_1x+\varepsilon
$$

解释会改变，因为：

$$
E[\log(Y)\mid x]\neq \log(E[Y\mid x])
$$

**具体解决 / Konkrete Loesung.**  
变换解释变量时，仍解释 $Y$ 的均值随变换后变量变化；变换响应变量时，要特别谨慎反变换预测和参数解释。

**R 练习 / R-Uebung.**

```r
data(cars)
fit1 <- lm(dist ~ speed, data = cars)
fit2 <- lm(log(dist) ~ speed, data = cars)

par(mfrow = c(1, 2))
plot(fitted(fit1), residuals(fit1))
abline(h = 0, col = "gray")
plot(fitted(fit2), residuals(fit2))
abline(h = 0, col = "gray")
par(mfrow = c(1, 1))
```

## 6.4 多项式与三角多项式 / Polynome und trigonometrische Polynome

**遇到的问题 / Problem.**  
有些关系有曲率或周期性。单一直线无法表达。

**本课采用的方法及好处 / Methode und Vorteil.**  
多项式模型：

$$
E(Y\mid x)=\beta_0+\beta_1x+\beta_2x^2+\dots+\beta_kx^k
$$

它关于 $x$ 非线性，但关于参数 $\beta$ 线性。周期性关系可用三角基函数：

$$
E(Y\mid t)=\beta_0+\sum_{j=1}^k
\left[a_j\sin(jt)+b_j\cos(jt)\right]
$$

好处是：基函数扩展后仍能用普通线性模型估计。

**具体解决 / Konkrete Loesung.**  
多项式阶数越高越灵活，但越容易过拟合和边界震荡。实际中要结合图形、交叉验证或后续模型选择。

**R 练习 / R-Uebung.**

```r
data(women)
fit_raw <- lm(weight ~ height + I(height^2), data = women)
fit_orth <- lm(weight ~ poly(height, 2), data = women)

summary(fit_raw)
summary(fit_orth)
```

练习 / Uebung: `poly(height, 2)` 和 `I(height^2)` 的拟合值是否通常相同？  
答案 / Antwort: 若模型空间相同，拟合值相同；参数解释不同，因为 `poly()` 默认使用正交多项式。

## 6.5 基函数与 trend model / Basisfunktionen

**遇到的问题 / Problem.**  
多项式、三角多项式、样条看似不同，其实都可以统一成基函数展开。

**本课采用的方法及好处 / Methode und Vorteil.**  
一般形式：

$$
E(Y\mid x)=\beta_0+\sum_{j=1}^K \beta_j b_j(x)
$$

其中 $b_j(x)$ 可以是：

- $x^j$；
- $\sin(jx)$ 和 $\cos(jx)$；
- 分段线性函数；
- spline basis。

好处是：只要基函数事先确定，估计仍是线性模型。

**具体解决 / Konkrete Loesung.**  
先选择一组能表达问题结构的基函数，再把它们作为设计矩阵的列。后续 spline 和 additive model 就是这个思想的系统扩展。

**R 练习 / R-Uebung.**

```r
data(women)
h <- women[["height"]]
B <- cbind(1, h, h^2)
head(B)

coef_manual <- solve(t(B) %*% B, t(B) %*% women[["weight"]])
coef_lm <- coef(lm(weight ~ height + I(height^2), data = women))

cbind(manual = as.vector(coef_manual), lm = coef_lm)
```

## 6.6 多重检验与同时置信区间 / Multiples Testen und simultane Konfidenzbereiche

**遇到的问题 / Problem.**  
如果对很多点、很多参数或很多模型同时做检验，单个 $95\%$ 区间不再保证整体错误率为 $5\%$。

**可能的解决路径 / Moegliche Wege.**

- Bonferroni：简单保守。
- Scheffe：适合所有线性组合的同时推断。
- 椭球置信域：从参数向量整体看不确定性。
- 最大统计量方法：基于最大偏差控制同时覆盖。

**本课采用的方法及好处 / Methode und Vorteil.**  
Bonferroni 思想：

$$
\alpha_{single}=\frac{\alpha}{m}
$$

若有 $m$ 个区间，每个用更小显著性水平，整体错误率可被控制在 $\alpha$ 以内。

**具体解决 / Konkrete Loesung.**  
若只是少数预先指定的比较，可用 Bonferroni；若要对整条曲线做同时置信带，需要 Scheffe 或最大统计量思想。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)

alpha <- 0.05
m <- 3
confint(fit, level = 1 - alpha / m)
```

## 6.7 选择题 / Multiple-Choice-Fragen

1. 交互项 $xz$ 的作用是什么？  
   A. 允许 $x$ 的作用随 $z$ 改变  
   B. 删除 $x$ 和 $z$  
   C. 强制线性关系不存在  
   D. 只用于二元响应  
   答案 / Antwort: A

2. 多项式回归为什么仍属于线性模型？  
   A. 因为它关于参数 $\beta$ 线性  
   B. 因为 $x^2$ 是直线  
   C. 因为不需要误差项  
   D. 因为只能有一个参数  
   答案 / Antwort: A

3. $E[\log(Y)\mid x]\neq \log(E[Y\mid x])$ 提醒什么？  
   A. 变换响应变量后解释要谨慎  
   B. 不能使用任何变换  
   C. 所有模型都无效  
   D. $R^2$ 必须为 $1$  
   答案 / Antwort: A

4. 基函数模型的核心形式是什么？  
   A. $E(Y\mid x)=\beta_0+\sum_j\beta_jb_j(x)$  
   B. $Y=0$  
   C. $X=Y$  
   D. $\varepsilon=1$  
   答案 / Antwort: A

5. Bonferroni 方法主要解决什么问题？  
   A. 多重检验错误率控制  
   B. 只估计截距  
   C. 删除所有异常值  
   D. logistic 回归分类阈值  
   答案 / Antwort: A

## 6.8 本章总结 / Zusammenfassung

中文总结：

- 连续变量之间可以有交互，交互表示一个变量的作用依赖另一个变量。
- 非线性关系可以通过变换、多项式、三角多项式和样条进入线性模型。
- 这些方法的统一语言是基函数 $b_j(x)$。
- 模型可以关于 $x$ 非线性，但只要关于参数 $\beta$ 线性，仍能用线性模型估计。
- 多重检验需要控制整体错误率，不能把许多单独区间当成一个整体区间。

Deutsche Zusammenfassung:

- Interaktionen bedeuten, dass ein Effekt vom Wert einer anderen Variablen abhaengt.
- Nichtlineare Effekte koennen durch Transformationen, Polynome und Basisfunktionen modelliert werden.
- Das Modell bleibt linear in den Parametern.
- Bei vielen Tests oder Intervallen braucht man simultane Verfahren.

## 6.9 专业德语单词汇总 / Fachwortschatz

| Deutsch | 中文 | 说明 |
|---|---|---|
| metrische Einflussgroesse | 连续解释变量 | quantitative covariate |
| Interaktion | 交互 | $xz$ |
| Effektmodifikation | 效应修饰 | 作用随另一变量变化 |
| Transformation | 变量变换 | $\log(x)$, $1/x$ |
| polynomiale Regression | 多项式回归 | $x,x^2,x^3$ |
| trigonometrisches Polynom | 三角多项式 | sine/cosine basis |
| Basisfunktion | 基函数 | $b_j(x)$ |
| Spline | 样条 | 分段平滑函数 |
| simultaner Konfidenzbereich | 同时置信区间/带 | family-wise control |
| Bonferroni | Bonferroni 校正 | $\alpha/m$ |
| Konfidenzellipsoid | 置信椭球 | 参数向量不确定性 |
| Scheffe | Scheffe 方法 | 所有线性组合 |
