# 第 2 章学生笔记：简单线性回归模型 / Kapitel 2: Das einfache lineare Regressionsmodell

材料 / Grundlage:

- 讲义 / Folien: `LinearModellingPPT.pdf`, Kapitel 2 `Das einfache lineare Regressionsmodell`, slides 17-50
- 主教材 / Hauptlehrbuch: Fahrmeir, Kneib, Lang, Marx, `Regression: Models, Methods and Applications`, Chapter 2.2.1 `Simple Linear Regression Model`, Chapter 3.1-3.3 中与经典线性模型、最小二乘、推断和预测相关的内容
- R 补充 / R-Ergaenzung: Faraway, `Linear Models with R`, Chapter 2 `Estimation`, Chapter 3 `Inference`

## 覆盖审核 / Abdeckungspruefung

| 来源 / Quelle | 内容点 / Inhaltspunkt | 笔记位置 / Wo |
|---|---|---|
| 讲义 2 | 模型 $Y_i=\beta_0+\beta_1x_i+\varepsilon_i$ 与基本假设 | 2.1 |
| 讲义 2 | absolute distance vs quadratic distance | 2.2 |
| 讲义 2 | KQ-Methode / 最小二乘估计 | 2.2 |
| 讲义 2 | 存在唯一性、估计公式、正规方程 | 2.3 |
| 讲义 2 | 参数解释：$\beta_0,\beta_1,\sigma$ | 2.4 |
| 讲义 2 | 无偏性、方差、正态假设下 LS = ML | 2.5 |
| 讲义 2 | $\widehat{\sigma}^2$ 与 $n-2$ 自由度 | 2.6 |
| 讲义 2 | $t$ 分布、$\beta_0,\beta_1$ 置信区间 | 2.7 |
| 讲义 2 | Wenigumstadt 例子中的解释重点和潜在问题 | 2.8 |
| 讲义 2 | 平方和分解、自由度、$R^2$ | 2.9 |
| 讲义 2 | 预测值、均值置信区间、预测区间 | 2.10 |
| 讲义 2 | 变换、$E[g(Y)]\neq g[E(Y)]$、线性近似与外推风险 | 2.11 |
| 主教材 Ch. 2.2.1 | 简单线性回归标准模型、Munich rent 例子、变量变换 | 2.1, 2.11 |
| 主教材 Ch. 3.1 | 经典线性模型假设、误差与残差、同方差/独立/正态性 | 2.1, 2.4, 2.12 |
| 主教材 Ch. 3.2 | LS 估计、误差方差估计、估计量性质 | 2.2-2.6 |
| 主教材 Ch. 3.3 | 参数置信区间与预测区间 | 2.7, 2.10 |

审核结论 / Ergebnis: 第 2 讲 slides 17-50 的知识点已覆盖；主教材中与本章直接相关的简单线性回归、最小二乘估计、参数推断和预测区间也已纳入。多元线性模型、系统模型选择和完整诊断属于后续讲义章节，这里只保留第 2 章已经需要的部分。

**本章逻辑链 / Logikkette.**
先设定 $Y_i=\beta_0+\beta_1x_i+\varepsilon_i$，再用最小二乘估计 $\beta_0,\beta_1$；在零均值、同方差、独立和正态误差假设下，继续得到无偏性、方差、$t$ 区间、预测区间；最后用 $R^2$ 和残差图检查模型是否合理。

## 2.1 模型设定 / Modellspezifikation

**遇到的问题 / Problem.**
我们有一个连续响应变量 $Y$ 和一个解释变量 $x$。散点图可能显示二者大致线性相关，但每个点不会完全落在一条直线上。问题是：如何同时表达平均趋势和随机波动？

Deutsch: Wir wollen eine Zielgroesse $Y$ durch eine Einflussgroesse $x$ erklaeren. Die Beziehung ist ungefaehr linear, aber zufaellig gestreut.

**可能的解决路径 / Moegliche Wege.**

- 只用总体均值 $\bar{Y}$：简单，但完全忽略 $x$。
- 使用任意曲线 $f(x)$：灵活，但第 2 章还无法系统推断。
- 使用简单线性回归：用直线表示条件均值，用误差项表示随机偏差。

**本课采用的方法及好处 / Methode und Vorteil.**
讲义采用简单线性回归模型：

$$
Y_i=\beta_0+\beta_1x_i+\varepsilon_i,\quad i=1,\dots,n
$$

基本假设：

$$
E(\varepsilon_i)=0
$$

$$
V(\varepsilon_i)=\sigma^2
$$

$$
\varepsilon_1,\dots,\varepsilon_n\ \text{sind stochastisch unabhaengig}
$$

做精确的 $t$ 推断时再加正态性：

$$
\varepsilon_i\sim N(0,\sigma^2)
$$

好处是 $\beta_0$、$\beta_1$ 和 $\sigma$ 都有明确解释；估计、推断和预测可以放在同一框架中。

**具体解决 / Konkrete Loesung.**
模型的条件均值是：

$$
E(Y_i\mid x_i)=\beta_0+\beta_1x_i
$$

这里先把 $x_i$ 看作固定已知值 / feste bekannte Einflussgroesse；随机性主要来自 $\varepsilon_i$。

符号含义:

- $Y_i$: 目标变量 / Zielgroesse，随机。
- $x_i$: 解释变量 / Einflussgroesse，本章暂作固定已知。
- $\varepsilon_i$: 随机误差 / Zufallsfehler，不可观测。
- $\beta_0,\beta_1,\sigma^2$: 未知参数 / unbekannte Parameter。
- $\widehat{\varepsilon}_i$: 残差 / Residuum，是误差的样本近似，不等于真实误差。

**R 练习 / R-Uebung.**

```r
data(cars)
fit <- lm(dist ~ speed, data = cars)
summary(fit)

plot(cars[["speed"]], cars[["dist"]],
     xlab = "speed",
     ylab = "stopping distance",
     main = "Simple linear regression")
abline(fit, col = "red", lwd = 2)
```

练习 / Uebung: 在这个模型中，$Y_i$、$x_i$ 和 $\varepsilon_i$ 分别对应什么？
答案 / Antwort: $Y_i$ 是刹车距离，$x_i$ 是速度，$\varepsilon_i$ 是直线趋势无法解释的随机偏差。

## 2.2 为什么选择最小二乘 / Warum kleinste Quadrate?

**遇到的问题 / Problem.**
散点图上可以画无数条直线。我们必须定义什么叫“最好”的直线。

Deutsch: Viele Geraden sind moeglich; wir brauchen ein objektives Kriterium.

**可能的解决路径 / Moegliche Wege.**

- 最小化绝对偏差：$\sum_{i=1}^n |Y_i-\beta_0-\beta_1x_i|$。它较稳健，但经典推导不如最小二乘直接。
- 最小化平方偏差：$\sum_{i=1}^n (Y_i-\beta_0-\beta_1x_i)^2$。它有闭式解，几何性质好，且正态误差下等价于最大似然。
- 使用稳健回归或分位数回归：适合离群点或分布不对称，属于后续扩展。

**本课采用的方法及好处 / Methode und Vorteil.**
讲义采用 KQ-Methode / Methode der kleinsten Quadrate:

$$
(\widehat{\beta}_0,\widehat{\beta}_1)
=
\arg\min_{\beta_0,\beta_1}
\sum_{i=1}^n (Y_i-\beta_0-\beta_1x_i)^2
$$

残差定义为：

$$
\widehat{\varepsilon}_i=Y_i-\widehat{\beta}_0-\widehat{\beta}_1x_i
$$

好处是：计算方便、有显式公式、可推广到多元线性模型，并且为平方和分解和 $R^2$ 打基础。

**具体解决 / Konkrete Loesung.**
对每一条候选直线，先计算所有残差，再计算残差平方和。残差平方和最小的直线，就是本章的拟合直线。

**R 练习 / R-Uebung.**

```r
data(cars)
x <- cars[["speed"]]
y <- cars[["dist"]]

fit_lm <- lm(y ~ x)
coef(fit_lm)

sum(residuals(fit_lm)^2)
sum(abs(residuals(fit_lm)))
```

练习 / Uebung: 为什么最小二乘对离群点敏感？
答案 / Antwort: 因为残差被平方，大残差会被放大。

## 2.3 存在唯一性、公式和正规方程 / Existenz, Formeln und Normalgleichungen

**遇到的问题 / Problem.**
最小化准则给出了目标函数，但我们还要知道估计量是否存在、是否唯一、如何计算。

Deutsch: Nach der Wahl des Kriteriums muessen Existenz, Eindeutigkeit und Berechnung geklaert werden.

**可能的解决路径 / Moegliche Wege.**

- 数值优化：通用，但对简单线性回归没有必要。
- 对目标函数求偏导：得到显式公式。
- 矩阵法：适合推广到多元线性回归。

**本课采用的方法及好处 / Methode und Vorteil.**
如果 $x$ 有变异，即：

$$
\sum_{i=1}^n (x_i-\bar{x})^2\neq 0
$$

则最小二乘估计存在且唯一：

$$
\widehat{\beta}_1=
\frac{\sum_{i=1}^n (x_i-\bar{x})(Y_i-\bar{Y})}
{\sum_{i=1}^n (x_i-\bar{x})^2}
$$

$$
\widehat{\beta}_0=\bar{Y}-\widehat{\beta}_1\bar{x}
$$

正规方程 / Normalgleichungen:

$$
\sum_{i=1}^n \widehat{\varepsilon}_i=0
$$

$$
\sum_{i=1}^n x_i\widehat{\varepsilon}_i=0
$$

这些公式的好处是：我们不仅能让软件给结果，还能理解拟合线为什么经过数据中心并平衡残差。

**具体解决 / Konkrete Loesung.**
先检查 $x_i$ 是否真的有变化。如果所有 $x_i$ 都相同，分母为零，斜率无法估计。若 $x_i$ 有变化，就先算 $\bar{x}$ 和 $\bar{Y}$，再用上面的公式计算 $\widehat{\beta}_1$ 和 $\widehat{\beta}_0$。

**R 练习 / R-Uebung.**

```r
data(cars)
x <- cars[["speed"]]
y <- cars[["dist"]]

b1 <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
b0 <- mean(y) - b1 * mean(x)
c(b0 = b0, b1 = b1)

coef(lm(dist ~ speed, data = cars))

res <- y - b0 - b1 * x
c(sum_residuals = sum(res),
  sum_x_residuals = sum(x * res))
```

练习 / Uebung: 如果所有 $x_i$ 都相同，为什么不能估计斜率？
答案 / Antwort: 因为 $\sum (x_i-\bar{x})^2=0$，没有关于 $x$ 变化的信息。

## 2.4 参数解释 / Interpretation der Parameter

**遇到的问题 / Problem.**
软件输出中有截距、斜率和残差标准误。如果只看显著性，不解释单位和语境，很容易误读。

Deutsch: Koeffizienten muessen inhaltlich und einheitenbezogen interpretiert werden.

**可能的解决路径 / Moegliche Wege.**

- 只看相关系数：只能看线性关联强弱，不能直接解释单位变化。
- 看斜率 $\beta_1$：解释 $x$ 增加一单位时 $Y$ 的平均变化。
- 看截距 $\beta_0$：只有 $x=0$ 有意义且在合理范围内时才解释。
- 看 $\sigma$：衡量点围绕回归线的典型波动。

**本课采用的方法及好处 / Methode und Vorteil.**
用条件均值解释：

$$
E(Y\mid x+1)-E(Y\mid x)=\beta_1
$$

这比“$Y$ 一定增加 $\beta_1$”更准确，因为单个观测还包含 $\varepsilon$。

**具体解决 / Konkrete Loesung.**

- $\beta_1$: Steigungsparameter / 斜率参数，表示 $x$ 增加一个单位时 $Y$ 的条件均值平均变化 $\beta_1$ 个单位。
- $\beta_0$: Achsenabschnitt / 截距，表示 $x=0$ 时的条件均值；如果 $x=0$ 不在合理范围内，就不要过度解释。
- $\sigma$: Standardabweichung des Stoerterms / 误差标准差，表示观测围绕回归线的典型距离。
- $\widehat{\varepsilon}_i$: 残差，是用估计模型得到的误差近似；真实误差 $\varepsilon_i$ 看不到。

常见误解：$\widehat{\beta}_1=3.9$ 不表示每个个体都增加 $3.9$，而表示条件均值平均增加 $3.9$。

**R 练习 / R-Uebung.**

```r
data(cars)
fit <- lm(dist ~ speed, data = cars)
coef(fit)
sigma(fit)

new_x <- data.frame(speed = c(10, 11))
preds <- predict(fit, newdata = new_x)
diff(preds)
coef(fit)[2]
```

练习 / Uebung: `speed` 从 $10$ 到 $11$ 时，预测均值变化应接近哪个数？
答案 / Antwort: 接近 $\widehat{\beta}_1$。

## 2.5 估计量性质：无偏性、方差与最大似然 / Eigenschaften: Erwartungstreue, Varianz, ML

**遇到的问题 / Problem.**
点估计只是一个数。我们还要知道它是否系统偏离真值，以及换样本时会有多大波动。

Deutsch: Punktschaetzer muessen hinsichtlich Bias und Varianz beurteilt werden.

**可能的解决路径 / Moegliche Wege.**

- 只报告 $\widehat{\beta}_0,\widehat{\beta}_1$：不够。
- 做 bootstrap：可行，但不是讲义当前主线。
- 通过模型假设推导估计量性质：能得到标准误、置信区间和检验。

**本课采用的方法及好处 / Methode und Vorteil.**
在 $E(\varepsilon_i)=0$ 下：

$$
E(\widehat{\beta}_0)=\beta_0,\quad E(\widehat{\beta}_1)=\beta_1
$$

在独立同方差误差下：

$$
V(\widehat{\beta}_1)=
\frac{\sigma^2}{\sum_{i=1}^n (x_i-\bar{x})^2}
$$

$$
V(\widehat{\beta}_0)=
\sigma^2\left(
\frac{1}{n}
+
\frac{\bar{x}^2}{\sum_{i=1}^n (x_i-\bar{x})^2}
\right)
$$

正态误差下，最小二乘估计同时是最大似然估计。因为：

$$
Y_i\sim N(\beta_0+\beta_1x_i,\sigma^2)
$$

其 log-likelihood 中与 $\beta_0,\beta_1$ 有关的部分等价于最小化：

$$
\sum_{i=1}^n (Y_i-\beta_0-\beta_1x_i)^2
$$

好处是：我们可以从同一模型推出标准误、置信区间、检验和预测区间。

**具体解决 / Konkrete Loesung.**
先确认误差期望为零，得到无偏性；再加独立同方差，得到方差公式；最后加正态性，把最小二乘与最大似然联系起来。注意：正态 ML 下 $\sigma^2$ 的估计分母是 $n$，但它有偏；讲义后面使用 $n-2$ 的无偏估计。

**R 练习 / R-Uebung.**

```r
set.seed(1)

simulate_slopes <- function(x, B = 1000) {
  beta0 <- 2
  beta1 <- 3
  sigma <- 2
  slopes <- numeric(B)
  for (b in seq_len(B)) {
    y <- beta0 + beta1 * x + rnorm(length(x), sd = sigma)
    slopes[b] <- coef(lm(y ~ x))[2]
  }
  sd(slopes)
}

x_narrow <- seq(0, 1, length.out = 30)
x_wide <- seq(0, 10, length.out = 30)

c(narrow = simulate_slopes(x_narrow),
  wide = simulate_slopes(x_wide))
```

练习 / Uebung: 哪个斜率模拟标准差通常更小？
答案 / Antwort: `wide` 更小，因为 $\sum (x_i-\bar{x})^2$ 更大。

## 2.6 误差方差估计与自由度 / Schaetzung von $\sigma^2$ und Freiheitsgrade

**遇到的问题 / Problem.**
$\sigma^2$ 是误差方差，但真实误差 $\varepsilon_i$ 看不到，只能用残差 $\widehat{\varepsilon}_i$ 近似。问题是分母该用 $n$ 还是 $n-2$？

Deutsch: Die Fehler sind unbekannt; deshalb verwenden wir Residuen zur Schaetzung von $\sigma^2$.

**可能的解决路径 / Moegliche Wege.**

- 使用 $\frac{1}{n}\sum\widehat{\varepsilon}_i^2$：这是正态模型下的 ML 形式，但有偏。
- 使用 $\frac{1}{n-2}\sum\widehat{\varepsilon}_i^2$：调整估计两个参数造成的自由度损失，是无偏估计。

**本课采用的方法及好处 / Methode und Vorteil.**

$$
\widehat{\sigma}^2=
\frac{1}{n-2}
\sum_{i=1}^n \widehat{\varepsilon}_i^2
$$

好处是：

- $E(\widehat{\sigma}^2)=\sigma^2$，即无偏；
- 与 $t$ 区间和检验直接衔接；
- 自由度解释清楚：简单线性回归估计了 $\beta_0$ 和 $\beta_1$ 两个参数，所以剩余自由度为 $n-2$。

**具体解决 / Konkrete Loesung.**
先拟合模型并得到残差，再计算残差平方和 $SSE=\sum\widehat{\varepsilon}_i^2$，最后除以 $n-2$。

**R 练习 / R-Uebung.**

```r
data(cars)
fit <- lm(dist ~ speed, data = cars)
res <- residuals(fit)
n <- length(res)

sigma2_unbiased <- sum(res^2) / (n - 2)
sigma2_ml <- sum(res^2) / n

c(unbiased = sigma2_unbiased,
  ml = sigma2_ml,
  sigma_from_lm = sigma(fit)^2)
```

练习 / Uebung: 哪个值更接近 `sigma(fit)^2`？
答案 / Antwort: `unbiased`，因为 `lm()` 使用 $n-2$。

## 2.7 置信区间 / Konfidenzintervalle

**遇到的问题 / Problem.**
我们想知道 $\beta_1$ 是否可能为 $0$，或者 $\beta_0,\beta_1$ 的合理范围是多少。点估计无法回答不确定性。

Deutsch: Punktschaetzungen reichen nicht; wir brauchen Unsicherheitsintervalle.

**可能的解决路径 / Moegliche Wege.**

- 只报告点估计：没有不确定性。
- 若 $\sigma$ 已知，用正态分布：实际很少直接可用。
- 用 $\widehat{\sigma}$ 代替 $\sigma$，并使用 $t$ 分布。

**本课采用的方法及好处 / Methode und Vorteil.**
对斜率：

$$
\frac{\widehat{\beta}_1-\beta_1}
{\widehat{\sigma}_{\widehat{\beta}_1}}
\sim t_{n-2}
$$

其中：

$$
\widehat{\sigma}_{\widehat{\beta}_1}
=
\frac{\widehat{\sigma}}
{\sqrt{\sum_{i=1}^n (x_i-\bar{x})^2}}
$$

因此 $1-\alpha$ 置信区间为：

$$
\widehat{\beta}_1
\pm
t_{n-2,1-\alpha/2}
\widehat{\sigma}_{\widehat{\beta}_1}
$$

对截距：

$$
\widehat{\sigma}_{\widehat{\beta}_0}
=
\widehat{\sigma}
\sqrt{
\frac{1}{n}
+
\frac{\bar{x}^2}{\sum_{i=1}^n (x_i-\bar{x})^2}
}
$$

$$
\widehat{\beta}_0
\pm
t_{n-2,1-\alpha/2}
\widehat{\sigma}_{\widehat{\beta}_0}
$$

好处是：区间直接表达参数估计的不确定性，比只看点估计更诚实。

**具体解决 / Konkrete Loesung.**
拟合 `lm()` 后读取 `confint(fit)`。如果 $\beta_1$ 的置信区间不含 $0$，说明在该模型和假设下，线性斜率有统计证据不为 $0$。

**R 练习 / R-Uebung.**

```r
data(cars)
fit <- lm(dist ~ speed, data = cars)
confint(fit)
summary(fit)[["coefficients"]]
```

练习 / Uebung: `speed` 的置信区间是否包含 $0$？
答案 / Antwort: 通常不包含，说明该线性模型下速度和平均刹车距离之间有明显线性关联。

## 2.8 讲义例子：Wenigumstadt 骨小梁厚度 / Beispiel: Wenigumstadt

**遇到的问题 / Problem.**
讲义中的 Wenigumstadt 例子研究死亡年龄与第四腰椎骨小梁厚度之间的关系。重点不是只画回归线，而是解释斜率、误差规模、拟合优度和潜在数据问题。

Deutsch: Im Beispiel stehen Interpretation der Steigung, Streuung und Anpassungsguete im Zentrum.

**可能的解决路径 / Moegliche Wege.**

- 只报告散点图：直观但不够量化。
- 只报告斜率：忽略残差波动。
- 同时报告 $\widehat{\beta}_1$、$\widehat{\sigma}$、$R^2$ 和假设问题：更完整。

**本课采用的方法及好处 / Methode und Vorteil.**
讲义特别提醒潜在问题：

- Messfehler / 测量误差；
- sehr kleine Werte / 很小的观测值；
- Zufallige $x_i$ / $x_i$ 可能也有随机性，而本章模型暂时把 $x_i$ 当固定值。

好处是：学生不会把 `lm()` 输出当成自动真理，而会同时检查数据质量和模型假设。

**具体解决 / Konkrete Loesung.**
读这种例子时按顺序做：确认 $Y$ 和 $x$ 的含义；解释 $\widehat{\beta}_1$ 的单位变化；报告 $\widehat{\sigma}$ 表示的典型残差规模；用 $R^2$ 说明解释变异比例；最后写下测量误差和假设疑点。

**R 练习 / R-Uebung.** 用内置 `women` 数据模拟同类“连续 $Y$ 对连续 $x$”解释过程。

```r
data(women)
fit <- lm(weight ~ height, data = women)
summary(fit)

c(slope = coef(fit)[2],
  sigma = sigma(fit),
  r_squared = summary(fit)[["r.squared"]])
```

练习 / Uebung: 斜率为正是否表示每个身高更高的人都一定更重？
答案 / Antwort: 不是。斜率描述的是平均趋势，不是个体层面的确定关系。

## 2.9 平方和分解与 $R^2$ / Quadratsummenzerlegung und Bestimmtheitsmass

**遇到的问题 / Problem.**
拟合出回归线后，我们想知道模型解释了多少 $Y$ 的变异。单看斜率不够，因为斜率受单位影响。

Deutsch: Wir brauchen ein Mass fuer die Anpassungsguete.

**可能的解决路径 / Moegliche Wege.**

- 看残差图：能诊断问题，但不是单一数值指标。
- 看 $SSE$：有单位，跨问题比较不直观。
- 看 $R^2$：给出解释变异占总变异的比例。

**本课采用的方法及好处 / Methode und Vorteil.**
平方和分解：

$$
\sum_{i=1}^n (Y_i-\bar{Y})^2
=
\sum_{i=1}^n (Y_i-\widehat{Y}_i)^2
+
\sum_{i=1}^n (\widehat{Y}_i-\bar{Y})^2
$$

即：

$$
SST=SSE+SSM
$$

其中：

$$
SST=\sum_{i=1}^n (Y_i-\bar{Y})^2
$$

$$
SSE=\sum_{i=1}^n (Y_i-\widehat{Y}_i)^2
$$

$$
SSM=\sum_{i=1}^n (\widehat{Y}_i-\bar{Y})^2
$$

判定系数 / Bestimmtheitsmass:

$$
R^2=\frac{SSM}{SST}=1-\frac{SSE}{SST}
$$

简单线性回归且含截距时：

$$
R^2=r_{xY}^2
$$

好处是：$R^2$ 把拟合优度变成一个比例语言。

**具体解决 / Konkrete Loesung.**
先计算总变异 $SST$，再计算残差变异 $SSE$，然后用 $1-SSE/SST$ 得到 $R^2$。注意：高 $R^2$ 不等于因果关系；$R^2=0$ 也只表示没有线性解释，不表示完全没有关系。

**R 练习 / R-Uebung.**

```r
data(cars)
fit <- lm(dist ~ speed, data = cars)
y <- cars[["dist"]]
yhat <- fitted(fit)

SST <- sum((y - mean(y))^2)
SSE <- sum((y - yhat)^2)
SSM <- sum((yhat - mean(y))^2)

c(SST = SST, SSE = SSE, SSM = SSM, SSE_plus_SSM = SSE + SSM)
c(R2_manual = SSM / SST,
  R2_summary = summary(fit)[["r.squared"]],
  cor_squared = cor(cars[["speed"]], cars[["dist"]])^2)
```

练习 / Uebung: 简单线性回归含截距时，$R^2$ 与相关系数有什么关系？
答案 / Antwort: $R^2=r_{xY}^2$。

## 2.10 预测与预测区间 / Prognose und Prognoseintervall

**遇到的问题 / Problem.**
给定新解释变量值 $x_{n+1}$，我们想预测未来的 $Y_{n+1}$。点预测容易给出，但单个新观测还包含随机误差，因此需要区间。

Deutsch: Eine Prognose braucht sowohl einen Punktwert als auch eine Unsicherheitsaussage.

**可能的解决路径 / Moegliche Wege.**

- 只给 $\widehat{Y}_{n+1}$：没有不确定性。
- 给均值置信区间：回答 $E(Y\mid x_{n+1})$ 的范围。
- 给预测区间：回答未来单个 $Y_{n+1}$ 的可能范围。

**本课采用的方法及好处 / Methode und Vorteil.**
点预测：

$$
\widehat{Y}_{n+1}
=
\widehat{\beta}_0+\widehat{\beta}_1x_{n+1}
$$

未来单个观测的预测区间为：

$$
\widehat{Y}_{n+1}
\pm
t_{n-2,1-\alpha/2}
\widehat{\sigma}
\sqrt{
1+\frac{1}{n}
+
\frac{(x_{n+1}-\bar{x})^2}{\sum_{i=1}^n (x_i-\bar{x})^2}
}
$$

好处是：预测区间同时包含均值估计不确定性和未来个体误差，因此比均值置信区间更宽。

**具体解决 / Konkrete Loesung.**
先问问题对象是均值 $E(Y\mid x_{n+1})$ 还是未来个体 $Y_{n+1}$。如果是未来个体，用 `interval = "prediction"`；如果是均值，用 `interval = "confidence"`。

**R 练习 / R-Uebung.**

```r
data(cars)
fit <- lm(dist ~ speed, data = cars)
new_speed <- data.frame(speed = c(10, 20))

predict(fit, newdata = new_speed, interval = "confidence")
predict(fit, newdata = new_speed, interval = "prediction")
```

练习 / Uebung: 哪个区间更宽？为什么？
答案 / Antwort: `prediction` 更宽，因为它还包含未来单个观测的随机误差。

## 2.11 变换、线性近似与外推风险 / Transformation, lineare Naeherung und Extrapolation

**遇到的问题 / Problem.**
真实关系可能不是直线。主教材的 Munich rent 例子和讲义的二次函数近似例子都提醒我们：线性模型有时只是局部近似。

Deutsch: Ein lineares Modell kann lokal sinnvoll sein, aber bei Nichtlinearitaet oder Extrapolation problematisch werden.

**可能的解决路径 / Moegliche Wege.**

- 继续用直线：简单，但可能错设。
- 变换解释变量，如 $x=1/\text{area}$ 或 $x=\log(z)$。
- 变换响应变量，如 $\log(Y)$：可处理乘性误差或右偏分布，但解释会变化。
- 加入多项式，如 $x^2$。
- 后续使用 splines 或非参数回归。

**本课采用的方法及好处 / Methode und Vorteil.**
第 2 章先用线性模型建立基础。优点是简单、可解释、可推断；缺点是必须检查线性关系是否合理，特别是不能随便外推。

讲义提醒：

$$
E[g(Y)]\neq g[E(Y)]
$$

因此变换 $Y$ 后，参数解释和预测解释都会改变。

**具体解决 / Konkrete Loesung.**
先画散点图和残差图。如果看到弯曲，尝试加入 $x^2$ 或做变量变换。若预测点远离观测范围，要明确标注外推风险。

例子：

$$
Y_i=\beta_0+\beta_1x_i+\beta_2x_i^2+\varepsilon_i
$$

这个模型关于 $x$ 是二次的，但关于参数 $\beta_0,\beta_1,\beta_2$ 仍然是线性的，所以仍属于线性模型框架。

**R 练习 / R-Uebung.**

```r
data(women)
fit_linear <- lm(weight ~ height, data = women)
fit_quad <- lm(weight ~ height + I(height^2), data = women)

plot(women[["height"]], women[["weight"]],
     xlab = "height",
     ylab = "weight",
     main = "Linear vs quadratic fit")
abline(fit_linear, col = "red", lwd = 2)

h <- seq(min(women[["height"]]), max(women[["height"]]), length.out = 100)
pred_quad <- predict(fit_quad, newdata = data.frame(height = h))
lines(h, pred_quad, col = "blue", lwd = 2)
legend("topleft",
       legend = c("linear", "quadratic"),
       col = c("red", "blue"),
       lwd = 2)
```

练习 / Uebung: 若残差图呈系统性弯曲，下一步可能做什么？
答案 / Antwort: 尝试变量变换、多项式项或后续的 spline 方法，而不是只看原线性模型的 $p$ 值。

## 2.12 模型假设检查的早期意识 / Fruehes Bewusstsein fuer Modellannahmen

**遇到的问题 / Problem.**
讲义第 2 章主要讲估计和推断，但主教材 Chapter 3.1 已经提醒：经典线性模型依赖误差均值、方差、相关性、正态性和线性结构等假设。若假设严重不成立，置信区间和预测区间会不可靠。

Deutsch: Die Inferenz des linearen Modells haengt von Modellannahmen ab.

**可能的解决路径 / Moegliche Wege.**

- 暂时忽略假设：能快速学习主线，但实践中危险。
- 用残差图提前观察：简单有效。
- 后续引入 WLS/GLS、稳健方法、变换、mixed model：用于修正具体问题。

**本课采用的方法及好处 / Methode und Vorteil.**
第 2 章先学最基本模型，但保留检查意识：

- 同方差 / Homoskedastizitaet: $V(\varepsilon_i)=\sigma^2$。
- 不相关 / unkorrelierte Fehler: $\operatorname{Cov}(\varepsilon_i,\varepsilon_j)=0$。
- 正态误差 / normalverteilte Fehler: $\varepsilon_i\sim N(0,\sigma^2)$。
- 线性均值 / Linearitaet: $E(Y_i\mid x_i)=\beta_0+\beta_1x_i$。

**具体解决 / Konkrete Loesung.**
拟合简单线性模型后，至少画残差对拟合值图和 QQ 图。若看到弯曲、扇形扩散、强离群点或明显非正态尾部，就把问题记录下来，并在后续章节用变换、加权最小二乘、相关误差模型、mixed model 或更灵活的非线性项处理。

**R 练习 / R-Uebung.**

```r
data(cars)
fit <- lm(dist ~ speed, data = cars)

par(mfrow = c(1, 2))
plot(fitted(fit), residuals(fit),
     xlab = "fitted values",
     ylab = "residuals")
abline(h = 0, col = "gray")
qqnorm(residuals(fit))
qqline(residuals(fit), col = "red")
par(mfrow = c(1, 1))
```

练习 / Uebung: 如果残差随 fitted values 增大而扩散，可能对应什么问题？
答案 / Antwort: 可能是异方差 / Heteroskedastizitaet，即 $V(\varepsilon_i)$ 不是常数。

## 2.13 选择题 / Multiple-Choice-Fragen

1. 简单线性回归中 $\beta_1$ 表示什么？
   A. $x$ 增加一单位时 $Y$ 条件均值的平均变化
   B. $Y$ 的总方差
   C. 样本量
   D. 残差平方和
   答案 / Antwort: A

2. 最小二乘法最小化什么？
   A. $\sum_{i=1}^n |Y_i-\beta_0-\beta_1x_i|$
   B. $\sum_{i=1}^n (Y_i-\beta_0-\beta_1x_i)^2$
   C. $\sum_{i=1}^n x_i$
   D. $\sum_{i=1}^n Y_i$
   答案 / Antwort: B

3. 为什么 $\widehat{\sigma}^2$ 的常用无偏估计分母是 $n-2$？
   A. 因为简单线性回归估计了 $\beta_0$ 和 $\beta_1$ 两个参数
   B. 因为样本量总是 $2$
   C. 因为 $R^2$ 有两个部分
   D. 因为误差没有方差
   答案 / Antwort: A

4. 预测区间通常比均值置信区间宽，原因是什么？
   A. 预测区间还包含未来单个观测的随机误差
   B. 预测区间不使用数据
   C. 预测区间要求 $R^2=1$
   D. 预测区间不使用 $\widehat{\sigma}$
   答案 / Antwort: A

5. 简单线性回归且含截距时，$R^2$ 等于什么？
   A. $r_{xY}^2$
   B. $\widehat{\beta}_0$
   C. $\widehat{\sigma}^2$
   D. $n-2$
   答案 / Antwort: A

6. 如果所有 $x_i$ 都相同，最直接的问题是什么？
   A. $\sum (x_i-\bar{x})^2=0$，斜率无法估计
   B. $Y_i$ 必定全相同
   C. $R^2$ 必定为 $1$
   D. $\sigma$ 必定为 $0$
   答案 / Antwort: A

7. $E[g(Y)]\neq g[E(Y)]$ 提醒我们什么？
   A. 变换响应变量后，解释和反变换预测要谨慎
   B. 所有线性模型都无效
   C. $R^2$ 必须为 $0$
   D. 最小二乘不能使用
   答案 / Antwort: A

8. 在简单线性回归中，如果 $R^2=0.95$，可以直接得出因果结论吗？
   A. 可以，因为 $R^2$ 很高
   B. 可以，因为模型一定正确
   C. 不可以，$R^2$ 只说明线性拟合解释了多少变异
   D. 不可以，因为 $R^2$ 总是无意义
   答案 / Antwort: C

## 2.14 本章总结 / Zusammenfassung

中文总结：

- 简单线性回归模型为 $Y_i=\beta_0+\beta_1x_i+\varepsilon_i$。
- $\beta_1$ 表示 $x$ 增加一个单位时 $Y$ 条件均值的平均变化；$\beta_0$ 只有在 $x=0$ 合理时才有直接解释；$\sigma$ 描述随机误差的典型规模。
- 最小二乘通过最小化 $\sum (Y_i-\beta_0-\beta_1x_i)^2$ 得到 $\widehat{\beta}_0,\widehat{\beta}_1$。
- 估计量存在唯一需要 $\sum (x_i-\bar{x})^2\neq 0$。
- 正规方程给出 $\sum\widehat{\varepsilon}_i=0$ 和 $\sum x_i\widehat{\varepsilon}_i=0$。
- 在 $E(\varepsilon_i)=0$ 下，$\widehat{\beta}_0,\widehat{\beta}_1$ 无偏；在独立同方差下可推导方差；在正态误差下 LS = ML。
- $\widehat{\sigma}^2=\frac{1}{n-2}\sum\widehat{\varepsilon}_i^2$ 是误差方差的无偏估计。
- $SST=SSE+SSM$，$R^2=1-\frac{SSE}{SST}$ 表示解释变异比例。
- 预测区间比均值置信区间宽，因为它包含未来单个观测的随机误差。
- 线性模型可能只是局部近似，变换、外推和残差诊断都必须谨慎。

Deutsche Zusammenfassung:

- Das einfache lineare Regressionsmodell lautet $Y_i=\beta_0+\beta_1x_i+\varepsilon_i$.
- $\beta_1$ ist der Steigungsparameter; $\beta_0$ ist der Achsenabschnitt; $\sigma$ beschreibt die typische Fehlerstreuung.
- Die KQ-Methode minimiert die Residuenquadratsumme $\sum (Y_i-\beta_0-\beta_1x_i)^2$.
- Eindeutigkeit verlangt Variation in $x$, also $\sum (x_i-\bar{x})^2\neq 0$.
- Die Normalgleichungen implizieren $\sum\widehat{\varepsilon}_i=0$ und $\sum x_i\widehat{\varepsilon}_i=0$.
- Unter $E(\varepsilon_i)=0$ sind die KQ-Schaetzer erwartungstreu; bei normalverteilten Fehlern gilt LS = ML.
- Das Bestimmtheitsmass $R^2$ misst den erklaerten Anteil der Gesamtstreuung.
- Prognoseintervalle sind breiter als Konfidenzintervalle fuer den Mittelwert.

## 2.15 专业德语单词汇总 / Fachwortschatz

| Deutsch | 中文 | 公式 / Hinweis |
|---|---|---|
| einfaches lineares Regressionsmodell | 简单线性回归模型 | $Y_i=\beta_0+\beta_1x_i+\varepsilon_i$ |
| Zielgroesse | 目标变量 | $Y_i$ |
| Einflussgroesse | 解释变量 | $x_i$ |
| feste bekannte Einflussgroesse | 固定已知解释变量 | 本章暂设 $x_i$ 固定 |
| Achsenabschnitt | 截距 | $\beta_0$ |
| Steigungsparameter | 斜率参数 | $\beta_1$ |
| Stoerterm / Fehlerterm | 误差项 | $\varepsilon_i$ |
| Residuum | 残差 | $\widehat{\varepsilon}_i=Y_i-\widehat{Y}_i$ |
| Methode der kleinsten Quadrate | 最小二乘法 | $\min\sum (Y_i-\beta_0-\beta_1x_i)^2$ |
| KQ-Schaetzer | 最小二乘估计量 | $\widehat{\beta}_0,\widehat{\beta}_1$ |
| Normalgleichungen | 正规方程 | $\sum\widehat{\varepsilon}_i=0$, $\sum x_i\widehat{\varepsilon}_i=0$ |
| erwartungstreu | 无偏 | $E(\widehat{\beta}_j)=\beta_j$ |
| Varianz des Schaetzers | 估计量方差 | $V(\widehat{\beta}_1)=\sigma^2/\sum (x_i-\bar{x})^2$ |
| Maximum-Likelihood-Schaetzer | 最大似然估计 | 正态误差下 LS = ML |
| Freiheitsgrade | 自由度 | 简单线性回归中 $SSE$ 的自由度为 $n-2$ |
| Standardfehler | 标准误 | $\widehat{\sigma}_{\widehat{\beta}_j}$ |
| Konfidenzintervall | 置信区间 | $\widehat{\beta}_j\pm t\widehat{\sigma}_{\widehat{\beta}_j}$ |
| Prognosewert | 预测值 | $\widehat{Y}_{n+1}=\widehat{\beta}_0+\widehat{\beta}_1x_{n+1}$ |
| Prognoseintervall | 预测区间 | 包含未来误差 $\varepsilon_{n+1}$ |
| Quadratsummenzerlegung | 平方和分解 | $SST=SSE+SSM$ |
| Bestimmtheitsmass | 判定系数 | $R^2=1-\frac{SSE}{SST}$ |
| Homoskedastizitaet | 同方差 | $V(\varepsilon_i)=\sigma^2$ |
| Heteroskedastizitaet | 异方差 | $V(\varepsilon_i)$ 不是常数 |
| Transformation | 变换 | 如 $x=\log(z)$ 或 $\log(Y)$ |
| Extrapolation | 外推 | 观测范围外预测需谨慎 |
