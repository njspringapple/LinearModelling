# 第 3 章学生笔记：多元线性回归模型 / Kapitel 3: Das multiple lineare Regressionsmodell

材料 / Grundlage:

- 讲义 / Folien: `chapter_pdfs/LinearModellingPPT/03_Das_multiple_lineare_Regressionsmodell.pdf`, slides 51-66
- 主教材 / Hauptlehrbuch: Fahrmeir, Kneib, Lang, Marx, `Regression: Models, Methods and Applications`, Chapter 2.2.2, Chapter 3.1, Chapter 3.2, Appendix A
- 辅助教材 / Hilfslehrbuch: Faraway, `Linear Models with R`, Chapter 2 `Estimation`
- 扩展教材 / Erweiterung: Faraway, `Extending The Linear Model With R`, Chapter 1 only as a short review, not the main source for this chapter

## 覆盖审核 / Abdeckungspruefung

| 来源 / Quelle | 内容点 / Inhaltspunkt | 笔记位置 / Wo |
|---|---|---|
| 讲义 3 | 多元线性模型 $Y_i=\beta_0+\beta_1x_{i1}+\dots+\beta_px_{ip}+\varepsilon_i$ | 3.1 |
| 讲义 3 | 矩阵形式 $Y=X\beta+\varepsilon$ 与设计矩阵 | 3.1, 3.2 |
| 讲义 3 | 模型假设 $E(\varepsilon)=0$, $V(\varepsilon)=\sigma^2I$, 正态误差 | 3.2 |
| 讲义 3 | 参数解释：固定其他变量时解释 $\beta_k$ | 3.3 |
| 讲义 3 | KQ-Schaetzer / 最小二乘估计 | 3.4 |
| 讲义 3 | $X^TX$ 可逆、正规方程 $X^T\widehat{\varepsilon}=0$ | 3.4, 3.5 |
| 讲义 3 | Produktsummenmatrix $X^TX$ | 3.5 |
| 讲义 3 | $\widehat{\beta}$ 的无偏性、方差和正态分布 | 3.6 |
| 讲义 3 | Hat-Matrix $P$ 与 Residualmatrix $Q$ | 3.7 |
| 讲义 3 | 投影几何、$P^T=P$, $P^2=P$, $Q^T=Q$, $Q^2=Q$, $PQ=0$ | 3.7, 3.8 |
| 讲义 3 | $V(\widehat{Y})=\sigma^2P$, $V(\widehat{\varepsilon})=\sigma^2Q$ | 3.8 |
| 讲义 3 | $\widehat{\sigma}^2=\frac{1}{n-(p+1)}\widehat{\varepsilon}^T\widehat{\varepsilon}$ | 3.9 |
| 讲义 3 | Lesetest 例子和多协变量解释 | 3.10 |
| 主教材 Ch. 2.2.2 | 多元线性回归的模型族位置 | 3.1, 3.3 |
| 主教材 Ch. 3.1 | full rank、误差与残差、条件于 $X$ 的模型假设 | 3.2, 3.7 |
| 主教材 Ch. 3.2 | LS 估计、误差方差估计、估计量性质 | 3.4, 3.6, 3.9 |
| 主教材 Appendix A | 矩阵、转置、逆、秩、迹、投影的基础 | 3.5, 3.8 |
| Faraway Ch. 2 | 矩阵表示、模型空间、hat matrix、Gauss-Markov、identifiability 和 R 提取 | 3.4, 3.6, 3.11 |

审核结论 / Ergebnis: 第 3 讲 slides 51-66 的知识点已逐项覆盖；主教材中与本讲直接相关的模型定义、参数估计、矩阵代数和估计量性质已纳入；Faraway 第 2 章作为 R 与几何直觉的辅助材料补充。后续统计推断、平方和分解和一般线性假设主要属于第 4 讲，这里只在必要处预告。

## 3.1 从简单到多元：为什么需要矩阵形式 / Vom einfachen zum multiplen Modell

**遇到的问题 / Problem.**  
第 2 章只有一个解释变量 $x$。实际问题中，目标变量 $Y$ 往往同时受多个因素影响。例如阅读测试错误数可能和性别、年级、阅读时间、电视时间等变量都有关系。如果一次只放一个变量，容易把其他变量的影响混到当前变量上。

Auf Deutsch: In vielen Anwendungen haengt die Zielgroesse $Y$ gleichzeitig von mehreren Einflussgroessen ab.

**可能的解决路径 / Moegliche Wege.**

- 分别做很多个简单回归：容易遗漏混杂因素。
- 只保留一个“最重要”的解释变量：模型太粗。
- 建立多元线性回归：同时控制多个变量，解释每个变量在其他变量固定时的边际作用。

**本课采用的方法及好处 / Methode und Vorteil.**  
讲义采用：

$$
Y_i=\beta_0+\beta_1x_{i1}+\beta_2x_{i2}+\dots+\beta_px_{ip}+\varepsilon_i,
\quad i=1,\dots,n
$$

其中 $p$ 是不含截距的解释变量个数，参数总数为 $q=p+1$。为了写得紧凑，定义：

$$
x_i^T=(1,x_{i1},\dots,x_{ip})
$$

则单个观测可写成：

$$
Y_i=x_i^T\beta+\varepsilon_i
$$

所有观测合在一起就是矩阵形式：

$$
Y=X\beta+\varepsilon
$$

好处是：多个变量、多个参数和多个正规方程可以用一条矩阵公式处理；第 4 章的推断和第 7 章的诊断也都建立在这个形式上。

**具体解决 / Konkrete Loesung.**  
把数据表翻译成矩阵：

$$
Y=
\begin{pmatrix}
Y_1\\
\vdots\\
Y_n
\end{pmatrix},
\quad
X=
\begin{pmatrix}
1 & x_{11} & \dots & x_{1p}\\
\vdots & \vdots & & \vdots\\
1 & x_{n1} & \dots & x_{np}
\end{pmatrix},
\quad
\beta=
\begin{pmatrix}
\beta_0\\
\beta_1\\
\vdots\\
\beta_p
\end{pmatrix}
$$

第一列全为 $1$，对应截距 $\beta_0$。如果模型不含截距，矩阵和自由度公式都会改变；本章默认含截距。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)

X <- model.matrix(fit)
y <- mtcars[["mpg"]]

dim(X)
head(X)
head(y)
```

练习 / Uebung: `model.matrix(fit)` 的第一列为什么全是 $1$？  
答案 / Antwort: 因为默认模型含截距，第一列对应 $\beta_0$。

## 3.2 多元线性模型假设 / Modellannahmen

**遇到的问题 / Problem.**  
写出 $Y=X\beta+\varepsilon$ 只是模型形式。要做估计量性质、方差和后续推断，还必须说明误差项如何 behave。

Auf Deutsch: Die Matrixschreibweise allein reicht nicht; wir brauchen Annahmen ueber den Fehlervektor.

**可能的解决路径 / Moegliche Wege.**

- 不写任何假设：只能做描述性拟合，无法推出标准误和区间。
- 只假设 $E(\varepsilon)=0$：能讨论无偏性。
- 再加同方差和独立：能得到方差公式和 Gauss-Markov 结论。
- 再加正态性：能做精确的 $t$ 和 $F$ 推断，主要进入第 4 章。

**本课采用的方法及好处 / Methode und Vorteil.**  
讲义列出：

$$
E(\varepsilon_i)=0
$$

等价地：

$$
E(\varepsilon)=0
$$

同方差：

$$
V(\varepsilon_i)=\sigma^2
$$

独立误差：

$$
\varepsilon_1,\dots,\varepsilon_n\ \text{sind unabhaengig}
$$

因此：

$$
V(\varepsilon)=\sigma^2I
$$

正态误差假设：

$$
\varepsilon\sim N(0,\sigma^2I)
$$

好处是：这些假设能分别支撑 $\widehat{\beta}$ 的无偏性、方差公式和正态分布。

**具体解决 / Konkrete Loesung.**  
做题时把假设分层写清楚：

- 只证明 $E(\widehat{\beta})=\beta$：需要 $E(\varepsilon)=0$。
- 推导 $V(\widehat{\beta})=\sigma^2(X^TX)^{-1}$：还需要 $V(\varepsilon)=\sigma^2I$。
- 写 $\widehat{\beta}\sim N(\beta,\sigma^2(X^TX)^{-1})$：还需要 $\varepsilon\sim N(0,\sigma^2I)$。


- **无偏性**：只需要误差均值为0
- **方差公式**：还需要误差同方差，不相关
- **正态推断**：还需要误差正态分布，

主教材补充：如果解释变量是随机抽样得到的，很多假设应理解为“条件于设计矩阵 $X$”，例如 $E(\varepsilon\mid X)=0$ 和 $V(\varepsilon\mid X)=\sigma^2I$。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)

par(mfrow = c(1, 2))
plot(fitted(fit), residuals(fit),
     xlab = "fitted values",
     ylab = "residuals",
     main = "Residuals vs fitted")
abline(h = 0, col = "gray")
qqnorm(residuals(fit))
qqline(residuals(fit), col = "red")
par(mfrow = c(1, 1))
```

练习 / Uebung: 残差图呈扇形时，哪条假设可能有问题？  
答案 / Antwort: 同方差假设 $V(\varepsilon_i)=\sigma^2$ 可能有问题。

## 3.3 参数解释：控制其他变量 / Interpretation bei Kontrolle anderer Variablen

**遇到的问题 / Problem.**  
多元回归中 $\beta_k$ 不再只是 $x_k$ 和 $Y$ 的简单关系。它表示在其他变量固定时，$x_k$ 增加一单位对应的条件均值变化。如果忘记“其他变量固定”，就会误读参数。

Auf Deutsch: $\beta_k$ beschreibt den Effekt von $x_k$, wenn alle anderen Kovariaten konstant gehalten werden.

**可能的解决路径 / Moegliche Wege.**

- 只看散点图 $Y$ vs $x_k$：会混入其他变量影响。
- 做简单回归 $Y\sim x_k$：仍然不能控制混杂因素。
- 做多元回归：在模型中同时放入相关协变量，解释“调整后的”作用。

**本课采用的方法及好处 / Methode und Vorteil.**  
条件均值为：

$$
E(Y_i\mid x_i)=\beta_0+\beta_1x_{i1}+\dots+\beta_px_{ip}
$$

若只让第 $k$ 个解释变量增加 $1$，其他变量保持不变，则：

$$
E(Y\mid x_k+1,\text{others})-E(Y\mid x_k,\text{others})=\beta_k
$$

好处是：$\beta_k$ 可以理解为在考虑其他变量后的影响，也就是讲义中说的 Confounder-Korrektur。

**具体解决 / Konkrete Loesung.**  
解释 $\beta_k$ 时必须写四件事：

1. $Y$ 的单位是什么。
2. $x_k$ 的单位是什么。
3. 哪些其他变量被固定。
4. 这是平均变化，不是每个个体的确定变化。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit_simple <- lm(mpg ~ wt, data = mtcars)
fit_multiple <- lm(mpg ~ wt + hp + cyl, data = mtcars)

coef(fit_simple)
coef(fit_multiple)
```

练习 / Uebung: 为什么 `wt` 的系数在两个模型中可能不同？  
答案 / Antwort: 因为多元模型把 `hp` 和 `cyl` 固定住后解释 `wt`，简单模型没有控制这些变量。

## 3.4 最小二乘估计 / KQ-Schaetzer

**遇到的问题 / Problem.**  
多元模型有 $p+1$ 个参数。我们需要一个统一标准来估计整组参数，而不是逐个手动调整。

Auf Deutsch: Wir brauchen einen Schaetzer fuer den gesamten Parametervektor $\beta$.

**可能的解决路径 / Moegliche Wege.**

- 数值搜索：可行，但看不清结构。
- 逐个变量做简单回归：不是多元回归估计。
- 最小化残差平方和：得到统一、可解释、可推广的估计公式。

**本课采用的方法及好处 / Methode und Vorteil.**  
最小二乘估计定义为：

$$
\widehat{\beta}
=
\arg\min_{\beta}
(Y-X\beta)^T(Y-X\beta)
$$

残差向量为：

$$
\widehat{\varepsilon}=Y-X\widehat{\beta}
$$

如果 $X^TX$ 可逆，则：

$$
\widehat{\beta}=(X^TX)^{-1}X^TY
$$

正规方程：

$$
X^T\widehat{\varepsilon}=0
$$

好处是：多元回归和简单回归使用同一个思想；残差与设计矩阵的每一列正交。

**具体解决 / Konkrete Loesung.**  
计算路径：

1. 构造 $X$ 和 $Y$。
2. 计算 $X^TX$。
3. 检查 $X^TX$ 是否可逆。
4. 用 $(X^TX)^{-1}X^TY$ 得到 $\widehat{\beta}$。
5. 用 $Y-X\widehat{\beta}$ 得到残差。

Faraway 补充：在实际 R 中通常不显式计算逆矩阵；`lm()` 使用更稳定的 QR 分解。手算矩阵公式主要用于理解。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)

X <- model.matrix(fit)
y <- mtcars[["mpg"]]

beta_manual <- solve(t(X) %*% X, t(X) %*% y)
beta_lm <- coef(fit)

cbind(manual = as.vector(beta_manual),
      lm = beta_lm)

res <- y - X %*% beta_manual
t(X) %*% res
```

练习 / Uebung: `t(X) %*% res` 为什么应接近零向量？  
答案 / Antwort: 因为正规方程 $X^T\widehat{\varepsilon}=0$ 表示残差与 $X$ 的每一列正交。

## 3.5 产品和矩阵、秩与可识别性 / Produktsummenmatrix, Rang und Identifizierbarkeit

**遇到的问题 / Problem.**  
公式 $\widehat{\beta}=(X^TX)^{-1}X^TY$ 要求 $X^TX$ 可逆。若解释变量之间存在线性依赖，模型参数就不能唯一确定。

Auf Deutsch: Eindeutige Schaetzung verlangt vollen Spaltenrang der Designmatrix.

**可能的解决路径 / Moegliche Wege.**

- 忽略线性依赖：软件可能给出 `NA` 或不稳定系数。
- 删除冗余变量：恢复 full rank。
- 重新编码分类变量：避免 dummy trap。
- 使用约束或正则化：属于后续扩展。

**本课采用的方法及好处 / Methode und Vorteil.**  
讲义称 $X^TX$ 为 Produktsummenmatrix。它包含变量的和、平方和与交叉乘积。若 $X$ 有满列秩：

$$
\operatorname{rank}(X)=p+1=q
$$

则：

$$
X^TX\ \text{ist invertierbar}
$$

好处是：满秩条件直接说明参数是否唯一可估计。

**具体解决 / Konkrete Loesung.**  
常见不可识别情形：

- 同一个变量用两个单位同时放入模型，如公斤和磅。
- 分类变量有截距时放入所有类别 dummy。
- 一个变量是其他变量的线性组合。

若 $X^TX$ 不可逆，正规方程可能有无穷多个解，$\beta$ 至少部分不可识别。

**R 练习 / R-Uebung.**

```r
data(mtcars)
df <- mtcars
df[["wt_kg_like"]] <- df[["wt"]] * 1000

fit_bad <- lm(mpg ~ wt + wt_kg_like + hp, data = df)
summary(fit_bad)

X_bad <- model.matrix(fit_bad)
c(columns = ncol(X_bad),
  rank = qr(X_bad)[["rank"]])
```

练习 / Uebung: 如果 `rank < columns`，说明什么？  
答案 / Antwort: 说明设计矩阵列线性相关，至少有一个参数不能唯一估计。

## 3.6 $\widehat{\beta}$ 的性质 / Eigenschaften des KQ-Schaetzers

**遇到的问题 / Problem.**  
得到 $\widehat{\beta}$ 后，还要知道它是否平均等于真值、波动多大，以及在正态假设下服从什么分布。

Auf Deutsch: Nach der Schaetzung fragen wir nach Erwartung, Varianz und Verteilung.

**可能的解决路径 / Moegliche Wege.**

- 只报告系数：不够。
- 用模拟或 bootstrap：有帮助，但不是讲义主线。
- 用模型假设推导估计量性质：给出标准误和后续推断基础。

**本课采用的方法及好处 / Methode und Vorteil.**  
在 $E(\varepsilon)=0$ 下：

$$
E(\widehat{\beta})=\beta
$$

在 $V(\varepsilon)=\sigma^2I$ 下：

$$
V(\widehat{\beta})=\sigma^2(X^TX)^{-1}
$$

若进一步有 $\varepsilon\sim N(0,\sigma^2I)$，则：

$$
\widehat{\beta}\sim N\left(\beta,\sigma^2(X^TX)^{-1}\right)
$$

好处是：每个系数的标准误来自矩阵 $\widehat{\sigma}^2(X^TX)^{-1}$ 的对角线。

**具体解决 / Konkrete Loesung.**  
第 $j$ 个系数的标准误为：

$$
\widehat{\operatorname{se}}(\widehat{\beta}_j)
=
\widehat{\sigma}
\sqrt{\left[(X^TX)^{-1}\right]_{jj}}
$$

Faraway 补充：Gauss-Markov theorem 说明，在 $E(\varepsilon)=0$ 且 $V(\varepsilon)=\sigma^2I$ 时，最小二乘在线性无偏估计量中方差最小，即 BLUE。但这不意味着任何情况下都必须用 OLS；若误差相关、异方差、重尾或强共线，后续会考虑 GLS、WLS、稳健估计或正则化。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)
X <- model.matrix(fit)

sigma_hat <- sigma(fit)
cov_unscaled <- solve(t(X) %*% X)
se_manual <- sigma_hat * sqrt(diag(cov_unscaled))

cbind(manual = se_manual,
      lm = summary(fit)[["coefficients"]][, "Std. Error"])
```

练习 / Uebung: 为什么标准误和 $(X^TX)^{-1}$ 有关？  
答案 / Antwort: 因为 $V(\widehat{\beta})=\sigma^2(X^TX)^{-1}$，设计矩阵的信息量决定估计量波动。

## 3.7 Hat-Matrix 与 Residualmatrix / Hat-Matrix und Residualmatrix

**遇到的问题 / Problem.**  
我们知道 $\widehat{Y}=X\widehat{\beta}$，但还想更直接地理解：拟合值如何由原始 $Y$ 得到？残差如何由原始 $Y$ 得到？

Auf Deutsch: Fitted values und Residuen lassen sich als lineare Transformationen von $Y$ schreiben.

**可能的解决路径 / Moegliche Wege.**

- 逐个观测计算拟合值：不利于理论推导。
- 用矩阵统一表示：得到投影、杠杆值和残差性质。

**本课采用的方法及好处 / Methode und Vorteil.**  
定义 Hat-Matrix：

$$
P=X(X^TX)^{-1}X^T
$$

则：

$$
\widehat{Y}=P Y
$$

定义 Residualmatrix：

$$
Q=I-P
$$

则：

$$
\widehat{\varepsilon}=Y-\widehat{Y}=(I-P)Y=QY
$$

好处是：拟合值是 $Y$ 在 $X$ 的列空间上的投影，残差是垂直于模型空间的剩余部分。

**具体解决 / Konkrete Loesung.**  
直觉图像：

- $X$ 的列向量张成模型空间 / Modellraum。
- $\widehat{Y}$ 是 $Y$ 投影到模型空间的点。
- $\widehat{\varepsilon}$ 是 $Y-\widehat{Y}$，与模型空间正交。

主教材提醒：残差 $\widehat{\varepsilon}_i$ 不是误差 $\varepsilon_i$。误差不可观测，残差只是基于拟合模型得到的估计或预测。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)

X <- model.matrix(fit)
y <- mtcars[["mpg"]]

P <- X %*% solve(t(X) %*% X) %*% t(X)
Q <- diag(nrow(X)) - P

yhat_manual <- P %*% y
res_manual <- Q %*% y

max(abs(as.vector(yhat_manual) - fitted(fit)))
max(abs(as.vector(res_manual) - residuals(fit)))
```

练习 / Uebung: 为什么 Faraway 说通常不想显式计算 $P$？  
答案 / Antwort: 因为 $P$ 是 $n\times n$ 矩阵，大数据中会很大；实际计算常用更稳定和高效的方法。

## 3.8 投影矩阵性质 / Eigenschaften von Projektionsmatrizen

**遇到的问题 / Problem.**  
Hat matrix 和 residual matrix 不只是计算工具。它们的投影性质解释了残差正交、自由度和方差结构。

Auf Deutsch: Die algebraischen Eigenschaften von $P$ und $Q$ erklaeren die Geometrie der Regression.

**可能的解决路径 / Moegliche Wege.**

- 只记公式：容易忘。
- 从投影理解：一次投影后，再投影不会改变；拟合部分和残差部分互相正交。

**本课采用的方法及好处 / Methode und Vorteil.**  
$P$ 和 $Q$ 都是投影矩阵：

$$
P^T=P,\quad P^2=P
$$

$$
Q^T=Q,\quad Q^2=Q
$$

并且：

$$
PQ=QP=0
$$

含义：

- $P^2=P$：对拟合值再做一次同样回归，结果不变。
- $PQ=0$：拟合部分和残差部分正交。
- 若含截距，残差和常数列正交，因此 $\sum_i\widehat{\varepsilon}_i=0$。

讲义还给出：

$$
V(\widehat{Y})=\sigma^2P
$$

$$
V(\widehat{\varepsilon})=\sigma^2Q
$$

**具体解决 / Konkrete Loesung.**  
用 trace 理解自由度：

$$
\operatorname{tr}(P)=\operatorname{rank}(P)=p+1=q
$$

$$
\operatorname{tr}(Q)=n-q=n-(p+1)
$$

这解释了为什么残差平方和的自由度是 $n-(p+1)$。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)
X <- model.matrix(fit)

P <- X %*% solve(t(X) %*% X) %*% t(X)
Q <- diag(nrow(X)) - P

c(sym_P = max(abs(P - t(P))),
  idem_P = max(abs(P %*% P - P)),
  sym_Q = max(abs(Q - t(Q))),
  idem_Q = max(abs(Q %*% Q - Q)),
  orthogonal = max(abs(P %*% Q)))

c(trace_P = sum(diag(P)),
  model_rank = qr(X)[["rank"]],
  trace_Q = sum(diag(Q)),
  df_residual = df.residual(fit))
```

练习 / Uebung: 为什么 $\operatorname{tr}(Q)$ 等于残差自由度？  
答案 / Antwort: 因为 $Q$ 投影到残差空间，残差空间维度是 $n$ 减去模型空间维度 $p+1$。

## 3.9 误差方差估计 / Erwartungstreue Schaetzung von $\sigma^2$

**遇到的问题 / Problem.**  
$\sigma^2$ 是真实误差的方差，但真实误差看不到。我们只能用残差平方和估计它，并且必须扣除已经估计的参数个数。

Auf Deutsch: Die Fehlervarianz wird aus der Residuenquadratsumme mit den passenden Freiheitsgraden geschaetzt.

**可能的解决路径 / Moegliche Wege.**

- 除以 $n$：类似 ML 形式，但会低估 $\sigma^2$。
- 除以 $n-(p+1)$：补偿估计 $p+1$ 个参数损失的自由度，是无偏估计。

**本课采用的方法及好处 / Methode und Vorteil.**  
讲义公式：

$$
\widehat{\sigma}^2
=
\frac{1}{n-(p+1)}
\widehat{\varepsilon}^T\widehat{\varepsilon}
=
\frac{1}{n-(p+1)}
\sum_{i=1}^n \widehat{\varepsilon}_i^2
$$

好处是：

$$
E(\widehat{\sigma}^2)=\sigma^2
$$

推导核心：

$$
E(\widehat{\varepsilon}^T\widehat{\varepsilon})
=
\sigma^2\operatorname{tr}(Q)
=
\sigma^2(n-(p+1))
$$

**具体解决 / Konkrete Loesung.**  
做题时先确认参数个数：

- 不含截距解释变量个数：$p$。
- 含截距参数总数：$p+1$。
- 残差自由度：$n-(p+1)$。
- 误差方差估计：$RSS/\text{df}_{res}$。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)

rss <- sum(residuals(fit)^2)
df_res <- df.residual(fit)

c(sigma2_manual = rss / df_res,
  sigma2_lm = sigma(fit)^2,
  df_residual = df_res)
```

练习 / Uebung: 这里为什么不是除以 $n$？  
答案 / Antwort: 因为已经用数据估计了 $p+1$ 个参数，残差空间只剩 $n-(p+1)$ 个自由度。

## 3.10 讲义例子：Lesetest / Beispiel Lesetest

**遇到的问题 / Problem.**  
讲义 Lesetest 例子研究小学生阅读测试错误数。潜在解释变量包括性别、年级、学校阅读时间、课外阅读、Gameboy、电视等。这里的目标不是只找一个变量，而是同时看多个因素和阅读错误数的关系。

Auf Deutsch: Im Lesetest-Beispiel wird die Fehlerzahl durch mehrere Kovariaten gleichzeitig modelliert.

**可能的解决路径 / Moegliche Wege.**

- 分别比较每个变量：无法控制其他因素。
- 建立多元线性模型：同时纳入多个协变量。
- 若班级结构进一步重要：后续可考虑 mixed model；第 3 讲先看普通多元线性模型。

**本课采用的方法及好处 / Methode und Vorteil.**  
讲义中的模型形式可写为：

$$
E(Y)=
\beta_0
+\beta_1GE
+\beta_2JG
+\beta_3LZ
+\beta_4WOL
+\beta_5WOG
+\beta_6WOTV
$$

变量含义：

- $Y$: 阅读测试错误数。
- $GE$: 性别指标，例如男孩为 $1$。
- $JG$: 年级指标，例如三年级为 $1$、四年级为 $0$。
- $LZ$: 学校阅读时间。
- $WOL$: 课外阅读频率。
- $WOG$: Gameboy 频率。
- $WOTV$: 看电视频率。

好处是：$\beta_1$ 可解释为在其他变量固定时的性别差异；$\beta_4$ 可解释为阅读频率增加一个等级时，平均错误数的变化。

**具体解决 / Konkrete Loesung.**  
解释这类模型时：

1. 明确 $Y$ 是计数型错误数，但本讲暂用线性模型近似。
2. 指标变量的系数表示组间平均差异。
3. 有序等级变量若当连续变量用，隐含“等级 1 到 2”和“等级 4 到 5”的平均差异相同。
4. 若学生来自不同班级，普通线性模型的独立性可能需要后续 mixed model 检查。

**R 练习 / R-Uebung.** 用内置 `swiss` 数据模拟“多个解释变量解释一个教育/社会指标”的过程。

```r
data(swiss)
fit <- lm(Fertility ~ Education + Examination + Agriculture,
          data = swiss)

summary(fit)
coef(fit)
model.matrix(fit)[1:5, ]
```

练习 / Uebung: `Education` 的系数应如何解释？  
答案 / Antwort: 在 `Examination` 和 `Agriculture` 固定时，`Education` 增加一单位，`Fertility` 的条件均值平均改变该系数个单位。

## 3.11 R 中如何把矩阵理论落地 / Matrix Theorie in R

**遇到的问题 / Problem.**  
公式很多：$X$, $X^TX$, $\widehat{\beta}$, $P$, $Q$, $RSS$, $\widehat{\sigma}$。如果只看公式，不知道它们在 R 的哪里，会很抽象。

Auf Deutsch: Die Matrixformeln sollen mit konkreten R-Objekten verbunden werden.

**可能的解决路径 / Moegliche Wege.**

- 只看 `summary(fit)`：快，但不理解来源。
- 手算所有矩阵：理解深，但实际中不稳定。
- 两者结合：用 `lm()` 拟合，再用矩阵对象验证关键公式。

**本课采用的方法及好处 / Methode und Vorteil.**  
Faraway 的辅助教材强调：`lm()` 对象中已经存了很多回归量，可以抽取出来理解矩阵公式。

常用 R 函数：

- `model.matrix(fit)`: 设计矩阵 $X$。
- `coef(fit)`: $\widehat{\beta}$。
- `fitted(fit)`: $\widehat{Y}$。
- `residuals(fit)`: $\widehat{\varepsilon}$。
- `deviance(fit)`: $RSS=\widehat{\varepsilon}^T\widehat{\varepsilon}$。
- `df.residual(fit)`: $n-(p+1)$。
- `sigma(fit)`: $\widehat{\sigma}$。
- `hatvalues(fit)`: $P$ 的对角线 $h_{ii}$。

**具体解决 / Konkrete Loesung.**  
把理论和 R 一一对应：

$$
\widehat{\beta}\leftrightarrow \texttt{coef(fit)}
$$

$$
\widehat{Y}\leftrightarrow \texttt{fitted(fit)}
$$

$$
\widehat{\varepsilon}\leftrightarrow \texttt{residuals(fit)}
$$

$$
\widehat{\varepsilon}^T\widehat{\varepsilon}\leftrightarrow \texttt{deviance(fit)}
$$

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)

c(rss_from_deviance = deviance(fit),
  rss_manual = sum(residuals(fit)^2))

c(df_from_lm = df.residual(fit),
  n_minus_rank = nrow(model.matrix(fit)) - qr(model.matrix(fit))[["rank"]])

head(hatvalues(fit))
summary(fit)[["coefficients"]]
```

练习 / Uebung: `hatvalues(fit)` 是哪个矩阵的对角线？  
答案 / Antwort: 它是 Hat-Matrix $P=X(X^TX)^{-1}X^T$ 的对角线。

## 3.12 选择题 / Multiple-Choice-Fragen

1. 多元线性回归的矩阵形式是什么？  
   A. $Y=X\beta+\varepsilon$  
   B. $Y=\beta+\varepsilon X$  
   C. $X=Y\beta+\varepsilon$  
   D. $Y=X+\beta$  
   答案 / Antwort: A

2. 在含截距的多元线性回归中，设计矩阵 $X$ 的第一列通常是什么？  
   A. 全是 $1$  
   B. 全是 $0$  
   C. 响应变量 $Y$  
   D. 残差  
   答案 / Antwort: A

3. $\widehat{\beta}=(X^TX)^{-1}X^TY$ 需要什么条件？  
   A. $X$ 满列秩，$X^TX$ 可逆  
   B. $Y$ 必须全为正  
   C. 所有残差都等于 $1$  
   D. $p>n$  
   答案 / Antwort: A

4. 正规方程 $X^T\widehat{\varepsilon}=0$ 的几何含义是什么？  
   A. 残差向量与设计矩阵列空间正交  
   B. 残差向量等于拟合值  
   C. 所有系数都为 $0$  
   D. $X$ 没有截距  
   答案 / Antwort: A

5. Hat-Matrix $P$ 满足哪条性质？  
   A. $P^T=P$ 且 $P^2=P$  
   B. $P^T=-P$  
   C. $P^2=0$  
   D. $P=Q$  
   答案 / Antwort: A

6. Residualmatrix $Q$ 如何定义？  
   A. $Q=I-P$  
   B. $Q=X^TX$  
   C. $Q=P^{-1}$  
   D. $Q=XY$  
   答案 / Antwort: A

7. 为什么 $\widehat{\sigma}^2$ 的分母是 $n-(p+1)$？  
   A. 因为含截距模型估计了 $p+1$ 个参数  
   B. 因为 $Y$ 有 $p+1$ 个观测  
   C. 因为 $P$ 总是零矩阵  
   D. 因为不需要自由度  
   答案 / Antwort: A

8. 多元回归中 $\beta_k$ 的正确解释是什么？  
   A. 其他变量固定时，$x_k$ 增加一单位，$Y$ 条件均值的平均变化  
   B. $x_k$ 与 $Y$ 的简单相关系数  
   C. $Y$ 的总方差  
   D. 第 $k$ 个残差  
   答案 / Antwort: A

9. 如果 $X$ 的列线性相关，会发生什么？  
   A. 参数可能不可唯一识别  
   B. 残差一定全为零  
   C. $R^2$ 一定为零  
   D. 正态性自动成立  
   答案 / Antwort: A

10. $P$ 的对角线在 R 中常用哪个函数查看？  
    A. `hatvalues(fit)`  
    B. `hist(fit)`  
    C. `mean(fit)`  
    D. `names(data)`  
    答案 / Antwort: A

## 3.13 本章总结 / Zusammenfassung

中文总结：

- 多元线性模型把简单线性回归扩展到多个解释变量：$Y_i=\beta_0+\beta_1x_{i1}+\dots+\beta_px_{ip}+\varepsilon_i$。
- 矩阵形式 $Y=X\beta+\varepsilon$ 是本章核心，$X$ 是设计矩阵，第一列通常对应截距。
- 多元模型中的 $\beta_k$ 表示在其他解释变量固定时，$x_k$ 增加一单位带来的条件均值平均变化。
- 最小二乘估计最小化 $(Y-X\beta)^T(Y-X\beta)$；若 $X^TX$ 可逆，则 $\widehat{\beta}=(X^TX)^{-1}X^TY$。
- 正规方程 $X^T\widehat{\varepsilon}=0$ 表示残差与设计矩阵每一列正交。
- Hat-Matrix $P=X(X^TX)^{-1}X^T$ 给出 $\widehat{Y}=PY$；Residualmatrix $Q=I-P$ 给出 $\widehat{\varepsilon}=QY$。
- $P$ 和 $Q$ 是投影矩阵，满足对称、幂等和互相正交。
- $\widehat{\sigma}^2=\frac{1}{n-(p+1)}\widehat{\varepsilon}^T\widehat{\varepsilon}$ 是 $\sigma^2$ 的无偏估计。
- full rank 是参数唯一估计的关键；若解释变量线性相关，会出现不可识别或软件删掉系数。

Deutsche Zusammenfassung:

- Das multiple lineare Regressionsmodell lautet $Y=X\beta+\varepsilon$.
- Die Designmatrix $X$ enthaelt die Einflussgroessen und meist eine Eins-Spalte fuer den Achsenabschnitt.
- $\beta_k$ beschreibt den Effekt von $x_k$, wenn alle anderen Kovariaten festgehalten werden.
- Die KQ-Schaetzung minimiert $(Y-X\beta)^T(Y-X\beta)$.
- Bei vollem Spaltenrang gilt $\widehat{\beta}=(X^TX)^{-1}X^TY$.
- Die Normalgleichungen $X^T\widehat{\varepsilon}=0$ beschreiben Orthogonalitaet der Residuen zum Modellraum.
- Die Hat-Matrix $P$ projiziert $Y$ auf den Modellraum; $Q=I-P$ projiziert auf den Residuenraum.
- Die Residuenquadratsumme wird durch $n-(p+1)$ geteilt, um $\sigma^2$ erwartungstreu zu schaetzen.

## 3.14 专业德语单词汇总 / Fachwortschatz

| Deutsch | 中文 | 公式 / Hinweis |
|---|---|---|
| multiples lineares Regressionsmodell | 多元线性回归模型 | $Y=X\beta+\varepsilon$ |
| Zielgroesse | 目标变量 | $Y$ |
| Einflussgroesse | 解释变量 | $x_j$ |
| Designmatrix | 设计矩阵 | $X$ |
| Regressionsparameter | 回归参数 | $\beta_0,\dots,\beta_p$ |
| Stoergroesse / Fehlerterm | 扰动项 / 误差项 | $\varepsilon_i$ |
| Kovariate | 协变量 | $x_j$ |
| Confounder-Korrektur | 混杂因素调整 | 固定其他 $x$ 后解释 $\beta_k$ |
| KQ-Schaetzer | 最小二乘估计量 | $\widehat{\beta}$ |
| Methode der kleinsten Quadrate | 最小二乘法 | $\min (Y-X\beta)^T(Y-X\beta)$ |
| Produktsummenmatrix | 乘积和矩阵 | $X^TX$ |
| Normalgleichungen | 正规方程 | $X^T\widehat{\varepsilon}=0$ |
| voller Spaltenrang | 满列秩 | $\operatorname{rank}(X)=p+1$ |
| Identifizierbarkeit | 可识别性 | 参数能唯一估计 |
| Hat-Matrix | 帽子矩阵 | $P=X(X^TX)^{-1}X^T$ |
| Residualmatrix | 残差矩阵 | $Q=I-P$ |
| Projektionsmatrix | 投影矩阵 | $P^2=P$, $Q^2=Q$ |
| orthogonal | 正交 | $PQ=0$ |
| Residuen | 残差 | $\widehat{\varepsilon}=Y-\widehat{Y}$ |
| angepasste Werte | 拟合值 | $\widehat{Y}=PY$ |
| Varianz-Kovarianz-Matrix | 方差-协方差矩阵 | $V(\widehat{\beta})=\sigma^2(X^TX)^{-1}$ |
| erwartungstreu | 无偏 | $E(\widehat{\beta})=\beta$ |
| Freiheitsgrade | 自由度 | $n-(p+1)$ |
| Spur | 迹 | $\operatorname{tr}(P)=p+1$ |
| Gauss-Markov-Theorem | Gauss-Markov 定理 | OLS 是 BLUE |
