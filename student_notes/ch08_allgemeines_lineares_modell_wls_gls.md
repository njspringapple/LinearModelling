# 第 8 章学生笔记：一般线性模型、WLS 与 GLS / Kapitel 8: Allgemeines lineares Modell

材料 / Grundlage:

- 讲义 / Folien: `chapter_pdfs/LinearModellingPPT/08_Das_allgemeine_lineare_Modell_Gewichtete_KQ-Methode_Autokorrelierte_und_heteroskedastische_Stor.pdf`
- 主教材 / Hauptlehrbuch: `Regression.pdf`, Chapter 4.1
- 辅助教材 / Hilfslehrbuch: Faraway, `Linear Models with R`, Chapter 6.1, Chapter 6.2

## 覆盖审核 / Abdeckungspruefung

| 来源 / Quelle | 内容点 / Inhaltspunkt | 笔记位置 / Wo |
|---|---|---|
| 讲义 8 | 一般线性模型 | 8.1 |
| 讲义 8 | 模型变换 | 8.2 |
| 讲义 8 | weighted least squares | 8.3 |
| 讲义 8 | generalized least squares | 8.4 |
| 讲义 8 | 一般 Gauss-Markov 定理 | 8.5 |
| 讲义 8 | 方差结构、AR(1) | 8.6 |
| 讲义 8 | 估计策略、ML/REML | 8.7 |
| 讲义 8 | 关于 $\beta$ 方差的推断 | 8.8 |
| 主教材 | WLS、GLS、异方差和自相关误差 | 全章 |
| Faraway | R 中 GLS/WLS 实操思想 | R 练习 |

## 8.1 为什么需要一般线性模型

**遇到的问题 / Problem.**  
前面经典线性模型假设 $V(\varepsilon)=\sigma^2I$。实际中误差可能异方差或相关，例如不同观测精度不同、时间序列误差相邻相关。

**可能的解决路径 / Moegliche Wege.**

- 继续用 OLS：系数可能仍可用，但标准误和效率有问题。
- 使用 WLS：处理独立但方差不同。
- 使用 GLS：处理一般协方差结构。

**本课采用的方法及好处 / Methode und Vorteil.**  
一般线性模型：

$$
Y=X\beta+\varepsilon,\quad E(\varepsilon)=0,\quad V(\varepsilon)=\sigma^2\Sigma
$$

其中 $\Sigma$ 不一定是单位矩阵。好处是：把异方差和相关误差纳入同一矩阵框架。

**R 练习 / R-Uebung.**

```r
data(cars)
fit_ols <- lm(dist ~ speed, data = cars)
plot(fitted(fit_ols), residuals(fit_ols))
abline(h = 0, col = "gray")
```

## 8.2 模型变换 / Transformation des Modells

**遇到的问题 / Problem.**  
若 $V(\varepsilon)=\sigma^2\Sigma$，OLS 的球形误差假设不成立。我们希望把模型变换成误差协方差为 $\sigma^2I$ 的形式。

**本课采用的方法及好处 / Methode und Vorteil.**  
若存在矩阵 $A$ 使：

$$
A\Sigma A^T=I
$$

则变换：

$$
AY=AX\beta+A\varepsilon
$$

变换后误差满足：

$$
V(A\varepsilon)=\sigma^2I
$$

好处是：可以在变换后的模型上使用普通最小二乘思想。

**具体解决 / Konkrete Loesung.**  
理论上可用 $\Sigma^{-1/2}$ 变换；实际中常直接使用 GLS 公式或软件函数。

## 8.3 加权最小二乘 / Gewichtete KQ-Methode

**遇到的问题 / Problem.**  
若误差独立但方差不同：

$$
V(\varepsilon_i)=\frac{\sigma^2}{w_i}
$$

高方差观测不应和低方差观测同等权重。

**本课采用的方法及好处 / Methode und Vorteil.**  
WLS 最小化：

$$
\sum_{i=1}^n w_i(Y_i-x_i^T\beta)^2
$$

矩阵形式：

$$
\widehat{\beta}_{WLS}
=(X^TWX)^{-1}X^TWY
$$

好处是：精度高的观测权重大，精度低的观测权重小。

**具体解决 / Konkrete Loesung.**  
权重 $w_i$ 应与误差方差成反比。若只知道方差结构的相对形式，可用估计的权重。

**R 练习 / R-Uebung.**

```r
data(cars)
fit_ols <- lm(dist ~ speed, data = cars)

weights <- 1 / cars[["speed"]]^2
fit_wls <- lm(dist ~ speed, data = cars, weights = weights)

coef(fit_ols)
coef(fit_wls)
```

## 8.4 广义最小二乘 / Verallgemeinerte KQ-Methode

**遇到的问题 / Problem.**  
WLS 只处理对角协方差矩阵。若误差之间相关，需要更一般的 GLS。

**本课采用的方法及好处 / Methode und Vorteil.**  
若：

$$
V(\varepsilon)=\sigma^2\Sigma
$$

且 $\Sigma$ 已知，则 GLS 估计为：

$$
\widehat{\beta}_{GLS}
=
(X^T\Sigma^{-1}X)^{-1}X^T\Sigma^{-1}Y
$$

好处是：当 $\Sigma$ 正确时，GLS 比 OLS 更有效。

**具体解决 / Konkrete Loesung.**  
GLS 的核心是用 $\Sigma^{-1}$ 重新定义距离：不是最小化普通残差平方和，而是最小化：

$$
(Y-X\beta)^T\Sigma^{-1}(Y-X\beta)
$$

**R 练习 / R-Uebung.** Base R 没有通用 `gls()`，下面手动构造一个简单 AR(1) 协方差示例。

```r
data(Nile)
y <- as.numeric(Nile)
t <- seq_along(y)
X <- model.matrix(~ t)
n <- length(y)
rho <- 0.6
Sigma <- outer(seq_len(n), seq_len(n), function(i, j) rho^abs(i - j))

beta_gls <- solve(t(X) %*% solve(Sigma, X),
                  t(X) %*% solve(Sigma, y))
coef(lm(y ~ t))
as.vector(beta_gls)
```

## 8.5 一般 Gauss-Markov 定理 / Allgemeines Gauss-Markov-Theorem

**遇到的问题 / Problem.**  
经典 Gauss-Markov 假设 $V(\varepsilon)=\sigma^2I$。若协方差不是单位矩阵，OLS 不再是最佳线性无偏估计。

**本课采用的方法及好处 / Methode und Vorteil.**  
在 $V(\varepsilon)=\sigma^2\Sigma$ 且 $\Sigma$ 已知时，GLS 是 BLUE。

这意味着：正确利用误差协方差结构，可以在保持无偏的同时降低估计方差。

**具体解决 / Konkrete Loesung.**  
先判断问题来自异方差还是相关性；再选择 WLS 或 GLS；最后用相应协方差矩阵计算标准误。

## 8.6 方差结构与 AR(1)

**遇到的问题 / Problem.**  
误差结构必须被参数化，否则 $\Sigma$ 太大，无法稳定估计。

**常见结构 / Hauefige Strukturen.**

- 独立异方差：

$$
\Sigma=\operatorname{diag}(\sigma_1^2,\dots,\sigma_n^2)
$$

- AR(1) 时间序列结构：

$$
\operatorname{Corr}(\varepsilon_t,\varepsilon_s)=\rho^{|t-s|}
$$

好处是：用少数参数描述误差结构。

**R 练习 / R-Uebung.**

```r
rho <- 0.7
n <- 8
Sigma_ar1 <- outer(seq_len(n), seq_len(n), function(i, j) rho^abs(i - j))
round(Sigma_ar1, 2)
```

## 8.7 估计策略、ML 与 REML

**遇到的问题 / Problem.**  
现实中 $\Sigma$ 或方差参数 $\theta$ 通常未知。要使用 GLS，必须先估计方差结构。

**可能的解决路径 / Moegliche Wege.**

- 先用 OLS 残差估计方差结构，再做 feasible GLS。
- 用 ML 同时估计 $\beta$ 和方差参数。
- 用 REML 估计方差参数，减少固定效应估计带来的偏差。

**本课采用的方法及好处 / Methode und Vorteil.**  
若 $\theta$ 已知，$\beta$ 的条件 MLE 就是加权/广义最小二乘。若 $\theta$ 未知，则需要先估计 $\theta$。REML 会在 mixed model 中再次出现。

**具体解决 / Konkrete Loesung.**  
记住三层：

1. $\Sigma$ 已知：直接 GLS。
2. $\Sigma(\theta)$ 已知形式但 $\theta$ 未知：估计 $\theta$ 后 GLS。
3. 方差结构复杂：使用 ML/REML 软件拟合。

## 8.8 关于 $\beta$ 的推断

**遇到的问题 / Problem.**  
使用 WLS/GLS 后，系数标准误也必须按照新的协方差结构计算。

**本课采用的方法及好处 / Methode und Vorteil.**  
若 $\Sigma$ 已知：

$$
V(\widehat{\beta}_{GLS})
=
\sigma^2(X^T\Sigma^{-1}X)^{-1}
$$

这替代了 OLS 的：

$$
\sigma^2(X^TX)^{-1}
$$

**具体解决 / Konkrete Loesung.**  
不要用 OLS 的标准误解释 GLS/WLS 的系数；估计方法变了，方差公式也必须变。

## 8.9 选择题 / Multiple-Choice-Fragen

1. 一般线性模型允许什么？  
   A. $V(\varepsilon)=\sigma^2\Sigma$  
   B. 只能有 $V(\varepsilon)=0$  
   C. 不能有解释变量  
   D. 只能做 logistic 回归  
   答案 / Antwort: A

2. WLS 权重应与什么相关？  
   A. 误差方差的倒数  
   B. 观测编号  
   C. 响应变量名字  
   D. 截距个数  
   答案 / Antwort: A

3. GLS 估计公式包含哪个矩阵？  
   A. $\Sigma^{-1}$  
   B. $Y^{-1}$  
   C. $\varepsilon^{-1}$  
   D. $R^2$  
   答案 / Antwort: A

4. AR(1) 结构常用于什么？  
   A. 时间序列相关误差  
   B. 只有分类变量的编码  
   C. 删除截距  
   D. 单因素 ANOVA  
   答案 / Antwort: A

5. REML 主要用于什么？  
   A. 方差参数估计  
   B. 把概率限制到 $[0,1]$  
   C. 画散点图  
   D. 计算 dummy 变量  
   答案 / Antwort: A

## 8.10 本章总结 / Zusammenfassung

中文总结：

- 一般线性模型把误差协方差从 $\sigma^2I$ 扩展为 $\sigma^2\Sigma$。
- WLS 适合独立但方差不同的观测，权重应与方差倒数相关。
- GLS 适合一般协方差结构，最小化加权二次型。
- 若协方差结构正确，GLS 是一般 Gauss-Markov 意义下的 BLUE。
- AR(1) 是常见时间序列误差相关结构。
- ML 和 REML 用于估计未知方差参数。

Deutsche Zusammenfassung:

- Das allgemeine lineare Modell erlaubt $V(\varepsilon)=\sigma^2\Sigma$.
- WLS behandelt heteroskedastische, aber unkorrelierte Fehler.
- GLS behandelt allgemeine Kovarianzstrukturen.
- Bei bekannter Kovarianzstruktur ist GLS der BLUE.
- ML und REML werden zur Schaetzung von Varianzparametern verwendet.

## 8.11 专业德语单词汇总 / Fachwortschatz

| Deutsch | 中文 | 公式 / Hinweis |
|---|---|---|
| allgemeines lineares Modell | 一般线性模型 | $V(\varepsilon)=\sigma^2\Sigma$ |
| gewichtete KQ-Methode | 加权最小二乘 | WLS |
| verallgemeinerte KQ-Methode | 广义最小二乘 | GLS |
| Varianzstruktur | 方差结构 | $\Sigma$ |
| heteroskedastische Stoerterme | 异方差误差 | 方差不同 |
| autokorrelierte Stoerterme | 自相关误差 | 误差相关 |
| AR(1)-Struktur | 一阶自回归结构 | $\rho^{|t-s|}$ |
| ML-Schaetzung | 最大似然估计 | likelihood |
| REML-Schaetzung | 限制最大似然 | variance components |
| allgemeines Gauss-Markov-Theorem | 一般 Gauss-Markov 定理 | GLS 是 BLUE |
