# 第 10 章学生笔记：Logistic 回归模型 / Kapitel 10: Das logistische Regressionsmodell

材料 / Grundlage:

- 讲义 / Folien: `chapter_pdfs/LinearModellingPPT/10_Das_logistische_Regressionsmodell.pdf`
- 主教材 / Hauptlehrbuch: `Regression.pdf`, Chapter 2.3, Chapter 5.1, Appendix B.4
- 辅助教材 / Hilfslehrbuch: Faraway, `Extending The Linear Model With R`, Chapter 2, Chapter 4, Chapter 8
- `Linear Models with R` 仅作为普通线性模型背景，不是本章主读材料

## 覆盖审核 / Abdeckungspruefung

| 来源 / Quelle | 内容点 / Inhaltspunkt | 笔记位置 / Wo |
|---|---|---|
| 讲义 10 | 二元响应问题、Bauteile 例子 | 10.1 |
| 讲义 10 | logistic 回归模型 | 10.2 |
| 讲义 10 | odds、log-odds、参数解释 | 10.3 |
| 讲义 10 | logistic regression as classification | 10.4 |
| 讲义 10 | ML 估计、MLE 性质 | 10.5 |
| 讲义 10 | MLE 存在唯一性 | 10.6 |
| 讲义 10 | Wald test、LRT、置信区间 | 10.7 |
| 讲义 10 | 模型拟合优度 | 10.8 |
| 讲义 10 | ROC、sensitivity、specificity、AUC | 10.9 |
| 讲义 10 | case-control study、propensity score | 10.10 |
| 主教材 | logit model、binary regression、likelihood inference | 全章 |
| Extending Faraway | binary response, logistic variants, GLM framework | 全章 |

## 10.1 为什么不能直接用普通线性回归

**遇到的问题 / Problem.**  
本章响应变量是二元的，例如合格/不合格、患病/未患病、成功/失败。此时 $Y\in\{0,1\}$，条件均值就是概率：

$$
E(Y\mid x)=P(Y=1\mid x)
$$

普通线性模型可能预测出小于 $0$ 或大于 $1$ 的值，不适合作为概率模型。

**可能的解决路径 / Moegliche Wege.**

- 线性概率模型：简单，但概率范围和异方差有问题。
- probit 模型：用正态分布函数连接概率。
- logistic 模型：用 logit 链接，解释为 odds ratio，课程主线采用。

**本课采用的方法及好处 / Methode und Vorteil.**  
用 logistic 函数把线性预测子映射到 $[0,1]$：

$$
\pi_i=P(Y_i=1\mid x_i)
$$

$$
\pi_i=\frac{\exp(\eta_i)}{1+\exp(\eta_i)}
$$

其中：

$$
\eta_i=\beta_0+\beta_1x_{i1}+\dots+\beta_px_{ip}
$$

好处是：概率合法，参数可以通过 odds 和 log-odds 解释。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- glm(vs ~ mpg + wt, data = mtcars, family = binomial)
summary(fit)
predict(fit, type = "response")[1:6]
```

## 10.2 Logistic 模型与 logit 链接

**遇到的问题 / Problem.**  
概率 $\pi$ 在 $[0,1]$ 内，而线性预测子 $\eta$ 在整个实数轴上。需要一个链接函数连接二者。

**本课采用的方法及好处 / Methode und Vorteil.**  
odds:

$$
\frac{\pi}{1-\pi}
$$

logit:

$$
\log\left(\frac{\pi}{1-\pi}\right)=\eta
$$

logistic 回归：

$$
\log\left(\frac{\pi_i}{1-\pi_i}\right)
=
\beta_0+\beta_1x_{i1}+\dots+\beta_px_{ip}
$$

好处是：log-odds 对解释变量线性，而概率自动限制在 $[0,1]$。

**具体解决 / Konkrete Loesung.**  
给定 $\eta$，概率为：

$$
\pi=\frac{\exp(\eta)}{1+\exp(\eta)}
$$

给定概率 $\pi$，logit 为：

$$
\operatorname{logit}(\pi)=\log\left(\frac{\pi}{1-\pi}\right)
$$

**R 练习 / R-Uebung.**

```r
eta <- seq(-5, 5, length.out = 100)
pi <- exp(eta) / (1 + exp(eta))
plot(eta, pi, type = "l",
     xlab = "eta",
     ylab = "probability")
```

## 10.3 参数解释：odds ratio

**遇到的问题 / Problem.**  
logistic 回归系数不是概率差。若直接说“$x$ 增加一单位概率增加 $\beta$”，就是错的。

**本课采用的方法及好处 / Methode und Vorteil.**  
若 $x_j$ 增加一单位，其他变量固定，则 log-odds 改变 $\beta_j$：

$$
\Delta \log\left(\frac{\pi}{1-\pi}\right)=\beta_j
$$

odds 乘以：

$$
\exp(\beta_j)
$$

这就是 odds ratio。好处是：参数解释清楚且不依赖当前概率尺度。

**具体解决 / Konkrete Loesung.**  
报告 logistic 系数时，通常同时报告 $\widehat{\beta}_j$ 和 $\exp(\widehat{\beta}_j)$。若 $\exp(\widehat{\beta}_j)>1$，表示 odds 增大；若小于 $1$，表示 odds 减小。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- glm(vs ~ mpg + wt, data = mtcars, family = binomial)

coef(fit)
exp(coef(fit))
confint.default(fit)
exp(confint.default(fit))
```

## 10.4 作为分类问题 / Klassifikation

**遇到的问题 / Problem.**  
logistic 回归输出概率，但实践中常需要分类，例如预测 $Y=1$ 或 $Y=0$。

**本课采用的方法及好处 / Methode und Vorteil.**  
给定阈值 $c$：

$$
\widehat{Y}=
\begin{cases}
1,& \widehat{\pi}\ge c\\
0,& \widehat{\pi}<c
\end{cases}
$$

常用 $c=0.5$，但阈值应根据错误成本调整。

**具体解决 / Konkrete Loesung.**  
分类性能不仅看准确率，还要看 sensitivity 和 specificity。若阳性很少，准确率可能误导。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- glm(vs ~ mpg + wt, data = mtcars, family = binomial)
p <- predict(fit, type = "response")
pred <- ifelse(p >= 0.5, 1, 0)
table(observed = mtcars[["vs"]], predicted = pred)
mean(pred == mtcars[["vs"]])
```

## 10.5 最大似然估计 / ML-Schaetzung

**遇到的问题 / Problem.**  
二元响应不是正态误差模型，不能用最小二乘作为主线。需要基于 Bernoulli 分布建立似然。

**本课采用的方法及好处 / Methode und Vorteil.**  
若：

$$
Y_i\sim Bernoulli(\pi_i)
$$

似然为：

$$
L(\beta)=\prod_{i=1}^n \pi_i^{Y_i}(1-\pi_i)^{1-Y_i}
$$

log-likelihood:

$$
\ell(\beta)=
\sum_{i=1}^n
\left[
Y_i\log(\pi_i)+(1-Y_i)\log(1-\pi_i)
\right]
$$

好处是：MLE 适用于非正态响应，是 GLM 的核心估计思想。

**具体解决 / Konkrete Loesung.**  
logistic 回归没有像 OLS 那样的闭式解，通常用数值迭代求 MLE。R 的 `glm(..., family = binomial)` 自动完成。

## 10.6 MLE 的存在唯一性

**遇到的问题 / Problem.**  
logistic 回归中，如果某个变量能完美分开 $Y=0$ 和 $Y=1$，MLE 可能不存在或趋向无穷大。

**可能的解决路径 / Moegliche Wege.**

- 检查 separation。
- 合并类别或删除导致完美分离的变量。
- 使用惩罚估计或 Firth 修正。

**本课采用的方法及好处 / Methode und Vorteil.**  
讲义提醒：logistic MLE 的存在唯一性不是自动成立。和线性模型满秩不同，二元响应还要考虑分离问题。

**R 练习 / R-Uebung.**

```r
x <- c(1, 2, 3, 4, 5, 6)
y <- c(0, 0, 0, 1, 1, 1)
fit_sep <- glm(y ~ x, family = binomial)
summary(fit_sep)
```

练习 / Uebung: 如果看到系数和标准误特别大，应想到什么？  
答案 / Antwort: 可能存在 separation 或近似 separation。

## 10.7 Wald test、LRT 与置信区间

**遇到的问题 / Problem.**  
logistic 回归没有正态线性模型中的 $F$ 检验主线。推断通常基于 MLE 的渐近性质。

**本课采用的方法及好处 / Methode und Vorteil.**  
MLE 渐近正态：

$$
\widehat{\beta}\approx N(\beta,\widehat{V}(\widehat{\beta}))
$$

Wald 统计量：

$$
z_j=\frac{\widehat{\beta}_j}{\widehat{se}(\widehat{\beta}_j)}
$$

LRT 比较完整模型和受限模型：

$$
LR=2(\ell_{full}-\ell_{restricted})
\approx \chi^2_r
$$

好处是：LRT 和 Wald test 都能检验线性假设；LRT 常更稳健，尤其小样本时。

**R 练习 / R-Uebung.**

```r
data(mtcars)
full <- glm(vs ~ mpg + wt, data = mtcars, family = binomial)
reduced <- glm(vs ~ mpg, data = mtcars, family = binomial)

summary(full)
anova(reduced, full, test = "Chisq")
```

## 10.8 模型拟合优度 / Modellguete

**遇到的问题 / Problem.**  
logistic 回归不能直接使用普通线性模型的残差平方和和 $R^2$ 解释。需要用 deviance、likelihood 和分类指标。

**本课采用的方法及好处 / Methode und Vorteil.**  
Deviance 衡量模型相对 saturated model 的差距。Null deviance 与 residual deviance 的比较类似“模型是否比只有截距好”。

伪 $R^2$ 有多种定义，不应和线性模型 $R^2$ 完全等同解释。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit <- glm(vs ~ mpg + wt, data = mtcars, family = binomial)

c(null_deviance = fit[["null.deviance"]],
  residual_deviance = fit[["deviance"]],
  aic = AIC(fit))
```

## 10.9 ROC、sensitivity、specificity 与 AUC

**遇到的问题 / Problem.**  
分类阈值改变会改变 true positive rate 和 false positive rate。单一阈值下的准确率不能完整评价模型。

**本课采用的方法及好处 / Methode und Vorteil.**  
sensitivity:

$$
\frac{TP}{TP+FN}
$$

specificity:

$$
\frac{TN}{TN+FP}
$$

ROC 曲线画的是不同阈值下：

$$
TPR=\text{sensitivity}
$$

对：

$$
FPR=1-\text{specificity}
$$

AUC 衡量随机抽一个阳性和一个阴性时，模型给阳性更高分数的概率。

**R 练习 / R-Uebung.** 不用额外包，手动计算 ROC 点。

```r
data(mtcars)
fit <- glm(vs ~ mpg + wt, data = mtcars, family = binomial)
p <- predict(fit, type = "response")
y <- mtcars[["vs"]]

thresholds <- seq(0, 1, by = 0.05)
roc <- sapply(thresholds, function(cut) {
  pred <- ifelse(p >= cut, 1, 0)
  TP <- sum(pred == 1 & y == 1)
  TN <- sum(pred == 0 & y == 0)
  FP <- sum(pred == 1 & y == 0)
  FN <- sum(pred == 0 & y == 1)
  c(sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP))
})
roc <- t(roc)
plot(1 - roc[, "specificity"], roc[, "sensitivity"],
     type = "b", xlab = "1 - specificity", ylab = "sensitivity")
abline(0, 1, lty = 2, col = "gray")
```

## 10.10 Case-control 与 propensity score

**遇到的问题 / Problem.**  
病例-对照研究中，样本中病例比例由抽样设计决定，不等于总体患病率。propensity score 则用 logistic 回归估计接受处理的概率。

**本课采用的方法及好处 / Methode und Vorteil.**

- Case-control: logistic 回归仍可用于估计 odds ratio，但截距和总体概率解释要谨慎。
- Propensity score:

$$
e(x)=P(T=1\mid X=x)
$$

可用于匹配、分层或加权，帮助平衡协变量。

**R 练习 / R-Uebung.**

```r
data(mtcars)
ps_fit <- glm(am ~ mpg + wt + hp, data = mtcars, family = binomial)
ps <- predict(ps_fit, type = "response")
summary(ps)
```

## 10.11 选择题 / Multiple-Choice-Fragen

1. logistic 回归中 $E(Y\mid x)$ 表示什么？  
   A. $P(Y=1\mid x)$  
   B. $x$ 的方差  
   C. 残差平方和  
   D. 样本量  
   答案 / Antwort: A

2. logit 链接是什么？  
   A. $\log(\pi/(1-\pi))$  
   B. $\pi^2$  
   C. $1-\pi$  
   D. $\log(x)$  
   答案 / Antwort: A

3. $\exp(\beta_j)$ 表示什么？  
   A. odds ratio  
   B. 残差  
   C. AIC  
   D. sensitivity  
   答案 / Antwort: A

4. LRT 比较什么？  
   A. 受限和完整模型的最大 log-likelihood  
   B. 两个变量名  
   C. 两个直方图  
   D. 两个样本量  
   答案 / Antwort: A

5. ROC 曲线横轴常是什么？  
   A. $1-\text{specificity}$  
   B. $\beta_0$  
   C. 残差  
   D. AIC  
   答案 / Antwort: A

## 10.12 本章总结 / Zusammenfassung

中文总结：

- logistic 回归用于二元响应变量，条件均值是概率 $P(Y=1\mid x)$。
- logit 链接把概率转成实数轴上的线性预测子。
- 系数解释为 log-odds 变化，$\exp(\beta)$ 是 odds ratio。
- 参数通过最大似然估计，通常没有闭式解。
- Wald test 和 likelihood-ratio test 是常用推断方法。
- 分类性能要结合阈值、sensitivity、specificity、ROC 和 AUC。
- case-control 研究中 odds ratio 可解释，但截距和总体概率要谨慎。
- propensity score 是接受处理的条件概率，可用于协变量平衡。

Deutsche Zusammenfassung:

- Logistische Regression modelliert binaere Zielgroessen.
- Der Logit-Link lautet $\log(\pi/(1-\pi))=\eta$.
- $\exp(\beta_j)$ ist ein Odds Ratio.
- Die Schaetzung erfolgt mit Maximum Likelihood.
- Wald-Tests und Likelihood-Quotienten-Tests sind zentrale Inferenzwerkzeuge.
- ROC-Kurve und AUC bewerten Klassifikationsleistung ueber Schwellenwerte hinweg.

## 10.13 专业德语单词汇总 / Fachwortschatz

| Deutsch | 中文 | 公式 / Hinweis |
|---|---|---|
| binaere Zielgroesse | 二元目标变量 | $Y\in\{0,1\}$ |
| logistisches Regressionsmodell | logistic 回归模型 | logit link |
| Wahrscheinlichkeit | 概率 | $\pi$ |
| Odds | 赔率 | $\pi/(1-\pi)$ |
| Log-Odds | 对数赔率 | $\log(\pi/(1-\pi))$ |
| Odds Ratio | 赔率比 | $\exp(\beta)$ |
| Maximum-Likelihood-Schaetzung | 最大似然估计 | MLE |
| Separation | 分离 | MLE 可能不存在 |
| Wald-Test | Wald 检验 | 渐近正态 |
| Likelihood-Quotienten-Test | 似然比检验 | LRT |
| Modellguete | 模型拟合优度 | deviance, AIC |
| Sensitivitaet | 灵敏度 | $TP/(TP+FN)$ |
| Spezifitaet | 特异度 | $TN/(TN+FP)$ |
| ROC-Kurve | ROC 曲线 | TPR vs FPR |
| AUC | 曲线下面积 | discrimination |
| Fall-Kontroll-Studie | 病例-对照研究 | case-control |
| Propensity Score | 倾向得分 | $P(T=1\mid X)$ |
