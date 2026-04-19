# 第 11 章学生笔记：变量选择 / Kapitel 11: Variablenselektion

材料 / Grundlage:

- 讲义 / Folien: `chapter_pdfs/LinearModellingPPT/11_Variablenselektion.pdf`
- 主教材 / Hauptlehrbuch: `Regression.pdf`, Chapter 3.4, Chapter 4.2, Chapter 4.3
- 辅助教材 / Hilfslehrbuch: Faraway, `Linear Models with R`, Chapter 8, Chapter 9, Chapter 10
- 扩展教材 / Erweiterung: Faraway, `Extending The Linear Model With R`, Chapter 2.5, Appendix A.3, Chapter 16-17 as outlook

## 覆盖审核 / Abdeckungspruefung

| 来源 / Quelle | 内容点 / Inhaltspunkt | 笔记位置 / Wo |
|---|---|---|
| 讲义 11 | 建模目标：解释 vs 预测 | 11.1 |
| 讲义 11 | 基本原则、过拟合、模型不确定性 | 11.2 |
| 讲义 11 | 模型优度指标 | 11.3 |
| 讲义 11 | $R^2$、调整 $R^2$、AIC、BIC | 11.3 |
| 讲义 11 | 交叉验证 | 11.4 |
| 讲义 11 | 前向、后向、逐步选择 | 11.5 |
| 讲义 11 | 变量选择后的推断后果 | 11.6 |
| 讲义 11 | explanatory vs predictive modeling | 11.7 |
| 讲义 11 | LASSO、boosting、树和神经网络展望 | 11.8 |
| 主教材 | model choice、regularization、boosting | 全章 |
| Faraway | variable selection、shrinkage、model uncertainty | 全章 |

## 11.1 变量选择首先取决于目标

**遇到的问题 / Problem.**  
“哪些变量应该进模型？”没有唯一答案。解释型建模和预测型建模的目标不同，变量选择标准也不同。

**可能的解决路径 / Moegliche Wege.**

- 解释型 / explanatory: 关注因果或科学解释，变量应由研究问题和混杂控制决定。
- 预测型 / predictive: 关注新数据预测误差，变量可为预测性能服务。
- 描述型 / descriptive: 关注总结数据结构。

**本课采用的方法及好处 / Methode und Vorteil.**  
讲义强调先明确目标，再选择模型。好处是避免用预测标准回答因果问题，或用单个 $p$ 值决定预测模型。

**R 练习 / R-Uebung.**

```r
data(mtcars)
fit_explain <- lm(mpg ~ wt + hp + cyl, data = mtcars)
fit_predict <- lm(mpg ~ ., data = mtcars)

summary(fit_explain)
summary(fit_predict)
```

## 11.2 过拟合、偏差-方差权衡与模型不确定性

**遇到的问题 / Problem.**  
加入更多变量通常能改善训练集拟合，但可能降低新数据预测能力。变量选择本身也会引入不确定性。

**本课采用的方法及好处 / Methode und Vorteil.**  
偏差-方差权衡：

- 模型太简单：高偏差，欠拟合。
- 模型太复杂：高方差，过拟合。
- 好模型：在偏差和方差之间折中。

模型不确定性：如果很多模型都差不多好，只报告一个“最终模型”会低估不确定性。

**具体解决 / Konkrete Loesung.**  
不要只看训练集 $R^2$；要结合调整 $R^2$、AIC/BIC、交叉验证和领域知识。

## 11.3 模型优度指标 / Masse fuer die Modellguete

**遇到的问题 / Problem.**  
不同指标偏好不同模型。需要知道每个指标在奖励拟合和惩罚复杂度之间如何取舍。

**本课采用的方法及好处 / Methode und Vorteil.**

普通 $R^2$：

$$
R^2=1-\frac{SSE}{SST}
$$

加入变量后 $R^2$ 不会下降，因此容易鼓励复杂模型。

调整 $R^2$：

$$
R^2_{adj}=1-\frac{SSE/(n-p-1)}{SST/(n-1)}
$$

AIC:

$$
AIC=-2\ell(\widehat{\theta})+2k
$$

BIC:

$$
BIC=-2\ell(\widehat{\theta})+\log(n)k
$$

BIC 对复杂度惩罚通常比 AIC 更强。

**R 练习 / R-Uebung.**

```r
data(mtcars)
m1 <- lm(mpg ~ wt, data = mtcars)
m2 <- lm(mpg ~ wt + hp + cyl, data = mtcars)
m3 <- lm(mpg ~ ., data = mtcars)

data.frame(
  model = c("m1", "m2", "m3"),
  r2 = c(summary(m1)[["r.squared"]],
         summary(m2)[["r.squared"]],
         summary(m3)[["r.squared"]]),
  adj_r2 = c(summary(m1)[["adj.r.squared"]],
             summary(m2)[["adj.r.squared"]],
             summary(m3)[["adj.r.squared"]]),
  AIC = c(AIC(m1), AIC(m2), AIC(m3)),
  BIC = c(BIC(m1), BIC(m2), BIC(m3))
)
```

## 11.4 交叉验证 / Kreuzvalidierung

**遇到的问题 / Problem.**  
训练误差低不代表新数据预测好。需要估计泛化误差。

**本课采用的方法及好处 / Methode und Vorteil.**  
$K$ 折交叉验证：

1. 把数据分成 $K$ 份。
2. 每次留一份作测试集，其余作训练集。
3. 计算测试误差。
4. 对 $K$ 次测试误差取平均。

好处是：更接近预测型目标。

**R 练习 / R-Uebung.** 手写 5-fold CV。

```r
set.seed(1)
data(mtcars)
K <- 5
fold <- sample(rep(seq_len(K), length.out = nrow(mtcars)))

cv_mse <- function(formula) {
  errs <- numeric(K)
  for (k in seq_len(K)) {
    train <- mtcars[fold != k, ]
    test <- mtcars[fold == k, ]
    fit <- lm(formula, data = train)
    pred <- predict(fit, newdata = test)
    errs[k] <- mean((test[["mpg"]] - pred)^2)
  }
  mean(errs)
}

c(wt = cv_mse(mpg ~ wt),
  wt_hp_cyl = cv_mse(mpg ~ wt + hp + cyl),
  all = cv_mse(mpg ~ .))
```

## 11.5 前向、后向和逐步选择

**遇到的问题 / Problem.**  
候选变量很多时，所有模型枚举成本高。逐步方法提供自动搜索，但也有风险。

**本课采用的方法及好处 / Methode und Vorteil.**

- Forward selection: 从空模型开始逐步加入变量。
- Backward elimination: 从完整模型开始逐步删除变量。
- Stepwise selection: 加入和删除都允许。

常用标准是 AIC 或 BIC。好处是快速；缺点是搜索路径依赖强，选择后的 $p$ 值不再像预先指定模型那样解释。

**R 练习 / R-Uebung.**

```r
data(mtcars)
full <- lm(mpg ~ ., data = mtcars)
null <- lm(mpg ~ 1, data = mtcars)

step_aic <- step(null,
                 scope = list(lower = null, upper = full),
                 direction = "both",
                 trace = 0)
formula(step_aic)
AIC(step_aic)
```

## 11.6 变量选择后的推断风险

**遇到的问题 / Problem.**  
如果先看数据选择模型，再把最终模型当作事先指定的模型来做 $p$ 值解释，会低估不确定性。

**本课采用的方法及好处 / Methode und Vorteil.**  
变量选择会带来：

- $p$ 值偏小；
- 置信区间过窄；
- 系数估计偏向显著结果；
- 模型不确定性被忽略。

**具体解决 / Konkrete Loesung.**

- 解释型研究应预先指定关键变量和混杂因素。
- 预测型研究应使用独立测试集或交叉验证。
- 报告选择过程，不只报告最终模型。
- 可考虑模型平均、bootstrap 或正则化。

## 11.7 explanatory vs predictive

**遇到的问题 / Problem.**  
解释模型和预测模型经常被混用。一个变量对预测有用，不代表它有因果解释；一个因果重要变量也可能预测增益不大。

**本课采用的方法及好处 / Methode und Vorteil.**

| 目标 | 关注点 | 常见工具 |
|---|---|---|
| explanatory | 参数解释、混杂控制、研究设计 | 预先设定模型、线性假设检验 |
| predictive | 新数据误差、泛化能力 | CV、测试集、AIC、regularization |

**具体解决 / Konkrete Loesung.**  
写报告时先说“这个模型用于解释还是预测”。不要用 AIC 选择出来的模型直接宣称因果关系。

## 11.8 正则化、LASSO、boosting 与展望

**遇到的问题 / Problem.**  
变量很多或共线性强时，逐步选择不稳定。正则化通过惩罚复杂度让估计更稳定。

**本课采用的方法及好处 / Methode und Vorteil.**  
Ridge:

$$
\min_\beta SSE+\lambda\sum_{j=1}^p\beta_j^2
$$

LASSO:

$$
\min_\beta SSE+\lambda\sum_{j=1}^p|\beta_j|
$$

LASSO 可以把部分系数压到 $0$，因此同时做 shrinkage 和变量选择。Boosting、trees 和 neural networks 是更偏预测的后续方向。

**具体解决 / Konkrete Loesung.**  
本课主线先理解 AIC/BIC/CV/逐步选择；正则化作为后续扩展，尤其适合高维、共线性和预测任务。

## 11.9 选择题 / Multiple-Choice-Fragen

1. 为什么普通 $R^2$ 不适合单独做变量选择？  
   A. 加变量后通常不会下降  
   B. 它总是负数  
   C. 它不依赖数据  
   D. 它只能用于 logistic 回归  
   答案 / Antwort: A

2. BIC 相比 AIC 通常怎样？  
   A. 对复杂度惩罚更强  
   B. 不惩罚变量数  
   C. 只看 $R^2$  
   D. 只能用于 ANOVA  
   答案 / Antwort: A

3. 交叉验证主要估计什么？  
   A. 泛化预测误差  
   B. 响应变量名字  
   C. 截距是否为零  
   D. 方差是否等于一  
   答案 / Antwort: A

4. 变量选择后的 $p$ 值有什么风险？  
   A. 可能过于乐观  
   B. 一定完全正确  
   C. 与选择过程无关  
   D. 不能计算  
   答案 / Antwort: A

5. LASSO 的惩罚项是什么？  
   A. $\lambda\sum|\beta_j|$  
   B. $\lambda\sum\beta_j^2$  
   C. $R^2$  
   D. AUC  
   答案 / Antwort: A

## 11.10 本章总结 / Zusammenfassung

中文总结：

- 变量选择必须先明确目标：解释还是预测。
- 普通 $R^2$ 会随变量增加而不降，不能单独作为选择标准。
- 调整 $R^2$、AIC、BIC 在拟合和复杂度之间做权衡。
- 交叉验证直接面向预测误差。
- 前向、后向和逐步选择方便但不稳定，会影响后续推断。
- 变量选择后的 $p$ 值和置信区间容易过于乐观。
- 解释型模型需要研究设计和混杂控制；预测型模型需要测试误差评估。
- Ridge、LASSO、boosting 是变量选择和预测建模的后续方向。

Deutsche Zusammenfassung:

- Variablenselektion haengt vom Modellierungsziel ab.
- $R^2$ allein ist kein gutes Selektionskriterium.
- AIC und BIC balancieren Anpassung und Komplexitaet.
- Kreuzvalidierung schaetzt die Vorhersageleistung.
- Stepwise-Verfahren sind praktisch, aber inferentiell riskant.
- LASSO kombiniert Shrinkage und Variablenselektion.

## 11.11 专业德语单词汇总 / Fachwortschatz

| Deutsch | 中文 | 说明 |
|---|---|---|
| Variablenselektion | 变量选择 | model selection |
| Modellwahl | 模型选择 | choice of model |
| Zielsetzung der Modellierung | 建模目标 | explanation vs prediction |
| Modellguete | 模型优度 | fit quality |
| angepasstes $R^2$ | 调整 $R^2$ | penalized for predictors |
| Akaike-Informationskriterium | AIC | $-2\ell+2k$ |
| Bayes-Informationskriterium | BIC | $-2\ell+\log(n)k$ |
| Kreuzvalidierung | 交叉验证 | CV |
| Vorwaertsauswahl | 前向选择 | forward selection |
| Rueckwaertselimination | 后向剔除 | backward elimination |
| schrittweise Auswahl | 逐步选择 | stepwise selection |
| Ueberanpassung | 过拟合 | overfitting |
| Modellunsicherheit | 模型不确定性 | selection uncertainty |
| Shrinkage | 收缩 | coefficient shrinkage |
| Ridge Regression | Ridge 回归 | $L_2$ penalty |
| LASSO | LASSO | $L_1$ penalty |
| Boosting | 提升法 | predictive method |
