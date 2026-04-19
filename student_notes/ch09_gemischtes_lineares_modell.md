# 第 9 章学生笔记：混合线性模型 / Kapitel 9: Das gemischte lineare Regressionsmodell

材料 / Grundlage:

- 讲义 / Folien: `chapter_pdfs/LinearModellingPPT/09_Das_gemischte_lineare_Regressionsmodell_Linear_mixed_Model.pdf`
- 主教材 / Hauptlehrbuch: `Regression.pdf`, Chapter 7
- 辅助教材 / Hilfslehrbuch: Faraway, `Linear Models with R`, Chapter 16 as nearby block/random-effect intuition
- 扩展教材 / Erweiterung: Faraway, `Extending The Linear Model With R`, Chapter 10, Chapter 11, Chapter 13 as advanced extension

## 覆盖审核 / Abdeckungspruefung

| 来源 / Quelle | 内容点 / Inhaltspunkt | 笔记位置 / Wo |
|---|---|---|
| 讲义 9 | 阅读促进研究例子 | 9.1 |
| 讲义 9 | 简单随机效应模型 | 9.2 |
| 讲义 9 | 边际模型 | 9.3 |
| 讲义 9 | 纵向数据和重复测量 | 9.4 |
| 讲义 9 | 表示为一般线性模型 | 9.5 |
| 讲义 9 | 混合模型推断 | 9.6 |
| 讲义 9 | 结果解释 | 9.7 |
| 讲义 9 | 固定效应模型 vs 混合模型 | 9.8 |
| 主教材 | LMM、clustered/longitudinal data、likelihood inference | 全章 |
| Extending Faraway | random effects、repeated measures、GLMM preview | 9.2-9.8 |

## 9.1 为什么需要混合模型

**遇到的问题 / Problem.**  
普通线性模型假设观测独立。但教育、医学和纵向数据中，观测常常按班级、学校、病人、地区或个体聚类。同一组内观测更相似，不能当作完全独立。

**可能的解决路径 / Moegliche Wege.**

- 忽略聚类：标准误可能偏小。
- 把所有组都当固定效应：可控制组差异，但参数很多，难以推广。
- 使用随机效应：把组差异建模为随机变量，得到混合模型。

**本课采用的方法及好处 / Methode und Vorteil.**  
混合线性模型同时包含：

- 固定效应 / fixed effects: 研究者关心的平均效应；
- 随机效应 / random effects: 聚类单位的随机偏差。

好处是：既能建模组内相关，又能估计总体平均效应。

## 9.2 简单随机截距模型 / Einfacher zufaelliger Effekt

**遇到的问题 / Problem.**  
不同班级或个体可能有不同基础水平。如果只放一个共同截距，会忽略组间差异。

**本课采用的方法及好处 / Methode und Vorteil.**  
随机截距模型：

$$
Y_{ij}=\beta_0+\beta_1x_{ij}+b_j+\varepsilon_{ij}
$$

其中：

$$
b_j\sim N(0,\tau^2)
$$

$$
\varepsilon_{ij}\sim N(0,\sigma^2)
$$

$b_j$ 表示第 $j$ 个组相对总体截距的随机偏差。好处是：同一组内观测共享 $b_j$，因此自然相关。

**具体解决 / Konkrete Loesung.**  
解释时：

- $\beta_0,\beta_1$ 是固定效应，描述总体平均关系。
- $b_j$ 是随机效应，描述组 $j$ 的偏离。
- $\tau^2$ 是组间方差。
- $\sigma^2$ 是组内残差方差。

**R 练习 / R-Uebung.** 使用 R 内置 `CO2` 数据，若 `nlme` 可用则拟合随机截距模型。

```r
data(CO2)
fit_lm <- lm(uptake ~ conc + Treatment + Plant, data = CO2)
summary(fit_lm)

if (requireNamespace("nlme", quietly = TRUE)) {
  fit_lme <- nlme::lme(uptake ~ conc + Treatment,
                       random = ~ 1 | Plant,
                       data = CO2)
  summary(fit_lme)
}
```

## 9.3 边际模型 / Marginales Modell

**遇到的问题 / Problem.**  
随机效应模型可从条件角度看，也可从边际协方差角度看。我们需要理解它如何导致组内相关。

**本课采用的方法及好处 / Methode und Vorteil.**  
条件模型：

$$
Y_{ij}\mid b_j=\beta_0+\beta_1x_{ij}+b_j+\varepsilon_{ij}
$$

边际均值：

$$
E(Y_{ij})=\beta_0+\beta_1x_{ij}
$$

同组内协方差：

$$
\operatorname{Cov}(Y_{ij},Y_{ik})=\tau^2
$$

组内相关系数：

$$
\rho=\frac{\tau^2}{\tau^2+\sigma^2}
$$

好处是：随机截距模型等价于一种特殊的相关结构。

## 9.4 纵向数据和重复测量 / Longitudinale Daten

**遇到的问题 / Problem.**  
同一个个体被重复测量时，观测不独立；而且每个个体可能有不同起点，甚至不同时间趋势。

**本课采用的方法及好处 / Methode und Vorteil.**  
随机截距 + 随机斜率模型：

$$
Y_{ij}=\beta_0+\beta_1t_{ij}+b_{0j}+b_{1j}t_{ij}+\varepsilon_{ij}
$$

其中 $b_{0j}$ 表示个体 $j$ 的随机起点，$b_{1j}$ 表示个体 $j$ 的随机趋势。

好处是：既允许平均趋势，又允许个体差异。

**R 练习 / R-Uebung.**

```r
data(Orange)
plot(Orange[["age"]], Orange[["circumference"]],
     col = as.numeric(Orange[["Tree"]]),
     xlab = "age",
     ylab = "circumference")

if (requireNamespace("nlme", quietly = TRUE)) {
  fit_orange <- nlme::lme(circumference ~ age,
                          random = ~ age | Tree,
                          data = Orange)
  summary(fit_orange)
}
```

## 9.5 表示为一般线性模型

**遇到的问题 / Problem.**  
混合模型看起来比第 8 章复杂，但从边际角度看，它仍然是一个带特殊协方差矩阵的一般线性模型。

**本课采用的方法及好处 / Methode und Vorteil.**  
矩阵形式：

$$
Y=X\beta+Zb+\varepsilon
$$

其中：

- $X$: 固定效应设计矩阵。
- $Z$: 随机效应设计矩阵。
- $b$: 随机效应向量。

通常假设：

$$
b\sim N(0,D)
$$

$$
\varepsilon\sim N(0,\sigma^2I)
$$

边际方差：

$$
V(Y)=ZDZ^T+\sigma^2I
$$

好处是：混合模型连接了第 8 章的 GLS 思想和分层数据结构。

## 9.6 推断：ML 与 REML

**遇到的问题 / Problem.**  
混合模型不仅要估计固定效应 $\beta$，还要估计方差分量，如 $\tau^2$ 和 $\sigma^2$。

**可能的解决路径 / Moegliche Wege.**

- ML：最大化完整似然。
- REML：在消除固定效应影响后估计方差分量。

**本课采用的方法及好处 / Methode und Vorteil.**  
REML 常用于方差分量估计，因为它考虑固定效应估计带来的自由度损失。ML 更适合比较固定效应结构不同的模型时使用同一框架，但模型比较要谨慎。

**具体解决 / Konkrete Loesung.**  
报告混合模型时至少说明：

- 固定效应估计。
- 随机效应结构。
- 方差分量。
- 使用 ML 还是 REML。
- 聚类层级和样本量。

## 9.7 结果解释 / Interpretation der Ergebnisse

**遇到的问题 / Problem.**  
混合模型输出比 `lm()` 更复杂，容易混淆固定效应、随机效应和方差分量。

**本课采用的方法及好处 / Methode und Vorteil.**

- 固定效应解释总体平均关系。
- 随机截距方差解释组间基线差异。
- 残差方差解释组内未解释波动。
- ICC 解释组内相似程度。

**具体解决 / Konkrete Loesung.**  
例如随机截距模型中，若 $\tau^2$ 很大，说明组之间基础水平差异明显；若 $\rho$ 高，说明同组内观测高度相似，不能忽略聚类。

## 9.8 固定效应模型 vs 混合模型

**遇到的问题 / Problem.**  
面对组别变量，可以把组别作为固定效应，也可以作为随机效应。选择不同，解释不同。

**可能的解决路径 / Moegliche Wege.**

- 固定效应模型：关心样本中这些特定组。
- 随机效应模型：把组看作来自更大总体的随机样本。

**本课采用的方法及好处 / Methode und Vorteil.**  
混合模型适合组数较多、组被看作随机抽样、并希望估计组间方差的场景。固定效应适合少数明确感兴趣的组。

**R 练习 / R-Uebung.**

```r
data(CO2)
fit_fixed_group <- lm(uptake ~ conc + Treatment + Plant, data = CO2)
fit_no_group <- lm(uptake ~ conc + Treatment, data = CO2)

anova(fit_no_group, fit_fixed_group)
```

## 9.9 选择题 / Multiple-Choice-Fragen

1. 混合模型中的 fixed effects 表示什么？  
   A. 总体平均效应  
   B. 随机噪声编号  
   C. 样本顺序  
   D. 只能是分类变量  
   答案 / Antwort: A

2. 随机截距 $b_j$ 表示什么？  
   A. 第 $j$ 组相对总体截距的偏差  
   B. 第 $j$ 个残差平方  
   C. 第 $j$ 个变量的单位  
   D. 所有组的共同斜率  
   答案 / Antwort: A

3. 随机截距模型为什么会产生组内相关？  
   A. 同组观测共享同一个 $b_j$  
   B. 因为所有 $Y$ 都相同  
   C. 因为没有误差项  
   D. 因为 $X$ 没有截距  
   答案 / Antwort: A

4. REML 常用于什么？  
   A. 方差分量估计  
   B. ROC 曲线  
   C. dummy coding  
   D. Bonferroni 校正  
   答案 / Antwort: A

5. 若组是从更大总体随机抽取的，常考虑什么模型？  
   A. 随机效应或混合模型  
   B. 只含截距的普通模型  
   C. 删除所有组信息  
   D. 必须 logistic 回归  
   答案 / Antwort: A

## 9.10 本章总结 / Zusammenfassung

中文总结：

- 混合线性模型用于聚类、纵向和重复测量数据。
- 模型同时包含固定效应和随机效应。
- 随机截距模型允许不同组有不同基线水平。
- 边际上，随机效应会产生组内相关。
- 一般矩阵形式为 $Y=X\beta+Zb+\varepsilon$。
- 方差分量通常用 ML 或 REML 估计。
- 固定效应模型和混合模型的选择取决于研究问题和组别的抽样含义。

Deutsche Zusammenfassung:

- Lineare gemischte Modelle behandeln gruppierte und longitudinale Daten.
- Feste Effekte beschreiben mittlere Effekte, zufaellige Effekte beschreiben gruppenspezifische Abweichungen.
- Ein zufaelliger Achsenabschnitt erzeugt Korrelation innerhalb derselben Gruppe.
- Die allgemeine Form lautet $Y=X\beta+Zb+\varepsilon$.
- Varianzkomponenten werden typischerweise mit ML oder REML geschaetzt.

## 9.11 专业德语单词汇总 / Fachwortschatz

| Deutsch | 中文 | 说明 |
|---|---|---|
| gemischtes lineares Modell | 混合线性模型 | LMM |
| fester Effekt | 固定效应 | population average |
| zufaelliger Effekt | 随机效应 | group deviation |
| zufaelliger Achsenabschnitt | 随机截距 | random intercept |
| zufaellige Steigung | 随机斜率 | random slope |
| marginales Modell | 边际模型 | integrated over random effects |
| longitudinale Daten | 纵向数据 | repeated over time |
| Wiederholungsmessungen | 重复测量 | repeated measures |
| Varianzkomponente | 方差分量 | $\tau^2$, $\sigma^2$ |
| Intraklassenkorrelation | 组内相关 | ICC |
| ML-Schaetzung | 最大似然估计 | likelihood |
| REML-Schaetzung | 限制最大似然 | variance components |
| Cluster | 聚类/组 | class, subject, school |
