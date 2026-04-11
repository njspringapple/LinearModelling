# 第 1 章学生笔记：导论与例子 / Kapitel 1: Einführung und Beispiele

材料 / Grundlage:

- 讲义 / Folien: `LinearModellingPPT.pdf`, Kapitel 1 `Einführung und Beispiele`, slides 7-16
- 主教材 / Hauptlehrbuch: Fahrmeir, Kneib, Lang, Marx, `Regression: Models, Methods and Applications`, Chapter 1 `Introduction`
- 相关主教材预告 / Überblick aus dem Hauptlehrbuch: Chapter 2.1 `Introduction` 和 2.10 `Models in a Nutshell`

## 覆盖审核 / Abdeckungsprüfung

| 来源 / Quelle | 内容点 / Inhaltspunkt | 笔记位置 / Wo in diesen Notizen |
|---|---|---|
| 讲义 1 | Galton 与 regression toward the mean | 1.1 |
| 讲义 1 | 回归概念：$Y$ 与 $x$ 的关系 | 1.1, 1.2 |
| 讲义 1 | Midlife crisis: 多元线性回归与 splines | 1.3 |
| 讲义 1 | 环境区 / Umweltzone: 协方差分析 | 1.3 |
| 讲义 1 | 教育例子：班级效应与 mixed model | 1.3 |
| 讲义 1 | 医学/流行病学：二元结果与 logistic model | 1.3 |
| 讲义 1 | 课程目标：理解、检验、预测、supervised learning | 1.2 |
| 主教材 Ch. 1 | 回归史、Galton 数据和回归线 | 1.1 |
| 主教材 Ch. 1.1 | 应用场景：租金、营养、医学、经济等 | 1.3 |
| 主教材 Ch. 1.2 | 初步数据分析：单变量分布、散点图、分组图 | 1.4 |
| 主教材 Ch. 1.3 | 记号约定：随机变量与观测值 | 1.5 |
| 主教材 Ch. 2 概览 | 线性、logit、mixed、additive、quantile 等模型族 | 1.3 |

审核结论 / Ergebnis: 第 1 章讲义内容已逐项覆盖；主教材 Chapter 1 中和本讲义最相关的历史、应用、EDA 和记号约定也已补入。主教材中大量具体案例细节没有逐字展开，因为它们是后续章节的例子素材；本笔记保留其建模用途与模型选择线索。

## 1.1 回归从哪里来：Galton 与均值回归 / Herkunft: Galton und Regression zur Mitte

**遇到的问题 / Problem.**
我们想研究孩子身高 $Y$ 与父母平均身高 $x$ 的关系。直觉上父母高，孩子也倾向于高；但这个关系并不是 $Y=x$，也不是确定性函数。Galton 发现极端高或极端矮的父母，其孩子平均上会更靠近总体均值，这就是 regression toward the mean。

Auf Deutsch: Galton untersuchte den Zusammenhang zwischen der Körpergröße erwachsener Kinder $Y$ und der mittleren Körpergröße der Eltern $x$. Der Zusammenhang zeigt eine Tendenz zur Mitte, also `Regression zur Mitte`.

**可能的解决路径 / Mögliche Wege.**

- 只列联表 / Kontingenztabelle: 能看频数，但不容易总结平均趋势。
- 只画散点图 / Streudiagramm: 能看方向，但不能量化趋势。
- 手画趋势线 / visuelle Trendlinie: 有直觉，但标准不统一。
- 建立回归模型 / Regressionsmodell: 用一个可估计、可解释的函数描述平均趋势。

**本课采用的方法及好处 / Methode der Vorlesung und Vorteil.**
把条件均值写成一条线：

$$
E(Y \mid x)=\beta_0+\beta_1 x
$$

Galton 例子中的斜率小于 $1$，直观含义是：父母身高偏离均值很远时，孩子的平均身高偏离均值较小。回归模型的好处是把“平均趋势”从随机波动中抽出来。

**具体解决 / Konkrete Lösung.**

1. 定义目标变量 / Zielgröße: $Y=$ 孩子身高。
2. 定义解释变量 / Einflussgröße: $x=$ 父母平均身高。
3. 画散点图，判断平均趋势是否近似线性。
4. 用 $\beta_0+\beta_1x$ 表示系统部分。
5. 用 $\varepsilon$ 表示未解释的个体差异。

**R 练习 / R-Übung.** R 内置没有 Galton 原始数据，这里用 `father.son` 不可靠，因为它不一定随 base R 安装。为了保证可操作，使用内置 `cars` 演示同一个思想：速度 $x$ 与刹车距离 $Y$ 的平均趋势。

```r
data(cars)
plot(cars[["speed"]], cars[["dist"]],
     xlab = "speed",
     ylab = "stopping distance",
     main = "Average trend plus random scatter")
fit <- lm(dist ~ speed, data = cars)
abline(fit, col = "red", lwd = 2)
coef(fit)
```

练习 / Übung: 红线是否表示每辆车的真实刹车距离？
答案 / Antwort: 不是。红线表示 $\widehat{E}(Y \mid x)$，即给定速度时刹车距离的估计平均趋势。

## 1.2 回归到底要解决什么：理解、检验、预测 / Ziele: Verstehen, Nachweisen, Prognostizieren

**遇到的问题 / Problem.**
面对同一个数据集，我们可能问三类问题：

- 哪些变量影响 $Y$？影响方向和大小如何？
- 某个影响是否有统计证据？
- 给定新个体 $x_0$，如何预测 $Y_0$？

Auf Deutsch: Dieselben Daten können für Erklärung, Inferenz oder Prognose verwendet werden.

**可能的解决路径 / Mögliche Wege.**

- 解释 / Verstehen: 关注参数 $\beta_j$ 的符号、大小、单位和可解释性。
- 检验 / Nachweisen: 关注假设如 $H_0:\beta_j=0$，以及置信区间和检验。
- 预测 / Prognostizieren: 关注 $\widehat{Y}_0$、预测误差和预测区间。

**本课采用的方法及好处 / Methode der Vorlesung und Vorteil.**
讲义把回归理解为 supervised learning：从观测对 $(Y_i,x_i)$ 学习 $x$ 与 $Y$ 的关系。形式上写为：

$$
Y = f(x_1,\dots,x_k)+\varepsilon
$$

或：

$$
E(Y \mid x_1,\dots,x_k)=f(x_1,\dots,x_k)
$$

好处是同一个模型可以同时服务解释、推断和预测，但不同目标会影响我们对模型好坏的判断。例如解释型模型强调可解释性；预测型模型强调泛化误差。

**具体解决 / Konkrete Lösung.**

- 先写研究问题：解释 / 检验 / 预测？
- 再写变量表：$Y$ 是什么，$x_j$ 是什么？
- 最后选择模型：线性、logistic、mixed、spline 等。

**R 练习 / R-Übung.**

```r
data(mtcars)
fit <- lm(mpg ~ wt, data = mtcars)
summary(fit)
confint(fit)
predict(fit, newdata = data.frame(wt = c(2.5, 3.5)),
        interval = "prediction")
```

直觉问题 / Intuitionsfrage: `summary(fit)` 和 `predict()` 分别更像在服务什么目标？
答案 / Antwort: `summary(fit)` 更偏解释与检验；`predict()` 更偏预测。

## 1.3 讲义中的应用例子与模型选择 / Beispiele und Modellwahl

**遇到的问题 / Problem.**
不同问题的数据结构不同。如果用同一个普通线性模型硬套所有场景，会误解问题。例如二元疾病状态不能直接当成正态连续响应来处理；班级内学生并不一定独立。

Auf Deutsch: Die Wahl der Modellklasse hängt von Zielgröße, Kovariaten und Datenstruktur ab.

**可能的解决路径 / Mögliche Wege.**

| 讲义例子 / Beispiel | 数据特征 / Datenstruktur | 可能模型 / Mögliches Modell |
|---|---|---|
| Galton 身高 | 连续 $Y$，连续 $x$ | 简单线性回归 |
| Midlife crisis | 连续满意度，多解释变量，年龄非线性 | 多元线性回归 + splines |
| Umweltzone | 连续污染物，政策阶段为分组变量，温度作协变量 | 协方差分析 / ANCOVA |
| 教育阅读能力 | 学生嵌套在班级中 | 混合线性模型 |
| 慢性支气管炎 CBR | $Y\in\{0,1\}$ | logistic 回归 |

**本课采用的方法及好处 / Methode der Vorlesung und Vorteil.**
教授先讲线性模型，然后逐步扩展。主教材 Chapter 2 也用统一框架介绍模型族：线性模型、logit 模型、Poisson 模型、GLM、mixed model、additive model、structured additive regression、quantile regression。这样后续模型都能看成对以下三部分之一的扩展：

- 目标变量 $Y$ 的类型变化；
- 误差项 $\varepsilon$ 的结构变化；
- 系统函数 $f(x)$ 的灵活性变化。

**具体解决 / Konkrete Lösung.**

如果 $Y$ 是连续的、误差大致独立同方差、平均关系近似线性，先用：

$$
Y_i=\beta_0+\beta_1x_{i1}+\dots+\beta_kx_{ik}+\varepsilon_i
$$

如果 $Y$ 是二元变量，转向：

$$
P(Y_i=1)=G(\beta_0+\beta_1x_{i1}+\dots+\beta_kx_{ik})
$$

其中 logistic 模型常用：

$$
G(\eta)=\frac{\exp(\eta)}{1+\exp(\eta)}
$$

如果有班级、个体、地区等聚类结构，就需要混合模型，例如：

$$
Y_{ij}=\beta_0+\beta_1x_{ij}+\alpha_j+\varepsilon_{ij}
$$

其中 $\alpha_j$ 表示第 $j$ 个班级或群组的随机效应。

**R 练习 / R-Übung.** 用内置 `iris` 展示不同响应类型导致不同模型。

```r
data(iris)

# 连续响应 + 连续解释变量
fit_lm <- lm(Sepal.Length ~ Petal.Length, data = iris)
summary(fit_lm)

# 连续响应 + 分类解释变量
fit_group <- lm(Sepal.Length ~ Species, data = iris)
summary(fit_group)

# 二元响应：是否 setosa
iris[["is_setosa"]] <- ifelse(iris[["Species"]] == "setosa", 1, 0)
fit_logit <- glm(is_setosa ~ Petal.Length,
                 data = iris,
                 family = binomial)
summary(fit_logit)
```

练习 / Übung: 为什么第三个模型不用 `lm()` 作为首选？
答案 / Antwort: 因为 $Y$ 是 $0/1$，条件均值是概率，必须在 $[0,1]$ 内。logistic 回归通过 $G(\eta)$ 把线性预测子转换成概率。

## 1.4 主教材补充：建模前的 first steps / Ergänzung: Erste Schritte vor dem Modell

**遇到的问题 / Problem.**
讲义第 1 章直接进入模型例子，但主教材强调：建模前要先看数据。否则我们可能在变量类型、异常值、非线性或分组结构还没理解时就套模型。

Auf Deutsch: Vor der Modellierung steht die explorative Datenanalyse.

**可能的解决路径 / Mögliche Wege.**

- 单变量分布 / univariate Verteilungen: 看 $Y$ 和每个 $x$ 的分布、范围、异常值。
- 连续变量关系 / grafische Assoziationsanalyse: 用散点图、平滑趋势、分组均值图。
- 分类变量关系 / kategoriale Kovariaten: 用箱线图、分组均值、密度图。
- 大样本散点图过密时 / große Stichprobe: 用分箱均值或分组标准差简化图像。

**本课采用的方法及好处 / Methode der Vorlesung und Vorteil.**
主教材的策略是先用图形理解变量，再决定模型族。这能提前发现三类问题：

- 关系不是线性：需要变换、polynomial 或 spline；
- 方差不恒定：后续可用 WLS/GLS；
- 分组明显：可能需要 dummy coding、ANOVA/ANCOVA 或 mixed model。

**具体解决 / Konkrete Lösung.**

建模前至少做：

1. `summary(data)` 检查范围和缺失。
2. `hist(Y)` 或 `boxplot(Y)` 看目标变量分布。
3. `plot(x, Y)` 看连续解释变量。
4. `boxplot(Y ~ group)` 看分类解释变量。

**R 练习 / R-Übung.** 用 `mtcars` 做一套建模前检查。

```r
data(mtcars)
summary(mtcars[, c("mpg", "wt", "hp", "cyl")])

hist(mtcars[["mpg"]], main = "Distribution of mpg", xlab = "mpg")
plot(mtcars[["wt"]], mtcars[["mpg"]],
     xlab = "weight",
     ylab = "mpg",
     main = "mpg vs weight")
boxplot(mpg ~ cyl, data = mtcars,
        xlab = "number of cylinders",
        ylab = "mpg",
        main = "mpg by cylinders")
```

直觉问题 / Intuitionsfrage: `boxplot(mpg ~ cyl)` 对应主教材中的哪类图形思想？
答案 / Antwort: 这是连续响应变量 $Y$ 对分类解释变量的分组比较。

## 1.5 主教材补充：记号约定 / Notation

**遇到的问题 / Problem.**
统计书有时用大写 $Y$ 表示随机变量，用小写 $y$ 表示观测值；回归书有时为简洁把二者都写成 $y$。如果不注意，容易混淆“随机量”和“已经观测到的数据”。

Auf Deutsch: Man muss zwischen Zufallsvariablen und Realisationen unterscheiden, auch wenn die Notation manchmal vereinfacht wird.

**可能的解决路径 / Mögliche Wege.**

- 严格写法：随机变量用 $Y$，观测值用 $y$。
- 回归教材常用写法：都写成 $y$，由上下文判断。
- 学生笔记建议：推导时保持 $Y_i$，实际数据和 R 输出可写 $y_i$。

**本课采用的方法及好处 / Methode der Vorlesung und Vorteil.**
本笔记采用折中：公式推导里用 $Y_i$ 强调随机性，用 $x_i$ 表示给定解释变量；R 代码和实际数据中用小写变量名。这有助于理解为什么同一个模型能谈期望、方差和置信区间。

**具体解决 / Konkrete Lösung.**

$$
Y_i = \text{random variable / Zufallsvariable}
$$

$$
y_i = \text{observed value / Realisation}
$$

$$
x_i = \text{fixed observed covariate / beobachtete Einflussgröße}
$$

## 1.6 选择题 / Multiple-Choice-Fragen

1. 回归分析中 $E(Y \mid x)$ 表示什么？
   A. 给定 $x$ 时 $Y$ 的平均趋势
   B. 每个观测点的精确值
   C. $x$ 的方差
   D. 样本量
   答案 / Antwort: A

2. 如果 $Y$ 是是否患病，取值为 $0/1$，最自然的模型起点是？
   A. 普通线性回归
   B. logistic 回归
   C. 只画直方图
   D. 不需要模型
   答案 / Antwort: B

3. 如果学生嵌套在班级中，且班级之间有差异，讲义建议的方向是？
   A. mixed model
   B. 只计算总体均值
   C. 忽略班级
   D. 只做直方图
   答案 / Antwort: A

4. 建模前的 exploratory data analysis 主要目的是什么？
   A. 替代所有模型
   B. 先理解变量分布、异常值和关系形状
   C. 保证 $R^2=1$
   D. 删除所有分类变量
   答案 / Antwort: B

5. `Zielgröße` 对应哪个概念？
   A. 目标变量 / response variable
   B. 残差
   C. 方差
   D. 样本编号
   答案 / Antwort: A

## 1.7 本章总结 / Zusammenfassung

中文总结：

- 回归研究目标变量 $Y$ 与解释变量 $x$ 的关系。
- 关系一般不是确定性的，所以写成 $Y=f(x)+\varepsilon$。
- 本课程关注三类目标：理解、检验、预测。
- 讲义例子展示了不同模型族：线性回归、协方差分析、混合模型、logistic 回归、splines。
- 主教材补充强调：正式建模前要做探索性数据分析，包括单变量分布、散点图、分组图和异常值检查。
- 记号上要区分随机变量 $Y_i$ 与观测值 $y_i$，虽然回归教材有时会简写。

Deutsche Zusammenfassung:

- Regression untersucht den Zusammenhang zwischen einer Zielgröße $Y$ und Einflussgrößen $x$.
- Der Zusammenhang wird als $Y=f(x)+\varepsilon$ mit systematischer und stochastischer Komponente formuliert.
- Die wichtigsten Ziele sind Verstehen, Nachweisen und Prognostizieren.
- Unterschiedliche Datenstrukturen führen zu unterschiedlichen Modellklassen.
- Vor der Modellierung ist explorative Datenanalyse wichtig.
- Die Notation unterscheidet konzeptionell zwischen Zufallsvariable $Y_i$ und beobachtetem Wert $y_i$.

## 1.8 专业德语单词汇总 / Fachwortschatz

| Deutsch | 中文 | 公式或说明 |
|---|---|---|
| Zielgröße | 目标变量 | $Y$ |
| Einflussgröße | 解释变量 / 影响变量 | $x$ |
| abhängige Variable | 因变量 | $Y$ |
| unabhängige Variable | 自变量 | $x$ |
| Kovariate | 协变量 | $x_j$ |
| Regressor | 回归变量 | $x_j$ |
| bedingter Erwartungswert | 条件期望 | $E(Y \mid x)$ |
| systematische Komponente | 系统部分 | $f(x)$ |
| stochastische Komponente | 随机部分 | $\varepsilon$ |
| Störterm / Fehlerterm | 扰动项 / 误差项 | $\varepsilon$ |
| Regression zur Mitte | 均值回归 | Galton 例子 |
| Streudiagramm | 散点图 | plot of $Y$ versus $x$ |
| explorative Datenanalyse | 探索性数据分析 | 建模前检查 |
| Verstehen | 理解 | 参数解释 |
| Nachweisen | 证明 / 检验 | $H_0:\beta_j=0$ |
| Prognostizieren | 预测 | $\widehat{Y}_0$ |
| supervised learning | 监督学习 | 从 $(Y_i,x_i)$ 学习关系 |
| logistisches Regressionsmodell | logistic 回归模型 | $P(Y=1)=G(\eta)$ |
| gemischtes lineares Modell | 混合线性模型 | $Y_{ij}=\beta_0+\beta_1x_{ij}+\alpha_j+\varepsilon_{ij}$ |
| Kovarianzanalyse | 协方差分析 | 连续 $Y$ + 分组变量 + 协变量 |
