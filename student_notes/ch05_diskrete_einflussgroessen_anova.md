# 第 5 章学生笔记：离散解释变量、编码与方差分析 / Kapitel 5: Diskrete Einflussgroessen

材料 / Grundlage:

- 讲义 / Folien: `chapter_pdfs/LinearModellingPPT/05_Diskrete_Einflugroen_Dummy-_und_Effektkodierung_Mehrfaktorielle_Varianzanalyse.pdf`
- 主教材 / Hauptlehrbuch: `Regression.pdf`, Chapter 3.1.3, Chapter 3.3, Chapter 3.4
- 辅助教材 / Hilfslehrbuch: Faraway, `Linear Models with R`, Chapter 13-16
- 扩展补充 / Erweiterung: Faraway, `Extending The Linear Model With R`, Chapter 10 only as preview for random effects

## 覆盖审核 / Abdeckungspruefung

| 来源 / Quelle | 内容点 / Inhaltspunkt | 笔记位置 / Wo |
|---|---|---|
| 讲义 5 | 离散解释变量模型 | 5.1 |
| 讲义 5 | dummy coding / reference coding | 5.2 |
| 讲义 5 | effect coding | 5.3 |
| 讲义 5 | 单因素方差分析 | 5.4 |
| 讲义 5 | 两因素方差分析 | 5.5 |
| 讲义 5 | 主效应与交互 | 5.5, 5.6 |
| 讲义 5 | 离散变量与连续变量结合 / ANCOVA | 5.7 |
| 讲义 5 | 不同类型平方和、前后测比较 | 5.8 |
| 主教材 | 分类协变量的线性模型表示、线性假设和模型比较 | 5.2-5.8 |
| Faraway | ANOVA、ANCOVA、factorial designs、block designs 的 R 实操 | R 练习 |

## 5.1 为什么分类变量也能进线性模型

**遇到的问题 / Problem.**  
线性模型里的解释变量不一定都是连续变量。性别、地区、处理组、班级、政策阶段等都是离散变量。问题是：这些文字类别如何进入 $X$ 矩阵？

**可能的解决路径 / Moegliche Wege.**

- 把类别随便编码成 $1,2,3$：可能错误暗示等距和线性趋势。
- 用 dummy variables：把类别转成若干 $0/1$ 变量。
- 用 effect coding：让参数表达相对总体均值的偏差。

**本课采用的方法及好处 / Methode und Vorteil.**  
把分类变量转成设计矩阵中的列。线性模型仍是：

$$
Y=X\beta+\varepsilon
$$

区别只在于 $X$ 的列现在表示类别对比。好处是：ANOVA、ANCOVA 和回归都能统一到同一个线性模型框架。

**具体解决 / Konkrete Loesung.**  
含 $K$ 个类别的因子，在含截距模型中通常只放 $K-1$ 个 dummy，否则会和截距线性相关。

**R 练习 / R-Uebung.**

```r
data(PlantGrowth)
fit <- lm(weight ~ group, data = PlantGrowth)
model.matrix(fit)[1:10, ]
summary(fit)
```

## 5.2 Dummy coding / Referenzkodierung

**遇到的问题 / Problem.**  
我们需要一个基准类别，其他类别的系数表示相对这个基准的差异。

**本课采用的方法及好处 / Methode und Vorteil.**  
若因子 $C$ 有 $K$ 类，选择第 $1$ 类为参考组：

$$
Y_i=\beta_0+\beta_2D_{i2}+\dots+\beta_KD_{iK}+\varepsilon_i
$$

其中：

$$
D_{ik}=
\begin{cases}
1,& C_i=k\\
0,& C_i\neq k
\end{cases}
$$

解释：

- $\beta_0$ 是参考组均值。
- $\beta_k$ 是第 $k$ 组与参考组的均值差。

**具体解决 / Konkrete Loesung.**  
解释 dummy coding 时必须说清楚参考类别。R 默认通常按 factor levels 的第一个水平作为参考组。

**R 练习 / R-Uebung.**

```r
data(PlantGrowth)
PlantGrowth[["group"]] <- relevel(PlantGrowth[["group"]], ref = "ctrl")
fit <- lm(weight ~ group, data = PlantGrowth)
coef(fit)
tapply(PlantGrowth[["weight"]], PlantGrowth[["group"]], mean)
```

练习 / Uebung: `grouptrt1` 的系数表示什么？  
答案 / Antwort: `trt1` 组均值与 `ctrl` 组均值的差。

## 5.3 Effect coding / Effektkodierung

**遇到的问题 / Problem.**  
有时我们不想把某一组当成唯一参考，而想把参数解释为相对总体均值的偏差。

**本课采用的方法及好处 / Methode und Vorteil.**  
Effect coding 通常让组效应满足约束：

$$
\sum_{k=1}^K \alpha_k=0
$$

模型可写为：

$$
Y_{ik}=\mu+\alpha_k+\varepsilon_{ik}
$$

其中 $\mu$ 是总体均值或平衡设计下的平均组均值，$\alpha_k$ 是第 $k$ 组相对总体均值的偏差。

**具体解决 / Konkrete Loesung.**  
Dummy coding 和 effect coding 表示同一个组均值模型，只是参数含义不同；拟合值和整体检验不会因为编码方式改变。

**R 练习 / R-Uebung.**

```r
data(PlantGrowth)
options(contrasts = c("contr.sum", "contr.poly"))
fit_effect <- lm(weight ~ group, data = PlantGrowth)
coef(fit_effect)
model.matrix(fit_effect)[1:10, ]

options(contrasts = c("contr.treatment", "contr.poly"))
```

## 5.4 单因素方差分析 / Einfache Varianzanalyse

**遇到的问题 / Problem.**  
若只有一个分类解释变量，我们要检验不同组的均值是否相同。

**本课采用的方法及好处 / Methode und Vorteil.**  
组均值模型：

$$
Y_{ij}=\mu+\alpha_j+\varepsilon_{ij}
$$

检验：

$$
H_0:\alpha_1=\dots=\alpha_K=0
$$

等价于所有组均值相同。好处是：ANOVA 可以看作带分类解释变量的线性模型。

**具体解决 / Konkrete Loesung.**  
看 `anova(fit)` 中因子对应的 $F$ 检验。如果显著，只说明至少两个组均值不同；哪两组不同需要进一步对比。

**R 练习 / R-Uebung.**

```r
data(PlantGrowth)
fit <- lm(weight ~ group, data = PlantGrowth)
anova(fit)
boxplot(weight ~ group, data = PlantGrowth)
```

## 5.5 两因素方差分析 / Zweifaktorielle Varianzanalyse

**遇到的问题 / Problem.**  
很多实验同时有两个分类因素，例如处理组和剂量、机器和工人。我们要区分主效应和交互效应。

**本课采用的方法及好处 / Methode und Vorteil.**  
无交互模型：

$$
E(Y\mid A,B)=\mu+\alpha_A+\gamma_B
$$

有交互模型：

$$
E(Y\mid A,B)=\mu+\alpha_A+\gamma_B+(\alpha\gamma)_{AB}
$$

好处是：主效应描述单个因素平均作用，交互描述一个因素的作用是否依赖另一个因素。

**具体解决 / Konkrete Loesung.**  
先画 interaction plot。若线条明显不平行，说明可能存在交互；模型中应加入 `A:B` 或 `A * B`。

**R 练习 / R-Uebung.**

```r
data(ToothGrowth)
ToothGrowth[["dose_f"]] <- factor(ToothGrowth[["dose"]])

fit_add <- lm(len ~ supp + dose_f, data = ToothGrowth)
fit_int <- lm(len ~ supp * dose_f, data = ToothGrowth)

anova(fit_add, fit_int)
interaction.plot(ToothGrowth[["dose_f"]],
                 ToothGrowth[["supp"]],
                 ToothGrowth[["len"]])
```

## 5.6 交互的解释 / Interpretation der Interaktion

**遇到的问题 / Problem.**  
有交互时，主效应不能孤立解释。例如处理 A 的效果可能在不同剂量下不同。

**本课采用的方法及好处 / Methode und Vorteil.**  
若模型为：

$$
Y=\beta_0+\beta_1x+\beta_2z+\beta_3xz+\varepsilon
$$

则 $x$ 的斜率为：

$$
\frac{\partial E(Y\mid x,z)}{\partial x}=\beta_1+\beta_3z
$$

对于两个分类变量，交互表示组均值差异不是简单相加。

**具体解决 / Konkrete Loesung.**  
有交互时优先解释“简单效应”：在某一水平下另一个变量的作用，而不是只读主效应系数。

## 5.7 ANCOVA：离散变量与连续变量结合

**遇到的问题 / Problem.**  
有时既有分组变量，又有连续协变量。只比较组均值会忽略连续协变量的影响。

**本课采用的方法及好处 / Methode und Vorteil.**  
协方差分析模型可写为：

$$
Y_i=\beta_0+\beta_1x_i+\gamma_1D_{i1}+\dots+\gamma_mD_{im}+\varepsilon_i
$$

若允许不同组有不同斜率，加入交互：

$$
Y_i=\beta_0+\beta_1x_i+\gamma D_i+\delta x_iD_i+\varepsilon_i
$$

好处是：可以比较在连续协变量调整后的组差异。

**R 练习 / R-Uebung.**

```r
data(mtcars)
mtcars[["cyl_f"]] <- factor(mtcars[["cyl"]])

fit_parallel <- lm(mpg ~ wt + cyl_f, data = mtcars)
fit_slopes <- lm(mpg ~ wt * cyl_f, data = mtcars)

anova(fit_parallel, fit_slopes)
```

## 5.8 平方和类型与前后测比较

**遇到的问题 / Problem.**  
不平衡设计中，不同类型平方和可能给出不同检验结果；前后测数据也不能简单当成独立两组。

**本课采用的方法及好处 / Methode und Vorteil.**

- Type I: 顺序平方和，依赖变量顺序。
- Type II: 在不含高阶交互的条件下检验主效应。
- Type III: 在包含其他所有项的条件下检验，依赖编码和模型设定。

前后测可使用差值：

$$
D_i=Y_{i,post}-Y_{i,pre}
$$

然后对 $D_i$ 建模；也可用包含个体效应的模型，后续 mixed model 会更系统处理。

**具体解决 / Konkrete Loesung.**  
报告 ANOVA 时必须说明平方和类型、模型公式和编码方式。对前后测，先确认观测是否配对。

**R 练习 / R-Uebung.**

```r
data(sleep)
wide <- reshape(sleep, idvar = "ID", timevar = "group", direction = "wide")
wide[["diff"]] <- wide[["extra.2"]] - wide[["extra.1"]]
t.test(wide[["diff"]])
```

## 5.9 选择题 / Multiple-Choice-Fragen

1. 含截距模型中，$K$ 个类别通常需要几个 dummy？  
   A. $K-1$  
   B. $K$  
   C. $n$  
   D. $1$  
   答案 / Antwort: A

2. Dummy coding 中截距通常表示什么？  
   A. 参考组均值  
   B. 所有残差之和  
   C. 样本量  
   D. 方差  
   答案 / Antwort: A

3. Effect coding 常用约束是什么？  
   A. $\sum_k \alpha_k=0$  
   B. 所有 $\alpha_k=1$  
   C. $Y=0$  
   D. $R^2=0$  
   答案 / Antwort: A

4. 有交互时，主效应应如何解释？  
   A. 结合另一个变量的水平解释  
   B. 完全忽略交互  
   C. 只看截距  
   D. 不能使用线性模型  
   答案 / Antwort: A

5. ANCOVA 主要处理什么结构？  
   A. 连续响应 + 分组变量 + 连续协变量  
   B. 只有一个二元响应  
   C. 只有时间序列误差  
   D. 只有随机森林  
   答案 / Antwort: A

## 5.10 本章总结 / Zusammenfassung

中文总结：

- 离散解释变量通过编码进入设计矩阵。
- Dummy coding 以参考组为基准，系数表示相对参考组差异。
- Effect coding 用约束让组效应围绕总体均值解释。
- ANOVA 是线性模型在分类解释变量场景下的特殊形式。
- 两因素模型要区分主效应和交互。
- ANCOVA 把分类变量和连续协变量放在同一个线性模型里。
- 不平衡设计中平方和类型会影响检验解释。

Deutsche Zusammenfassung:

- Diskrete Einflussgroessen werden durch Kodierung in die Designmatrix aufgenommen.
- Dummykodierung interpretiert Effekte relativ zu einer Referenzkategorie.
- Effektkodierung interpretiert Abweichungen vom Gesamtmittel.
- Varianzanalyse ist ein Spezialfall des linearen Modells.
- Interaktionen bedeuten, dass ein Effekt vom Niveau eines anderen Faktors abhaengt.

## 5.11 专业德语单词汇总 / Fachwortschatz

| Deutsch | 中文 | 说明 |
|---|---|---|
| diskrete Einflussgroesse | 离散解释变量 | 分类变量 |
| Dummykodierung | dummy 编码 | 参考组对比 |
| Referenzkategorie | 参考类别 | 截距对应 |
| Effektkodierung | 效应编码 | 组效应和为零 |
| einfache Varianzanalyse | 单因素方差分析 | one-way ANOVA |
| zweifaktorielle Varianzanalyse | 双因素方差分析 | two-way ANOVA |
| Haupteffekt | 主效应 | main effect |
| Interaktion | 交互 | effect modification |
| Kovarianzanalyse | 协方差分析 | ANCOVA |
| Typ-I-Quadratsumme | I 型平方和 | 顺序平方和 |
| Vorher-Nachher-Vergleich | 前后测比较 | paired structure |
