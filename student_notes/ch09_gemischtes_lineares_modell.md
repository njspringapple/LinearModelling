# 9. Das gemischte lineare Regressionsmodell（混合线性回归模型）

## 章节导航

| 主题 | 对应内容 | 关键领悟 |
|------|---------|---------|
| 1 | 动机：为什么需要混合模型？ | 聚类数据违反独立性假设 |
| 2 | Random-Intercept 模型 | 每组有一个随机截距 γᵢ |
| 3 | 边际模型（Marginal Model） | 将随机效应"积分掉"，归结为一般线性模型 |
| 4 | 纵向数据的扩展 | 个体内还有时间相关结构 |
| 5 | 推断：ML 与 REML | 与第 8 章（一般线性模型）完全对接 |
| 6 | 结果的解释 | ICC、方差分量、固定效应 |
| 7 | 固定效应模型 vs. 混合模型 | 什么时候用哪个？ |
| 8 | 扩展与补充 | 随机斜率、多个随机效应等 |

---

## 1 动机：为什么需要混合模型？

### 1.1 一个具体的例子：阅读促进研究（Leseförderung）

> **研究设计：**
> - **因变量 (Zielgröße):** 阅读能力的改善
> - **自变量 (Einflussgröße):** 是否接受了特殊促进
> - **控制变量 (Störgröße):** 初始水平（Ausgangsniveau）
> - **关键特点:** 实验是以**班级为单位**进行的（klassenweise durchgeführt）

### 1.2 问题出在哪里？

在经典线性模型 $ Y = X\beta + \epsilon $ 中，我们假设：

$$
\epsilon_i \sim N(0, \sigma^2), \quad \text{i.i.d.（独立同分布）}
$$

但如果数据是**按班级聚类**的（Cluster-Daten），同一个班级的学生之间**不独立**！

**直觉理解：**
> 想象两个班级：A 班老师很严格，B 班老师很松。
> 即使没有特殊促进，A 班学生的表现可能系统性地高于 B 班。
> 这种**班级层面的系统差异**使得同班学生的残差是相关的。

### 1.3 解决思路的演进

| 步骤 | 方法 | 问题 |
|------|------|------|
| 第一步 | 忽略聚类结构 | ❌ 违反独立性，标准误偏小，p 值偏小 |
| 第二步 | 引入班级固定效应（Klasseneffekt als fester Effekt） | ❌ 参数太多（每个班级一个参数） |
| 第三步 | 引入班级随机效应（Klasseneffekt als zufälliger Effekt） | ✅ 只增加一个方差参数 σ²_γ |

**关键洞察：**
> 固定效应 = 每个班级估计一个截距（g 个参数）
> 随机效应 = 假设班级截距来自一个正态分布 N(0, σ²_γ)，只需估计一个参数 σ²_γ
> 这是一种**参数节省**的策略，也称为**部分池化（partial pooling）**

---

## 2 带简单随机效应的模型（Random-Intercept-Modell）

### 2.1 模型定义

也称为**方差分量模型（Varianzkomponenten-Modell）**。

考虑分组数据，组别索引为 i：

$$
Y_{ij} = x'_{ij}\beta + \gamma_i + \epsilon_{ij}, \quad i = 1, \ldots, g; \quad j = 1, \ldots, n_i
$$

其中：
- $i$：组别（如班级、医院、个体）
- $j$：组内观测（如学生、患者、时间点）
- $g$：组的总数
- $n_i$：第 i 组的观测数
- $x'_{ij}\beta$：**固定效应部分**（与经典线性模型相同）
- $\gamma_i$：第 i 组的**随机截距**（Random Intercept）
- $\epsilon_{ij}$：个体层面的**残差**

### 2.2 分布假设

$$
\epsilon_{ij} \sim N(0, \sigma^2) \quad \text{i.i.d.}
$$

$$
\gamma_i \sim N(0, \sigma^2_\gamma) \quad \text{i.i.d.}
$$

$$
\gamma_i \text{ 与 } \epsilon \text{ 相互独立}
$$

**三个独立性假设要记清：**
> 1. 残差之间 i.i.d.
> 2. 随机效应之间 i.i.d.
> 3. 残差与随机效应之间独立

### 2.3 模型的直觉图像

想象每个组有自己的"基线水平"，但这个基线水平不是完全自由的参数，
而是从一个正态分布中"抽取"的：

```
组 1: Y₁ⱼ = x'₁ⱼβ + γ₁ + ε₁ⱼ    （γ₁ 偏高 → 这组整体偏高）
组 2: Y₂ⱼ = x'₂ⱼβ + γ₂ + ε₂ⱼ    （γ₂ 偏低 → 这组整体偏低）
组 3: Y₃ⱼ = x'₃ⱼβ + γ₃ + ε₃ⱼ    （γ₃ ≈ 0 → 这组接近总体平均）
...
所有 γᵢ ~ N(0, σ²_γ)
```

> 如果 σ²_γ 很大 → 组间差异大
> 如果 σ²_γ ≈ 0 → 各组几乎没有差异，退化为经典线性模型

### 2.4 参数清单

| 参数 | 类型 | 说明 |
|------|------|------|
| β（p 维） | 固定效应 | 与经典线性模型相同 |
| σ² | 方差参数 | 个体层面残差方差 |
| σ²_γ | 方差参数 | 随机效应的方差（方差分量） |

**总共需要估计的参数：** p + 2（而不是 p + g）

这就是随机效应的核心优势——**参数节省**。

---

## 3 边际模型（Das marginale Modell）

### 3.1 从条件模型到边际模型

原始模型是"条件于 γᵢ"的：

$$
Y_{ij} | \gamma_i = x'_{ij}\beta + \gamma_i + \epsilon_{ij}
$$

我们可以将 γᵢ 和 εᵢⱼ 合并为一个新的误差项 δᵢⱼ：

$$
Y_{ij} = x'_{ij}\beta + \delta_{ij}
$$

其中：

$$
\delta_{ij} = \epsilon_{ij} + \gamma_i
$$

### 3.2 边际误差项的性质

**方差：**

$$
\text{Var}(\delta_{ij}) = \text{Var}(\epsilon_{ij}) + \text{Var}(\gamma_i) = \sigma^2 + \sigma^2_\gamma
$$

> 因为 εᵢⱼ 和 γᵢ 独立，方差直接相加。

**同组内的协方差（j₁ ≠ j₂，同一组 i）：**

$$
\text{Cov}(\delta_{ij_1}, \delta_{ij_2}) = \text{Cov}(\epsilon_{ij_1} + \gamma_i, \; \epsilon_{ij_2} + \gamma_i)
$$

展开：

$$
= \text{Cov}(\epsilon_{ij_1}, \epsilon_{ij_2}) + \text{Cov}(\epsilon_{ij_1}, \gamma_i) + \text{Cov}(\gamma_i, \epsilon_{ij_2}) + \text{Cov}(\gamma_i, \gamma_i)
$$

$$
= 0 + 0 + 0 + \sigma^2_\gamma = \sigma^2_\gamma
$$

> **同组内的观测有正相关！** 相关性来源就是共享的随机效应 γᵢ。

**不同组之间的协方差（i₁ ≠ i₂）：**

$$
\text{Cov}(\delta_{i_1 j_1}, \delta_{i_2 j_2}) = 0
$$

> 不同组之间仍然是独立的。

### 3.3 协方差结构的矩阵形式

考虑第 i 组的 nᵢ 个观测，它们的协方差矩阵为：

$$
\text{Var}(\delta_i) = \sigma^2 I_{n_i} + \sigma^2_\gamma e_i e'_i
$$

其中 $ e_i $ 是长度为 $ n_i $ 的全 1 向量。

**具体展开（以 nᵢ = 3 为例）：**

$$
\text{Var}(\delta_i) = \sigma^2 \begin{pmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & 1 \end{pmatrix} + \sigma^2_\gamma \begin{pmatrix} 1 \\ 1 \\ 1 \end{pmatrix} \begin{pmatrix} 1 & 1 & 1 \end{pmatrix}
$$

$$
= \begin{pmatrix} \sigma^2 + \sigma^2_\gamma & \sigma^2_\gamma & \sigma^2_\gamma \\ \sigma^2_\gamma & \sigma^2 + \sigma^2_\gamma & \sigma^2_\gamma \\ \sigma^2_\gamma & \sigma^2_\gamma & \sigma^2 + \sigma^2_\gamma \end{pmatrix}
$$

> **对角线：** σ² + σ²_γ（总方差）
> **非对角线：** σ²_γ（组内协方差 = 组内相关性的来源）

这就是所谓的**复合对称结构（Compound Symmetry）**。

### 3.4 整体向量形式

$$
Y = X\beta + \delta
$$

$$
\delta \sim N\Big(0, \; \sigma^2 I + \text{diag}[\sigma^2_\gamma e_i e'_i]\Big)
$$

这里 $\text{diag}[\sigma^2_\gamma e_i e'_i]$ 表示将各组的 $\sigma^2_\gamma e_i e'_i$ 块对角拼接。

**整体协方差矩阵 V 的结构：**

$$
V = \begin{pmatrix} \sigma^2 I_{n_1} + \sigma^2_\gamma e_1 e'_1 & 0 & \cdots & 0 \\ 0 & \sigma^2 I_{n_2} + \sigma^2_\gamma e_2 e'_2 & \cdots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \cdots & \sigma^2 I_{n_g} + \sigma^2_\gamma e_g e'_g \end{pmatrix}
$$

> **块对角矩阵！** 组间独立，组内有 Compound Symmetry 结构。

**核心认识：** 边际模型就是一个**一般线性模型（Allgemeines lineares Modell）**
$ Y = X\beta + \delta, \; \delta \sim N(0, V) $，其中 V 不是 σ²I，而是有特殊结构的矩阵。

→ 因此，第 8 章关于一般线性模型的所有推断方法都可以直接使用！

---

## 4 纵向数据（Longitudinale Daten / Wiederholungsmessungen）

### 4.1 为什么需要进一步扩展？

在聚类数据（如班级中的学生）中，同一组内的残差是 i.i.d. 的（给定 γᵢ）。

但在**纵向数据**（同一个人在不同时间点的重复测量）中，
即使给定 γᵢ，**相邻时间点的残差可能仍然相关**！

例如：今天的血压和明天的血压比今天和一个月后的血压更相似。

### 4.2 扩展模型

$$
Y_{ij} = x'_{ij}\beta + \gamma_i + \epsilon_{ij}, \quad i = 1, \ldots, g; \quad j = 1, \ldots, n_i
$$

$$
\epsilon_i \sim N(0, \sigma^2 \Sigma_i) \quad \text{（注意：不再是 } \sigma^2 I \text{！）}
$$

$$
\gamma_i \sim N(0, \sigma^2_\gamma) \quad \text{i.i.d.}
$$

$$
\gamma_i \text{ 与 } \epsilon_i \text{ 独立}
$$

关键变化：
- $ \epsilon_i $ 是第 i 个人的**残差向量**（不再是标量）
- $ \Sigma_i $ 是 $ n_i \times n_i $ 的**相关矩阵**，描述了个体内的时间相关性

### 4.3 边际模型

$$
Y = X\beta + \delta
$$

$$
\delta \sim N\Big(0, \; \text{diag}(\sigma^2 \Sigma_i) + \text{diag}[\sigma^2_\gamma e_i e'_i]\Big)
$$

> **两层相关性叠加：**
> 1. 随机效应 γᵢ 带来的组内等相关（Compound Symmetry 部分）
> 2. 个体内残差的时间相关性（Σᵢ 部分）

### 4.4 常见的相关结构 Σᵢ

| 结构名称                    | 描述                 | 参数数量     |
| ----------------------- | ------------------ | -------- |
| 独立（Identity）            | Σ = I              | 0        |
| 复合对称（Compound Symmetry） | 所有 off-diagonal 相等 | 1        |
| AR(1)                   | 相邻相关 ρ，距离 k 的相关 ρᵏ | 1        |
| 无结构（Unstructured）       | 所有元素自由估计           | n(n-1)/2 |

---

## 5 推断（Inferenz）

### 5.1 与一般线性模型的对接

边际模型 $Y = X\beta + \delta, \; \delta \sim N(0, V(\vartheta))$ 就是第 8 章的一般线性模型。

其中 $\vartheta$ 是所有描述方差结构的参数，例如 $\vartheta = (\sigma^2, \sigma^2_\gamma)$。

### 5.2 Log-Likelihood

$$
\ell(\beta, \vartheta) = -\frac{1}{2}\Big(\ln|V(\vartheta)| + (Y - X\beta)' V^{-1}(\vartheta)(Y - X\beta)\Big)
$$

（忽略常数项）

### 5.3 估计方法

**ML 估计：** 同时关于 β 和 ϑ 最大化 ℓ(β, ϑ)

**REML 估计：**
> - 先消除 β 的影响
> - 只基于残差信息估计 ϑ
> - 然后用估计的 ϑ̂ 回代得到 β̂

回顾第 8 章：

$$
\hat{\beta}(\vartheta) = (X' V^{-1}(\vartheta) X)^{-1} X' V^{-1}(\vartheta) Y
$$

先把 β 用上式替换掉，得到只关于 ϑ 的 profile likelihood（或 restricted likelihood），
然后最大化得到 ϑ̂_REML，最后代回得到 β̂。

**REML 的优势：**
> REML 考虑了估计 β 时损失的自由度，因此对方差参数的估计是**无偏的**（或偏差更小）。
> 类似于用 n-p 而不是 n 作为 σ² 估计量的分母。

### 5.4 实际求解

由于 V(ϑ) 依赖于 ϑ，log-likelihood 通常没有解析解，
需要**数值优化**（如 Newton-Raphson、Fisher Scoring、EM 算法）。

在 R 中，`lme4::lmer()` 和 `nlme::lme()` 都实现了这些算法。

---

## 6 结果的解释（Interpretation der Ergebnisse）

### 6.1 固定效应 βₖ

与经典线性模型完全一样解释：
- 在其他变量不变的情况下，xₖ 增加一个单位，Y 的期望变化 βₖ
- 可以包含交互作用等所有经典变体

### 6.2 残差方差 σ²_ε

个体层面（组内）的残差方差，即"纯随机噪声"的大小。

### 6.3 随机效应的方差 σ²_γ（方差分量）

描述组间差异的大小。

> σ²_γ 大 → 各组的基线水平差异很大
> σ²_γ ≈ 0 → 各组几乎没有差异

### 6.4 ⭐ Intra-Class Correlation (ICC)

$$
\text{ICC} = \frac{\sigma^2_\gamma}{\sigma^2_\gamma + \sigma^2_\epsilon}
$$

**含义：**
> ICC 是**同组内两个观测之间的相关系数**。

**推导：**

$$
\text{Cor}(Y_{ij_1}, Y_{ij_2}) = \frac{\text{Cov}(\delta_{ij_1}, \delta_{ij_2})}{\sqrt{\text{Var}(\delta_{ij_1}) \cdot \text{Var}(\delta_{ij_2})}} = \frac{\sigma^2_\gamma}{(\sigma^2_\gamma + \sigma^2_\epsilon)} = \text{ICC}
$$

**直觉解读：**
> ICC = 组间方差占总方差的比例
>
> - ICC ≈ 0：几乎没有组效应，经典线性模型就够了
> - ICC ≈ 1：几乎所有变异都来自组间差异
> - ICC = 0.2：20% 的变异由组效应解释，是很常见的经验值

**考试高频考点！** 经常要求解释 ICC 的含义并计算。

### 6.5 相关参数

如果模型中包含个体内的相关结构 Σᵢ，还需要解释相关参数（如 AR(1) 中的 ρ）。

---

## 7 固定效应模型 vs. 混合模型

### 7.1 争议背景

在纵向数据（Panel-Daten）中，经济学家和统计学家对于
**是用固定效应模型还是混合模型**一直存在争论。

### 7.2 支持随机效应模型的论据

**（1）参数数量大大减少**
> 固定效应模型：每个个体一个截距（g 个参数）
> 随机效应模型：只需一个方差参数 σ²_γ

**（2）可以纳入时间恒定的个体特征**
> 固定效应模型中，个体截距"吸收"了所有个体层面变量的效应
>（如性别、种族），这些变量的系数无法估计。
> 随机效应模型中可以！

**（3）随机效应有合理的可解释性**
> "个体之间存在未观测的异质性，且这种异质性服从正态分布"

### 7.3 反对随机效应模型的论据

**（1）关键假设可能被违反**
> 随机效应模型假设 γᵢ 与 xᵢⱼ **独立**。
> 如果不独立（例如：不可观测的个体特征与自变量相关），
> 则固定效应的估计会有**偏差**（omitted variable bias）！
> 这在经济学中被称为 **"内生性问题"**。

**（2）只关心组内效应**
> 如果研究者只对"一个人随时间变化的效应"感兴趣（within effect），
> 固定效应模型更合适，因为它完全消除了个体间差异。

**（3）模型更简单、假设更少**
> 固定效应模型不需要对 γᵢ 做分布假设。

### 7.4 实际建议

> **Hausman 检验**可以检验随机效应模型的关键假设是否成立。
> 如果 Hausman 检验拒绝了原假设（γᵢ 与 x 独立），则应使用固定效应模型。
> 如果不拒绝，则随机效应模型更高效。

---

## 8 扩展与补充

### 8.1 多个随机效应

例如在语言学实验中：

$$
Y_{ijk} = x'_{ijk}\beta + \gamma_i^{(\text{Person})} + \gamma_j^{(\text{Wort})} + \epsilon_{ijk}
$$

- 受试者（Person）和词语（Wort）都是随机效应
- 这被称为**交叉随机效应（crossed random effects）**

### 8.2 随机斜率（Random Slopes）

不仅截距因组而异，**斜率也因组而异**：

$$
Y_{ij} = (\beta_0 + \gamma_{0i}) + (\beta_1 + \gamma_{1i}) x_{ij} + \epsilon_{ij}
$$

其中：

$$
\begin{pmatrix} \gamma_{0i} \\ \gamma_{1i} \end{pmatrix} \sim N\left(\begin{pmatrix} 0 \\ 0 \end{pmatrix}, \begin{pmatrix} \sigma^2_{\gamma_0} & \sigma_{\gamma_{01}} \\ \sigma_{\gamma_{01}} & \sigma^2_{\gamma_1} \end{pmatrix}\right)
$$

> 截距和斜率的随机效应可以是相关的（通过 σ_γ₀₁）。

### 8.3 随机效应的预测（BLUP）

虽然 γᵢ 不是参数（而是随机变量），但我们可以"预测"它们。

**Best Linear Unbiased Predictor (BLUP)：**

$$
\hat{\gamma}_i = \sigma^2_\gamma e'_i V_i^{-1}(Y_i - X_i \hat{\beta})
$$

> BLUP 是一种**收缩估计量（Shrinkage Estimator）**：
> 它将每组的"样本偏差"向总体均值 0 拉回，
> 样本量小的组收缩更多，样本量大的组收缩更少。

### 8.4 不同类型的残差

| 残差类型 | 定义                                             | 用途            |
| ---- | ---------------------------------------------- | ------------- |
| 边际残差 | $Y_{ij} - x'_{ij}\hat{\beta}$                  | 包含随机效应的部分     |
| 条件残差 | $Y_{ij} - x'_{ij}\hat{\beta} - \hat{\gamma}_i$ | 排除随机效应后的"纯残差" |

### 8.5 不同类型的 R²

在混合模型中，R² 的定义不唯一：
- **边际 R²（marginal R²）：** 只考虑固定效应解释了多少方差
- **条件 R²（conditional R²）：** 固定效应 + 随机效应一起解释了多少方差

> Nakagawa & Schielzeth (2013) 提出了广泛使用的定义。

---

## 9 总结与考试要点

### 9.1 核心公式速查

| 内容             | 公式                                                           |     |                                 |
| -------------- | ------------------------------------------------------------ | --- | ------------------------------- |
| 条件模型           | $Y_{ij} = x'_{ij}\beta + \gamma_i + \epsilon_{ij}$           |     |                                 |
| 边际模型           | $Y_{ij} = x'_{ij}\beta + \delta_{ij}$                        |     |                                 |
| 边际方差           | $\text{Var}(\delta_{ij}) = \sigma^2 + \sigma^2_\gamma$       |     |                                 |
| 组内协方差          | $\text{Cov}(\delta_{ij_1}, \delta_{ij_2}) = \sigma^2_\gamma$ |     |                                 |
| ICC            | $\frac{\sigma^2_\gamma}{\sigma^2_\gamma + \sigma^2}$         |     |                                 |
| Log-Likelihood | $\ell(\beta, \vartheta) = -\frac{1}{2}(\ln                   | V   | + (Y-X\beta)'V^{-1}(Y-X\beta))$ |

### 9.2 概念辨析

| 概念 | 固定效应 | 随机效应 |
|------|---------|---------|
| 本质 | 未知常数 | 随机变量 |
| 估计 | 直接估计每一个 | 估计其分布的参数（σ²_γ） |
| 参数数量 | 每组一个 | 全局一个方差参数 |
| 对 x 的假设 | 可与 x 相关 | 假设与 x 独立 |
| 可以估计组间变量效应？ | 否 | 是 |

### 9.3 典型考试题型

**题型 1：** 给定场景，判断应该用随机效应还是固定效应。
> 提示：看分组变量是"随机抽样的"还是"研究者刻意选定的"。

**题型 2：** 写出模型方程，推导边际方差和协方差。
> 关键：利用 Var(A+B) = Var(A) + Var(B) + 2Cov(A,B)，以及独立性假设。

**题型 3：** 计算并解释 ICC。
> 公式 + 含义（组间方差占比 = 组内相关系数）。

**题型 4：** 比较 ML 和 REML 估计。
> REML 对方差参数估计更优（类似 n-p 校正）。

**题型 5：** 解释 R 输出中的 random effects 和 fixed effects 部分。

---

## 10 与前面章节的联系

```
第 1-6 章：经典线性模型 Y = Xβ + ε, ε ~ N(0, σ²I)
    ↓ 放宽 Var(ε) = σ²I
第 8 章：一般线性模型 Y = Xβ + ε, ε ~ N(0, V)
    ↓ V 具有特殊的"随机效应"结构
第 9 章：混合线性模型 Y = Xβ + Zγ + ε
    ↓ γ ~ N(0, G), ε ~ N(0, R)
    ↓ 边际化后：Y ~ N(Xβ, ZGZ' + R) → 回到第 8 章
```

> **混合模型是一般线性模型的特例**，
> 其中协方差矩阵 V = ZGZ' + R 具有来自随机效应的特殊结构。
> 所有第 8 章的推断工具（GLS、ML、REML）直接适用。