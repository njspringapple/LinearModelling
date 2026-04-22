# 7. Modelldiagnose（模型诊断）

**Lineare Modelle SoSe 2025 — Sabine Hoffmann, LMU**

理论推导 + 直觉理解 + 每一步都搞懂

---

## 章节导航

| 主题 | 对应内容 | 关键领悟 |
|------|---------|---------|
| 1 | 残差的类型 | 普通残差 → 标准化 → 学生化 → 交叉验证 → 递归 |
| 2 | 正态性违反 | KQ 估计仍无偏，但置信区间和预测区间失效 |
| 3 | 异方差 | 估计无偏但不再最优，检验和 CI 不正确 |
| 4 | 误差相关性 | Durbin-Watson 检验 |
| 5 | 异常值与高杠杆点 | Leverage、Cook 距离 |
| 6 | 模型设定错误 | Partial Leverage Plot |
| 7 | 多重共线性 | 条件数、VIF |
| 8 | 测量误差 | 向零偏误（attenuation bias） |

---

## 0 核心出发点

模型为：

$$
\mathbf{Y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon}, \quad \boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{I})
$$

模型诊断的目的是检查这些假设是否成立。由于假设针对的是**误差项 $\varepsilon$**（而非 $Y$ 的边际分布！），所以我们主要通过**残差**来诊断。

> ⚠️ 重要提醒：模型假设**不是关于 $Y$ 的边际分布**的，而是关于**给定 $\mathbf{X}$ 后 $\varepsilon$ 的条件分布**的。所以不能简单看 $Y$ 是否正态。

---

## 1 各类残差（Verschiedene Typen von Residuen）

### 1.1 普通残差（Gewöhnliche Residuen）

$$
\hat{\varepsilon} = \mathbf{Q}\boldsymbol{\varepsilon} = (\mathbf{I} - \mathbf{P})\boldsymbol{\varepsilon}
$$

其中 $\mathbf{P} = \mathbf{X}(\mathbf{X}'\mathbf{X})^{-1}\mathbf{X}'$ 是帽子矩阵（Hat-Matrix），$\mathbf{Q} = \mathbf{I} - \mathbf{P}$ 是残差矩阵。

**关键事实**：残差的方差**不相等**！

$$
\text{Var}(\hat{\varepsilon}_i) = q_{ii} \cdot \sigma^2
$$

其中 $q_{ii} = 1 - h_{ii}$ 是 $\mathbf{Q}$ 的第 $i$ 个对角元素。

这意味着：即使真正的误差项 $\varepsilon_i$ 全部等方差（$\sigma^2$），残差 $\hat{\varepsilon}_i$ 的方差也**随观测点不同而不同**。原因在于帽子矩阵 $\mathbf{P}$ 的作用——离数据中心越远的点，$h_{ii}$ 越大，对应的残差方差 $q_{ii}\sigma^2$ 就越小。

**例子**：课件给了一个 $5 \times 2$ 的设计矩阵，$\mathbf{Q}$ 的对角元素分别为 $2/5, 7/10, 4/5, 7/10, 2/5$——中间的点方差最大，两端的点方差最小。

---

### 1.2 标准化残差（Standardisierte Residuen）

为了让残差可比较，除以它自己的标准差估计：

$$
r_i := \frac{\hat{\varepsilon}_i}{\hat{\sigma}\sqrt{q_{ii}}}, \quad i = 1, \ldots, n
$$

其中 $q_{ii} = 1 - h_{ii}$。

**关于 $h_{ii}$ 和 $q_{ii}$ 的性质**：

$$
\frac{1}{n} \leq h_{ii} \leq 1, \quad \sum h_{ii} = p'
$$

$$
0 \leq q_{ii} \leq 1 - \frac{1}{n}, \quad \sum q_{ii} = n - p'
$$

如果模型正确，标准化残差 $r_i$ 应该近似服从 $N(0,1)$。

---

### 1.3 学生化残差（Studentisierte Residuen）

**问题**：标准化残差中，$\hat{\sigma}$ 的计算用到了**所有残差**（包括 $\hat{\varepsilon}_i$ 本身）。如果第 $i$ 个观测是异常值，它的大残差会**膨胀** $\hat{\sigma}$，从而**掩盖**自己的异常性。

**解决方案**：用**删除第 $i$ 个观测后**的 $\hat{\sigma}_{(i)}$：

$$
r_i^* := \frac{\hat{\varepsilon}_i}{\hat{\sigma}_{(i)}\sqrt{q_{ii}}}, \quad i = 1, \ldots, n
$$

这就是学生化残差（也叫外部学生化残差）。它服从 $t(n - p' - 1)$ 分布，因此可以直接用 $t$ 分布的临界值来判断异常值。

---

### 1.4 交叉验证残差（Kreuzvalidierungs-Residuen）

**定义**：用**删除第 $i$ 个观测后**的估计 $\hat{\boldsymbol{\beta}}_{(i)}$ 来预测第 $i$ 个观测：

$$
e_{(i)} := y_i - \mathbf{x}_i'\hat{\boldsymbol{\beta}}_{(i)}
$$

**PRESS 统计量**（Predicted Residual Sum of Squares）：

$$
\text{PRESS} := \sum_{i=1}^n e_{(i)}^2
$$

这是一种 leave-one-out 交叉验证的度量，可用于模型比较。

**核心公式**（极其重要）：

$$
e_{(i)} = \frac{\hat{\varepsilon}_i}{q_{ii}} = \frac{\hat{\varepsilon}_i}{1 - h_{ii}}
$$

这个公式意味着：**不需要真的跑 $n$ 次回归**！只需拟合一次完整模型，就能算出所有交叉验证残差。

**方差**：

$$
\text{Var}(e_{(i)}) = \frac{1}{q_{ii}} \sigma^2 > \sigma^2
$$

交叉验证残差的方差总是大于 $\sigma^2$，这是合理的——用少一个观测来预测，预测误差自然更大。

---

### 1.5 公式 $e_{(i)} = \hat{\varepsilon}_i / q_{ii}$ 的证明

这个证明分几步，核心工具是 Sherman-Morrison-Woodbury 公式（课件中叫 Lemma）。

**Sherman-Morrison 公式**：对于可逆矩阵 $A$、向量 $b$、标量 $\lambda$：

$$
(A + \lambda bb')^{-1} = A^{-1} - \frac{\lambda}{1 + \lambda b'A^{-1}b} A^{-1}bb'A^{-1}
$$

**步骤 1**：删除第 $i$ 个观测相当于从 $\mathbf{X}'\mathbf{X}$ 中减去 $\mathbf{x}_i\mathbf{x}_i'$：

$$
(\mathbf{X}_{(i)}'\mathbf{X}_{(i)})^{-1} = (\mathbf{X}'\mathbf{X} - \mathbf{x}_i\mathbf{x}_i')^{-1}
$$

用 Sherman-Morrison 公式（取 $\lambda = -1$）：

$$
(\mathbf{X}_{(i)}'\mathbf{X}_{(i)})^{-1}\mathbf{x}_i = \frac{1}{1 - h_{ii}}(\mathbf{X}'\mathbf{X})^{-1}\mathbf{x}_i
$$

**步骤 2**：从正规方程出发，推导出：

$$
\hat{\boldsymbol{\beta}} - \hat{\boldsymbol{\beta}}_{(i)} = \frac{\hat{\varepsilon}_i}{1 - h_{ii}} (\mathbf{X}'\mathbf{X})^{-1}\mathbf{x}_i
$$

**直觉**：删除一个观测后参数估计的变化，与该观测的残差成正比，与其杠杆值成正比。

**步骤 3**：

$$
e_{(i)} = y_i - \mathbf{x}_i'\hat{\boldsymbol{\beta}}_{(i)} = y_i - \mathbf{x}_i'\hat{\boldsymbol{\beta}} + \mathbf{x}_i'(\hat{\boldsymbol{\beta}} - \hat{\boldsymbol{\beta}}_{(i)}) = \hat{\varepsilon}_i + \frac{h_{ii}}{1 - h_{ii}}\hat{\varepsilon}_i = \frac{\hat{\varepsilon}_i}{1 - h_{ii}}
$$

证毕。■

---

### 1.6 递归残差（Rekursive Residuen）

用于**时间序列**数据。用前 $i-1$ 个观测估计参数，预测第 $i$ 个观测：

$$
\omega_i := \frac{y_i - \mathbf{x}_i'\hat{\boldsymbol{\beta}}_{[i-1]}}{\sqrt{1 - \mathbf{x}_i'(\mathbf{X}_{[i-1]}'\mathbf{X}_{[i-1]})^{-1}\mathbf{x}_i}}, \quad i = p'+1, \ldots, n
$$

如果模型正确且参数稳定，递归残差应该是独立同分布的 $N(0, \sigma^2)$。如果参数发生了结构性变化（structural break），递归残差会在变化点附近偏离零。

---

### 1.7 各类残差对比总结

| 残差类型                         | 公式                                                          | 方差                     | 用途         |
| ---------------------------- | ----------------------------------------------------------- | ---------------------- | ---------- |
| 普通残差 $\hat{\varepsilon}_i$ | $y_i - \mathbf{x}_i'\hat{\boldsymbol{\beta}}$             | $q_{ii}\sigma^2$（不等） | 基本诊断图      |
| 标准化残差 $r_i$                | $\hat{\varepsilon}_i / (\hat{\sigma}\sqrt{q_{ii}})$       | ≈1                     | 可比较的残差图    |
| 学生化残差 $r_i^*$              | $\hat{\varepsilon}_i / (\hat{\sigma}_{(i)}\sqrt{q_{ii}})$ | 服从 $t(n-p'-1)$         | 异常值检测      |
| 交叉验证残差 $e_{(i)}$           | $\hat{\varepsilon}_i / q_{ii}$                            | $\sigma^2/q_{ii}$    | 模型预测能力评估   |
| 递归残差 $\omega_i$            | 见公式 (6.10)                                                  | $\sigma^2$（如果正确）     | 时间序列结构变化检测 |

---

## 2 误差项非正态（Die Störterme sind nicht normalverteilt）

**原因**：$Y$ 是计数数据、生存时间、比例等本质上非正态的量。

**后果**（分轻重）：

| 影响对象 | 严重程度 | 说明 |
|---------|---------|------|
| $\hat{\boldsymbol{\beta}}$ 无偏性 | ✅ 不受影响 | Gauss-Markov 定理不需要正态假设 |
| $\hat{\boldsymbol{\beta}}$ BLUE 性 | ✅ 不受影响 | 同上 |
| F-检验 | ⚠️ 通常稳健 | 大样本下几乎不受影响 |
| 置信区间 | ❌ 小样本下可能不准 | 依赖正态假设 |
| **预测区间** | ❌❌ **尤其受影响** | 直接依赖误差分布的形状 |

**诊断方法**：残差的偏度和峰度、Normal Q-Q Plot

**治疗方法**：对 $Y$ 做变换（如取对数）；使用广义线性模型（GLM）

---

## 3 异方差（Heterogene Varianzen）

**定义**：$\text{Var}(\varepsilon_i)$ 依赖于 $i$（即不同观测的误差方差不同）。

**原因**：计数数据（方差随均值增大）、比例数据、分组数据、乘性误差结构（$\sigma_i$ 与 $Y_i$ 的大小相关）。

**后果**：

| 影响对象 | 状态 |
|---------|------|
| $\hat{\boldsymbol{\beta}}$ 无偏性 | ✅ 仍然无偏 |
| $\hat{\boldsymbol{\beta}}$ 效率 | ❌ 不再是 BLUE（方差不是最小的） |
| 置信区间和检验 | ❌ **不再正确** |

**诊断方法**：残差图（$\hat{\varepsilon}_i$ 对 $\hat{Y}_i$）。典型的异方差模式是"喇叭形"——随着 $\hat{Y}_i$ 增大，残差的散布也增大。

**治疗方法**：对 $Y$ 做变换（如 $\ln(Y)$ 或 $\sqrt{Y}$）；加权最小二乘（WLS）。

---

## 4 误差项相关（Korrelation zwischen den Störtermen）

**定义**：$\text{Cov}(\varepsilon_i, \varepsilon_j) \neq 0$，对某些 $i \neq j$。

**原因**：时间序列结构（相邻时刻的误差正相关）、空间结构、未建模的分组效应。

**后果**：

| 影响对象 | 状态 |
|---------|------|
| $\hat{\boldsymbol{\beta}}$ 无偏性 | ✅ 仍然无偏 |
| $\hat{\boldsymbol{\beta}}$ 效率 | ❌ 不再最优 |
| $\hat{\sigma}^2$ | ❌ **有偏** |
| 置信区间和 F-检验 | ❌❌ **严重失效**（通常低估标准误） |

> 正相关的误差 → $\hat{\sigma}^2$ 被**低估** → 标准误被低估 → p 值过小 → **虚假显著**！

**诊断方法**：残差对时间的图、$\hat{\varepsilon}_i$ 对 $\hat{\varepsilon}_{i-1}$ 的散点图、Durbin-Watson 检验。

**治疗方法**：时间序列方法、加入趋势和季节项、加权最小二乘。

---

### 4.1 Durbin-Watson 检验

检验相邻误差项的一阶自相关 $\rho$：

$$
d := \frac{\sum_{i=2}^n (\hat{\varepsilon}_i - \hat{\varepsilon}_{i-1})^2}{\sum_{i=1}^n \hat{\varepsilon}_i^2} \approx 2(1 - \hat{\rho})
$$

**$H_0: \rho = 0$**（无自相关）

| $d$ 的值 | 含义 |
|---------|------|
| $d \approx 2$ | 无自相关 ✅ |
| $d$ 接近 0 | 强正自相关（$\rho > 0$） |
| $d$ 接近 4 | 强负自相关（$\rho < 0$） |

**判定规则**（$d_l, d_u$ 是查表得到的临界值）：

| $d$ 的范围 | 决策 |
|----------|------|
| $[0, d_l]$ 或 $[4-d_l, 4]$ | 拒绝 $H_0$ |
| $[d_u, 4-d_u]$ | 不拒绝 $H_0$ |
| $(d_l, d_u)$ 或 $(4-d_u, 4-d_l)$ | 不确定区域 |

---

## 5 异常值与高杠杆点（Ausreißer und High Leverage Points）

### 5.1 两个概念的区别

| 概念 | 定义 | 空间 |
|------|------|------|
| **高杠杆点（High leverage point）** | $\mathbf{x}_i$ 远离数据中心 | **$X$ 空间**中的异常 |
| **异常值（Ausreißer）** | $\varepsilon_i$ 的绝对值很大 | **$Y$ 方向**的异常 |

最危险的是：**同时是异常值又是高杠杆点**的观测——它既有能力拉动回归线（高杠杆），又偏离了正确的模式（大残差）。

---

### 5.2 Leverage（杠杆值）

$$
h_{ii} := \mathbf{x}_i'(\mathbf{X}'\mathbf{X})^{-1}\mathbf{x}_i
$$

这是帽子矩阵 $\mathbf{P}$ 的第 $i$ 个对角元素。

**性质**：

$$
\frac{1}{n} \leq h_{ii} \leq 1
$$

**判定标准**：

| $h_{ii}$ 的值 | 判断 |
|-------------|------|
| $h_{ii} = p'/n$ | 平均水平 |
| $h_{ii} > 2p'/n$ | **高杠杆点** ⚠️ |

**直觉理解**：$h_{ii}$ 衡量的是第 $i$ 个观测点在 $X$ 空间中距离数据中心的 Mahalanobis 距离。$h_{ii}$ 越大，这个点对回归线的"拉力"越大。

另一个理解：$\hat{Y}_i = \sum_j h_{ij} Y_j$，所以 $h_{ii}$ 是第 $i$ 个拟合值中来自 $Y_i$ 自身的权重。如果 $h_{ii} = 1$，拟合值完全由自己决定——回归线被迫穿过这个点。

---

### 5.3 Cook 距离

$$
D_i = \frac{(\hat{\boldsymbol{\beta}}_{(i)} - \hat{\boldsymbol{\beta}})'(\mathbf{X}'\mathbf{X})(\hat{\boldsymbol{\beta}}_{(i)} - \hat{\boldsymbol{\beta}})}{\hat{\sigma}^2 p'}
$$

这衡量的是：**删除第 $i$ 个观测后，所有拟合值整体变化了多少**。

**等价分解**（极其优美）：

$$
D_i = \underbrace{\frac{r_i^2}{p'}}_{\text{残差大小}} \cdot \underbrace{\frac{h_{ii}}{1 - h_{ii}}}_{\text{杠杆大小}}
$$

Cook 距离 = **残差分量** × **杠杆分量**。

也就是说，一个观测对回归结果影响大，当且仅当它**同时**有大残差和高杠杆。

另一个等价形式：

$$
D_i = \frac{(\hat{\mathbf{Y}}_{(i)} - \hat{\mathbf{Y}})'(\hat{\mathbf{Y}}_{(i)} - \hat{\mathbf{Y}})}{p'\hat{\sigma}^2}
$$

即所有拟合值向量因删除第 $i$ 个观测而产生的变化的平方和。

**解释**：Cook 距离没有固定的临界值，主要通过**与其他观测比较**来解读（看哪些点的 $D_i$ 明显大于其余的）。常用的经验规则是 $D_i > 1$ 或 $D_i > 4/n$ 值得关注。

---

## 6 模型设定错误（Regressionsgleichung ist nicht korrekt）

**原因**：遗漏了重要变量、包含了不必要的变量、非线性关系未捕捉、忽略了交互作用。

**后果**：参数估计和预测有**系统性偏差**。但模型仍可能提供有用的近似。

**诊断方法**：残差图（$\hat{\varepsilon}_i$ 对 $\hat{y}_i$）、F-检验检查新变量的效果、检查交互项和高阶多项式项。

**治疗方法**：扩展模型、变量变换、变量选择。

### 6.1 Partial Leverage Plot（偏杠杆图）

这是一种非常有用的图形工具，用于**在控制其他变量后**，可视化 $y$ 与某个特定变量 $x_k$ 之间的关系。

**做法**：

1. 令 $\mathbf{Q}_{(k)}$ 为不包含变量 $x_k$ 的模型的残差矩阵
2. 计算 $\mathbf{y}^* := \mathbf{Q}_{(k)}\mathbf{y}$（从 $y$ 中去除其他变量的影响）
3. 计算 $\mathbf{x}_k^* := \mathbf{Q}_{(k)}\mathbf{x}_k$（从 $x_k$ 中去除其他变量的影响）
4. 画 $\mathbf{y}^*$ 对 $\mathbf{x}_k^*$ 的散点图

这个图展示的是：**在去除所有其他变量的影响后**，$y$ 和 $x_k$ 之间剩余的关系。如果这个关系是非线性的，说明需要对 $x_k$ 做变换。

---

## 7 多重共线性（Kollinearität）

**定义**：$\mathbf{X}$ 的列（近似）线性相关。

**原因**：自变量之间高度相关、实验设计不当、离散变量的编码方式。

**后果**：

| 影响对象 | 状态 |
|---------|------|
| $\hat{\boldsymbol{\beta}}$ 无偏性 | ✅ 仍然无偏 |
| $\hat{\boldsymbol{\beta}}$ 精度 | ❌ 方差极大，估计不稳定 |
| 符号 | ❌ 经常出现**错误的正负号** |
| 置信区间 | ✅ 形式上正确，但**非常宽** |

> **直觉**：想象两个几乎平行的变量 $x_1$ 和 $x_2$。数据无法区分"是 $x_1$ 的效应还是 $x_2$ 的效应"，所以两者的系数估计都极不稳定。但它们的**线性组合**（如 $\beta_1 + \beta_2$）的估计可能仍然精确。

### 7.1 诊断工具

**a) 条件数（Konditionszahl）**

$$
K(\mathbf{X}) := \sqrt{\frac{\lambda_{\max}}{\lambda_{\min}}}
$$

其中 $\lambda_{\max}, \lambda_{\min}$ 是 $\mathbf{X}'\mathbf{X}$ 的最大和最小特征值。

| $K(\mathbf{X})$ 的值 | 判断 |
|---------------------|------|
| $< 10$ | 共线性不严重 |
| $10 \sim 30$ | 中等程度的共线性 |
| $> 30$ | 严重的共线性 ⚠️ |

**b) 方差膨胀因子（VIF）**

$$
\text{VIF}_j := \frac{1}{1 - R_j^2}
$$

其中 $R_j^2$ 是用**所有其他自变量**回归 $x_j$ 得到的 $R^2$。

**与参数方差的关系**：

$$
\sigma_{\hat{\beta}_j} = \frac{\sigma^2}{(\mathbf{x}_j - \bar{\mathbf{x}}_j)'(\mathbf{x}_j - \bar{\mathbf{x}}_j)} \cdot \text{VIF}_j
$$

所以 VIF 直接告诉你：**由于共线性，$\beta_j$ 的方差被膨胀了多少倍**。

| $\text{VIF}_j$ 的值 | 含义 |
|--------------------|------|
| 1 | 与其他变量完全不相关 |
| 5 | 方差是无共线性时的 5 倍 |
| 10 | 一般认为是临界值 ⚠️ |
| $\infty$ | 完全共线性（$R_j^2 = 1$） |

**注意**：VIF 也适用于分类变量（类别变量）。

---

## 8 测量误差（Messfehler）

**原因**：测量仪器误差、调查问卷中的不正确回答等。

**后果（对 $x$ 的加性测量误差）**：

| 影响 | 说明 |
|------|------|
| $\beta_i$ 被系统性**低估**（绝对值） | 这叫做 **attenuation bias（向零偏误）** |
| F-检验的 Power 降低 | 效应被削弱了 |

> **向零偏误的直觉**：如果 $x$ 有测量误差，观测到的 $x$ 的变异中有一部分是纯噪声。回归将效应"稀释"到了更大的 $x$ 变异范围上，所以斜率估计偏小。

**诊断**：重复测量

**治疗**：测量误差校正方法 → **Errors-in-Variables 模型**理论

---

## 9 诊断问题全景总结

| 问题 | $\hat{\boldsymbol{\beta}}$ 无偏？ | $\hat{\boldsymbol{\beta}}$ 最优？ | CI/检验正确？ | 主要诊断工具 | 主要治疗 |
|------|:----:|:----:|:----:|------|------|
| 非正态 | ✅ | ✅ | ⚠️ 小样本 | Q-Q Plot | 变换 $Y$、GLM |
| 异方差 | ✅ | ❌ | ❌ | 残差 vs $\hat{Y}$ | WLS、变换 $Y$ |
| 误差相关 | ✅ | ❌ | ❌❌ | Durbin-Watson | 时序方法、WLS |
| 异常值/高杠杆 | ⚠️ | ⚠️ | ⚠️ | $h_{ii}$、Cook's $D$、$r_i^*$ | 稳健回归、敏感性分析 |
| 模型设定错误 | ❌ | ❌ | ❌ | 残差图、Partial Leverage | 扩展模型、变换 |
| 共线性 | ✅ | ✅（但方差大） | ✅（但很宽） | VIF、条件数 | 合并/删除变量、Ridge |
| 测量误差($x$) | ❌ | ❌ | ❌ | 重复测量 | EIV 模型 |