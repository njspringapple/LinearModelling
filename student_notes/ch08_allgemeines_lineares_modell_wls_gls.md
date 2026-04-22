# 8. Das allgemeine lineare Modell（一般线性模型）

**Lineare Modelle SoSe 2025 — Sabine Hoffmann, LMU**

加权最小二乘、自相关与异方差误差项

理论推导 + 直觉理解 + 每一步都搞懂

---

## 章节导航

| 主题 | 对应内容 | 关键领悟 |
|------|---------|---------|
| 1 | 异方差模型（对角 $V$） | 每个观测的误差方差不同 |
| 2 | 模型变换的核心思想 | 除以 $\sqrt{v_i}$ 使方差齐性化 |
| 3 | 加权 KQ 估计量 | 用 $V^{-1}$ 加权，方差小的观测权重大 |
| 4 | 一般协方差结构（任意 $V$） | 误差可以相关（如时间序列） |
| 5 | 广义 Gauss-Markov 定理 | $\hat{\beta}_W$ 是 BLUE |
| 6 | 协方差结构的典型例子 | 分组异方差、AR(1)、纵向数据 |
| 7 | ML 与 REML 估计 | 当 $V$ 未知时如何估计协方差参数 |
| 8 | 关于 $\beta$ 的推断 | 近似 $t$-检验与自由度估计 |

---

## 1 出发点：为什么需要一般化？

标准线性模型假设：

$$
\mathbf{Y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon}, \quad \boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{I})
$$

即所有误差项**等方差**且**不相关**。但现实中经常违反：

| 情境 | 违反了什么 | 例子 |
|------|---------|------|
| 分组数据，各组方差不同 | 等方差 | 不同医院的治疗效果 |
| 计数数据 | 等方差（方差随均值变） | 事故次数 |
| 时间序列 | 不相关 | 月度经济数据 |
| 纵向数据（同一个体重复测量） | 不相关 | 临床试验中的重复测量 |

**解决方案**：允许误差的协方差矩阵不再是 $\sigma^2 \mathbf{I}$，而是 $\sigma^2 \mathbf{V}$。

---

## 2 异方差模型（Das allgemeine lineare Modell — 对角情形）

### 2.1 模型设定

$$
\mathbf{Y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon}
$$

$$
\boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{V})
$$

$$
\mathbf{V} = \text{diag}(v_1, v_2, \ldots, v_n)
$$

这里 $\mathbf{V}$ 是**已知的**对角矩阵，描述各观测的相对方差。

第 $i$ 个观测的方差为 $\text{Var}(\varepsilon_i) = v_i \sigma^2$。

**权重矩阵（Gewichtsmatrix）**：

$$
\mathbf{W} = \mathbf{V}^{-1} = \text{diag}(v_1^{-1}, v_2^{-1}, \ldots, v_n^{-1})
$$

> **直觉**：$v_i$ 大的观测方差大、不可靠，应该给低权重 $w_i = 1/v_i$；$v_i$ 小的观测方差小、可靠，给高权重。

---

### 2.2 变换的核心思想（Transformation）

**基本想法**：把异方差模型**变换**成等方差模型，然后用标准 OLS。

原始模型：

$$
y_i = \beta_0 \cdot 1 + \beta_1 x_{i1} + \cdots + \beta_p x_{ip} + \varepsilon_i, \quad \text{Var}(\varepsilon_i) = v_i \sigma^2
$$

两边除以 $\sqrt{v_i}$：

$$
\frac{y_i}{\sqrt{v_i}} = \beta_0 \cdot \frac{1}{\sqrt{v_i}} + \beta_1 \frac{x_{i1}}{\sqrt{v_i}} + \cdots + \beta_p \frac{x_{ip}}{\sqrt{v_i}} + \frac{\varepsilon_i}{\sqrt{v_i}}
$$

此时：

$$
\text{Var}\!\left(\frac{\varepsilon_i}{\sqrt{v_i}}\right) = \frac{v_i \sigma^2}{v_i} = \sigma^2 \quad \checkmark
$$

变换后的误差项**等方差**了！

> ⚠️ **注意**：变换后截距项的"自变量"变成了 $1/\sqrt{v_i}$，不再是常数 1。所以变换后的模型**没有通常意义上的截距**。

---

### 2.3 矩阵形式的变换

定义 $\mathbf{W}^{1/2}$ 为 $\mathbf{W}$ 的（对角）平方根矩阵，即 $\mathbf{W}^{1/2}(\mathbf{W}^{1/2})' = \mathbf{W}$。

变换后的量：

$$
\mathbf{Y}^* := \mathbf{W}^{1/2}\mathbf{Y}
$$

$$
\mathbf{X}^* := \mathbf{W}^{1/2}\mathbf{X}
$$

$$
\boldsymbol{\varepsilon}^* := \mathbf{W}^{1/2}\boldsymbol{\varepsilon}
$$

则：

$$
\mathbf{Y}^* = \mathbf{X}^*\boldsymbol{\beta} + \boldsymbol{\varepsilon}^*
$$

$$
\boldsymbol{\varepsilon}^* \sim N(\mathbf{0}, \sigma^2 \mathbf{I})
$$

> **关键领悟**：变换后就是一个**标准线性模型**！所有之前学过的理论（OLS、F-检验、$t$-检验、置信区间……）都**直接适用于变换后的模型**。参数 $\boldsymbol{\beta}$ 不变！

---

## 3 加权 KQ 估计量（Der gewichtete KQ-Schätzer）

对变换后的模型应用标准 OLS：

$$
\hat{\boldsymbol{\beta}}_W := (\mathbf{X}^{*'}\mathbf{X}^*)^{-1}\mathbf{X}^{*'}\mathbf{Y}^*
$$

展开：

$$
\hat{\boldsymbol{\beta}}_W = (\mathbf{X}'\mathbf{W}^{1/2'}\mathbf{W}^{1/2}\mathbf{X})^{-1}\mathbf{X}'\mathbf{W}^{1/2'}\mathbf{W}^{1/2}\mathbf{Y}
$$

$$
= (\mathbf{X}'\mathbf{V}^{-1}\mathbf{X})^{-1}\mathbf{X}'\mathbf{V}^{-1}\mathbf{Y}
$$

**方差估计**：

$$
\hat{\sigma}^2 = \frac{\hat{\boldsymbol{\varepsilon}}'\mathbf{V}^{-1}\hat{\boldsymbol{\varepsilon}}}{n - p'}
$$

### 3.1 $\hat{\boldsymbol{\beta}}_W$ 最小化什么？

$\hat{\boldsymbol{\beta}}_W$ 最小化的是**加权残差平方和**：

$$
(\mathbf{Y} - \mathbf{X}\boldsymbol{\beta})'\mathbf{V}^{-1}(\mathbf{Y} - \mathbf{X}\boldsymbol{\beta}) = \sum_{i=1}^n \frac{1}{v_i}(y_i - \mathbf{x}_i'\boldsymbol{\beta})^2
$$

**直觉解读**：这就是给每个观测的残差平方乘以权重 $w_i = 1/v_i$。方差大的观测（$v_i$ 大）权重小，方差小的观测权重大。

> **类比**：想象你在做加权平均。你会更信任方差小的测量值，给它们更大的权重。加权最小二乘做的就是同样的事情。

### 3.2 与普通 OLS 的对比

| | 普通 OLS | 加权 KQ |
|---|---------|--------|
| 最小化 | $\sum (y_i - \mathbf{x}_i'\boldsymbol{\beta})^2$ | $\sum \frac{1}{v_i}(y_i - \mathbf{x}_i'\boldsymbol{\beta})^2$ |
| 估计量 | $(\mathbf{X}'\mathbf{X})^{-1}\mathbf{X}'\mathbf{Y}$ | $(\mathbf{X}'\mathbf{V}^{-1}\mathbf{X})^{-1}\mathbf{X}'\mathbf{V}^{-1}\mathbf{Y}$ |
| 方差 | $\sigma^2(\mathbf{X}'\mathbf{X})^{-1}$ | $\sigma^2(\mathbf{X}'\mathbf{V}^{-1}\mathbf{X})^{-1}$ |
| $\sigma^2$ 估计 | $\hat{\boldsymbol{\varepsilon}}'\hat{\boldsymbol{\varepsilon}}/(n-p')$ | $\hat{\boldsymbol{\varepsilon}}'\mathbf{V}^{-1}\hat{\boldsymbol{\varepsilon}}/(n-p')$ |
| 适用条件 | $\text{Var}(\boldsymbol{\varepsilon}) = \sigma^2\mathbf{I}$ | $\text{Var}(\boldsymbol{\varepsilon}) = \sigma^2\mathbf{V}$ |

当 $\mathbf{V} = \mathbf{I}$ 时，加权 KQ 退化为普通 OLS。

---

## 4 一般协方差结构（Verallgemeinerte KQ-Methode）

### 4.1 模型

$$
\mathbf{Y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon}
$$

$$
\boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{V})
$$

这里 $\mathbf{V} \in \mathbb{R}^{n \times n}$ 是**任意的已知正定矩阵**，不再要求是对角矩阵。这意味着误差项之间**可以相关**。

### 4.2 变换

定义 $\mathbf{W} := \mathbf{V}^{-1}$。存在可逆矩阵 $\mathbf{W}^{1/2}$（例如通过 Cholesky 分解），使得：

$$
\mathbf{W}^{1/2}(\mathbf{W}^{1/2})' = \mathbf{W}
$$

同样定义变换后的量：

$$
\mathbf{Y}^* := \mathbf{W}^{1/2}\mathbf{Y}, \quad \mathbf{X}^* := \mathbf{W}^{1/2}\mathbf{X}, \quad \boldsymbol{\varepsilon}^* := \mathbf{W}^{1/2}\boldsymbol{\varepsilon}
$$

验证：

$$
\text{Var}(\boldsymbol{\varepsilon}^*) = \mathbf{W}^{1/2} \cdot \sigma^2\mathbf{V} \cdot (\mathbf{W}^{1/2})' = \sigma^2 \mathbf{W}^{1/2} \mathbf{V} (\mathbf{W}^{1/2})' = \sigma^2 \mathbf{W}^{1/2} \mathbf{W}^{-1} (\mathbf{W}^{1/2})' 
$$

由于 $\mathbf{V} = \mathbf{W}^{-1}$，而 $\mathbf{W}^{1/2}(\mathbf{W}^{1/2})' = \mathbf{W}$，所以 $\mathbf{W}^{1/2}\mathbf{W}^{-1}(\mathbf{W}^{1/2})' = \mathbf{W}^{1/2}(\mathbf{W}^{1/2})'^{-1} \cdot (\mathbf{W}^{1/2})^{-1} \cdot (\mathbf{W}^{1/2})' $。更直接地：

$$
\text{Var}(\boldsymbol{\varepsilon}^*) = \sigma^2 \mathbf{W}^{1/2}\mathbf{V}(\mathbf{W}^{1/2})' = \sigma^2 \mathbf{I} \quad \checkmark
$$

### 4.3 广义 KQ 估计量

$$
\hat{\boldsymbol{\beta}}_W = (\mathbf{X}'\mathbf{V}^{-1}\mathbf{X})^{-1}\mathbf{X}'\mathbf{V}^{-1}\mathbf{Y}
$$

$$
\hat{\sigma}^2 = \frac{\hat{\boldsymbol{\varepsilon}}'\mathbf{V}^{-1}\hat{\boldsymbol{\varepsilon}}}{n - p'}
$$

形式上和异方差情形完全相同。唯一的区别是 $\mathbf{V}$ 现在可以有非零的非对角元素（即误差之间可以相关）。

---

## 5 性质与广义 Gauss-Markov 定理

### 5.1 基本性质

$$
E(\hat{\boldsymbol{\beta}}_W) = \boldsymbol{\beta}
$$

$$
\text{Var}(\hat{\boldsymbol{\beta}}_W) = \sigma^2(\mathbf{X}'\mathbf{V}^{-1}\mathbf{X})^{-1}
$$

$\hat{\boldsymbol{\beta}}_W$ 同时也是 **ML 估计量**（在正态假设下）。

所有标准的检验方法和平方和分解都可以在变换后的模型 $\mathbf{Y}^* = \mathbf{X}^*\boldsymbol{\beta} + \boldsymbol{\varepsilon}^*$ 中进行，从而**归结为等方差情形**。

### 5.2 广义 Gauss-Markov 定理（Allgemeines Gauss-Markov-Theorem）

**定理**：给定模型

$$
\mathbf{Y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon}, \quad \text{rg}(\mathbf{X}) = p', \quad E(\boldsymbol{\varepsilon}) = \mathbf{0}, \quad \text{Var}(\boldsymbol{\varepsilon}) = \sigma^2\mathbf{V}
$$

其中 $\mathbf{V}$ 已知。则 $\hat{\boldsymbol{\beta}}_W = (\mathbf{X}'\mathbf{V}^{-1}\mathbf{X})^{-1}\mathbf{X}'\mathbf{V}^{-1}\mathbf{Y}$ 是所有线性无偏估计量中**方差最小的**（即 BLUE）。

> **注意**：这个定理**不需要正态假设**，只需要 $E(\boldsymbol{\varepsilon}) = \mathbf{0}$ 和 $\text{Var}(\boldsymbol{\varepsilon}) = \sigma^2\mathbf{V}$。正态假设只是为了推导精确的 $t$- 和 $F$-分布。

> **直觉**：如果你知道了误差的协方差结构但忽略它（仍然用 OLS），你的估计仍然无偏，但**不再是最有效的**。加权 KQ 利用了协方差结构的信息，因此效率更高。

### 5.3 如果用错了权重会怎样？

| 情况 | 无偏？ | BLUE？ | CI/检验正确？ |
|------|:----:|:----:|:----:|
| 用正确的 $\mathbf{V}$ 做 WLS | ✅ | ✅ | ✅ |
| 忽略 $\mathbf{V}$，用普通 OLS | ✅ | ❌ | ❌ |
| 用错误的 $\mathbf{V}$ 做 WLS | ✅ | ❌ | ❌ |

OLS 在异方差/自相关下仍然无偏（Gauss-Markov 不需要等方差），但标准误和检验都会出错。

---

## 6 协方差结构的典型例子

### 6.1 分组异方差

各组的方差不同，但组内等方差且观测之间独立：

$$
\mathbf{V} = \text{diag}(\underbrace{v_1, \ldots, v_1}_{n_1 \text{ 个}}, \underbrace{v_2, \ldots, v_2}_{n_2 \text{ 个}}, \ldots)
$$

实践中 $v_k$ 可以通过各组的残差方差来估计。

### 6.2 AR(1) 结构（时间序列）

误差项服从一阶自回归过程：

$$
\varepsilon_i = \rho \varepsilon_{i-1} + \eta_i
$$

其中 $\eta_i$ 独立同分布，$\text{Var}(\eta_i) = \sigma^2$，$\text{Var}(\varepsilon_1) = \sigma^2/(1-\rho^2)$。

协方差矩阵为：

$$
\mathbf{V} = \frac{\sigma^2}{1 - \rho^2} \begin{pmatrix} 1 & \rho & \rho^2 & \cdots & \rho^{n-1} \\ \rho & 1 & \rho & \cdots & \rho^{n-2} \\ \rho^2 & \rho & 1 & \cdots & \rho^{n-3} \\ \vdots & & & \ddots & \vdots \\ \rho^{n-1} & \rho^{n-2} & \rho^{n-3} & \cdots & 1 \end{pmatrix}
$$

**性质**：

| 特征 | 说明 |
|------|------|
| $\text{Corr}(\varepsilon_i, \varepsilon_{i+k}) = \rho^k$ | 相关性随距离指数衰减 |
| $\|\rho\| < 1$ | 平稳性条件 |
| $\rho > 0$ | 正自相关（最常见） |
| $\rho = 0$ | 退化为标准模型 |

> **直觉**：如果这个月的误差偏高，下个月的误差也倾向于偏高（$\rho > 0$）。随着时间间隔增大，相关性衰减。

### 6.3 纵向数据（块对角结构）

每个个体有 $T$ 次重复测量，个体之间独立，但同一个体内部的测量相关：

$$
\mathbf{V} = \begin{pmatrix} \mathbf{V}_1 & \mathbf{0} & \cdots & \mathbf{0} \\ \mathbf{0} & \mathbf{V}_2 & \cdots & \mathbf{0} \\ \vdots & & \ddots & \vdots \\ \mathbf{0} & \mathbf{0} & \cdots & \mathbf{V}_m \end{pmatrix}
$$

每个 $\mathbf{V}_k$ 是 $T \times T$ 的协方差矩阵。

### 6.4 对称结构 / 复合对称（Compound Symmetry）

常见于混合效应模型。同一组内任意两个观测的相关系数相同：

$$
\mathbf{V}_k = \begin{pmatrix} 1 & \rho & \rho & \cdots & \rho \\ \rho & 1 & \rho & \cdots & \rho \\ \vdots & & & \ddots & \vdots \\ \rho & \rho & \rho & \cdots & 1 \end{pmatrix}
$$

这等价于假设组内有一个共同的随机效应（random intercept）。

---

## 7 当 $\mathbf{V}$ 未知时：ML 与 REML 估计

### 7.1 问题

之前假设 $\mathbf{V}$ 已知，但实践中 $\mathbf{V}$ 通常依赖于未知参数 $\boldsymbol{\vartheta}$（如 AR(1) 中的 $\rho$，或分组异方差中的各组方差）。

写作 $\mathbf{V} = \mathbf{V}(\boldsymbol{\vartheta})$，需要同时估计 $\boldsymbol{\beta}$ 和 $\boldsymbol{\vartheta}$。

### 7.2 对数似然函数

$$
l(\boldsymbol{\beta}, \boldsymbol{\vartheta}) = -\frac{1}{2}\left(\ln|\mathbf{V}(\boldsymbol{\vartheta})| + (\mathbf{Y} - \mathbf{X}\boldsymbol{\beta})'\mathbf{V}^{-1}(\boldsymbol{\vartheta})(\mathbf{Y} - \mathbf{X}\boldsymbol{\beta})\right)
$$

（省略了与参数无关的常数项）

### 7.3 两步策略：Profile Likelihood

**步骤 1**：如果 $\boldsymbol{\vartheta}$ 已知，$\boldsymbol{\beta}$ 的 MLE 就是加权 KQ 估计量：

$$
\hat{\boldsymbol{\beta}}(\boldsymbol{\vartheta}) = \left(\mathbf{X}'\mathbf{V}(\boldsymbol{\vartheta})^{-1}\mathbf{X}\right)^{-1}\mathbf{X}'\mathbf{V}^{-1}(\boldsymbol{\vartheta})\mathbf{Y}
$$

**步骤 2**：将 $\hat{\boldsymbol{\beta}}(\boldsymbol{\vartheta})$ 代入似然函数，得到**Profile Log-Likelihood**（只关于 $\boldsymbol{\vartheta}$ 的函数）：

$$
l(\boldsymbol{\vartheta}) = -\frac{1}{2}\left(\ln|\mathbf{V}(\boldsymbol{\vartheta})| + (\mathbf{Y} - \mathbf{X}\hat{\boldsymbol{\beta}}(\boldsymbol{\vartheta}))'\mathbf{V}^{-1}(\boldsymbol{\vartheta})(\mathbf{Y} - \mathbf{X}\hat{\boldsymbol{\beta}}(\boldsymbol{\vartheta}))\right)
$$

**步骤 3**：最大化 $l(\boldsymbol{\vartheta})$ 得到 $\hat{\boldsymbol{\vartheta}}_{\text{ML}}$，然后 $\hat{\boldsymbol{\beta}} = \hat{\boldsymbol{\beta}}(\hat{\boldsymbol{\vartheta}}_{\text{ML}})$。

### 7.4 REML（Restricted Maximum Likelihood）

**问题**：ML 估计 $\boldsymbol{\vartheta}$（特别是方差参数）时往往**有偏**——系统性地低估方差。

> **类比**：就像普通线性模型中 $\hat{\sigma}^2_{\text{ML}} = \text{RSS}/n$（有偏），而 $\hat{\sigma}^2 = \text{RSS}/(n-p')$（无偏）。ML 没有考虑到估计 $\boldsymbol{\beta}$ 时"用掉了" $p'$ 个自由度。

**REML 的解决方案**：最大化一个修正后的似然函数：

$$
L_R(\boldsymbol{\vartheta}) = l(\boldsymbol{\vartheta}) - \frac{1}{2}\ln|\mathbf{X}'\mathbf{V}(\boldsymbol{\vartheta})^{-1}\mathbf{X}|
$$

第二项 $-\frac{1}{2}\ln|\mathbf{X}'\mathbf{V}(\boldsymbol{\vartheta})^{-1}\mathbf{X}|$ 是一个**修正项**，补偿了因估计 $\boldsymbol{\beta}$ 而损失的自由度。

**在简单线性模型中**（$\mathbf{V} = \mathbf{I}$），REML 估计恰好给出 $\hat{\sigma}^2 = \text{RSS}/(n-p')$，即**无偏估计**。

### 7.5 ML vs REML 对比

| | ML | REML |
|---|------|------|
| 估计 $\boldsymbol{\vartheta}$ | 有偏（低估方差） | 近似无偏 |
| 小样本表现 | 较差 | **更好** |
| 模型比较（不同固定效应） | ✅ 可以用似然比检验 | ❌ 不可以（不同 $\mathbf{X}$ 时 REML 不可比） |
| 模型比较（不同协方差结构） | ✅ | ✅ |
| 实践中 | 较少使用 | **更常用**（推荐） |

> **经验法则**：
> - 比较**不同的协方差结构**（相同固定效应）→ 用 REML
> - 比较**不同的固定效应**（相同协方差结构）→ 用 ML

---

## 8 关于 $\boldsymbol{\beta}$ 的推断

### 8.1 $\hat{\boldsymbol{\beta}}$ 的分布

如果 $\boldsymbol{\vartheta}$ **已知**：

$$
\hat{\boldsymbol{\beta}}(\boldsymbol{\vartheta}) \sim N\left(\boldsymbol{\beta},\; \sigma^2(\mathbf{X}'\mathbf{V}^{-1}\mathbf{X})^{-1}\right)
$$

这是**精确的**正态分布，所有标准推断（$t$-检验、$F$-检验、置信区间）都**精确成立**。

### 8.2 $\boldsymbol{\vartheta}$ 未知时的困难

当 $\mathbf{V}$ 用 (RE)ML 估计量 $\mathbf{V}(\hat{\boldsymbol{\vartheta}})$ 替代时：

$$
\text{Var}(\hat{\boldsymbol{\beta}}) \approx \sigma^2(\mathbf{X}'\mathbf{V}(\hat{\boldsymbol{\vartheta}})^{-1}\mathbf{X})^{-1}
$$

但这只是**近似的**，因为 $\hat{\boldsymbol{\vartheta}}$ 本身也有不确定性。

**后果**：

$$
\frac{\hat{\beta}_j - \beta_j}{\widehat{\text{s.e.}}(\hat{\beta}_j)}
$$

**不再精确服从** $t$ 分布。

### 8.3 实践中的解决方案：近似 $t$-检验

在实践中，人们仍然用 $t$-检验的形式：

$$
T_j = \frac{\hat{\beta}_j - \beta_j}{\widehat{\text{s.e.}}(\hat{\beta}_j)} \;\dot{\sim}\; t(\nu)
$$

但需要**估计自由度 $\nu$**。常用的方法包括 Satterthwaite 近似和 Kenward-Roger 近似。

| 方法 | 说明 | 适用 |
|------|------|------|
| Satterthwaite | 基于方差估计的矩匹配 | R 中 `lmerTest` 包 |
| Kenward-Roger | 对协方差矩阵做小样本修正 | 更精确，但计算更复杂 |

> **注意**：对于特殊模型（如平衡的方差分析），渐近正态性已被严格证明。但**一般情况下**，渐近正态性的严格证明尚不完整。实践中需要依赖近似方法。

---

## 9 全章逻辑总结

整个第 8 章的逻辑链条：

**问题**：标准模型假设 $\text{Var}(\boldsymbol{\varepsilon}) = \sigma^2\mathbf{I}$ 太强，实际中经常违反。

**解决方案**：

1. 允许 $\text{Var}(\boldsymbol{\varepsilon}) = \sigma^2\mathbf{V}$

2. 如果 $\mathbf{V}$ **已知**：
   - 用 $\mathbf{W}^{1/2}$ 变换模型 → 标准模型
   - 加权 KQ 估计量 $\hat{\boldsymbol{\beta}}_W = (\mathbf{X}'\mathbf{V}^{-1}\mathbf{X})^{-1}\mathbf{X}'\mathbf{V}^{-1}\mathbf{Y}$
   - 广义 Gauss-Markov 定理保证 BLUE
   - 所有标准推断（$t$, $F$, CI）在变换后的模型中精确成立

3. 如果 $\mathbf{V}$ **未知**（$\mathbf{V} = \mathbf{V}(\boldsymbol{\vartheta})$）：
   - 用 ML 或 REML 估计 $\boldsymbol{\vartheta}$
   - REML 在小样本中更可靠
   - 推断是**近似的**（需要估计自由度）

| 步骤 | $\mathbf{V}$ 已知 | $\mathbf{V}$ 未知 |
|------|:----:|:----:|
| 估计 $\boldsymbol{\beta}$ | $(\mathbf{X}'\mathbf{V}^{-1}\mathbf{X})^{-1}\mathbf{X}'\mathbf{V}^{-1}\mathbf{Y}$ | 同左，但用 $\mathbf{V}(\hat{\boldsymbol{\vartheta}})$ 替代 |
| $\hat{\boldsymbol{\beta}}$ 是 BLUE？ | ✅（精确） | ⚠️ 近似 |
| 推断精确？ | ✅ | ⚠️ 近似（需要 Satterthwaite 等） |
| $\hat{\sigma}^2$ 无偏？ | ✅ | ML: ❌ / REML: ≈ ✅ |