# 6. metrische Einflussgrößen（连续型自变量的处理）

## 章节导航

| 主题 | 对应内容 | 关键领悟 |
|------|---------|---------|
| 1 | 连续变量的交互作用 | 效果不再是常数，而是随另一个变量变化 |
| 2 | 连续变量的 8 种处理方式 | 从简单线性到 spline，灵活度逐步递增 |
| 3 | 通用基函数框架 | 所有方法都可以写成 $\beta_0 + \sum \beta_j B_j(x)$ |
| 4 | 多重检验与同时置信区间 | Bonferroni、Scheffé、Tukey/Hothorn |
| 5 | 置信椭球体 | 多维参数的联合置信域 |

---

## 1 连续变量的交互作用（Interaktionen bei metrischen Variablen）

### 1.1 模型

两个连续变量 $x_1, x_2$ 加上交互项：

$$
E(Y) = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \beta_3 x_1 \cdot x_2
$$

### 1.2 怎么理解？两种等价视角

**视角 A**：固定 $x_1$，看 $x_2$ 的效果

$$
E(Y) = \underbrace{(\beta_0 + \beta_1 x_1)}_{\text{截距随 } x_1 \text{ 变}} + \underbrace{(\beta_2 + \beta_3 x_1)}_{\text{斜率随 } x_1 \text{ 变}} \cdot x_2
$$

→ $x_2$ 对 $Y$ 的边际效应**不是常数**，而是 $\beta_2 + \beta_3 x_1$，取决于 $x_1$ 的值。

**视角 B**：固定 $x_2$，看 $x_1$ 的效果

$$
E(Y) = (\beta_0 + \beta_2 x_2) + (\beta_1 + \beta_3 x_2) \cdot x_1
$$

→ $x_1$ 对 $Y$ 的边际效应是 $\beta_1 + \beta_3 x_2$。

### 1.3 ⚠️ 解释陷阱

$\beta_1$ 的含义是"**当 $x_2 = 0$ 时**，$x_1$ 的斜率"。$\beta_2$ 同理是"**当 $x_1 = 0$ 时**，$x_2$ 的斜率"。

问题：$x_2 = 0$ 或 $x_1 = 0$ 往往没有实际意义！

> **解决方案：中心化（Zentrierung）**
>
> 将变量替换为 $\tilde{x}_1 = x_1 - \bar{x}_1$，$\tilde{x}_2 = x_2 - \bar{x}_2$。
> 这样 $\beta_1$ 就是"当 $x_2$ 取**均值**时，$x_1$ 的斜率"——有意义得多！

### 1.4 直觉图示

想象一个三维曲面 $E(Y)$ 对 $(x_1, x_2)$：

| 模型 | 曲面形状 |
|------|---------|
| 无交互项（$\beta_3 = 0$） | 平面，各方向斜率固定 |
| 有交互项（$\beta_3 \neq 0$） | 扭曲的曲面，斜率随位置变化 |

---

## 2 连续型自变量的 8 种处理方式

这是本章最核心的内容。从简单到复杂，列出了**8 种方法**来建模 $x$ 对 $Y$ 的非线性效应。

### 方式 ① 简单线性

$$
E(Y) = \beta_0 + \beta_1 x
$$

最简单，假设效应是一条直线。适用于效应确实是线性的情况。

---

### 方式 ② 变量变换

$$
E(Y) = \beta_0 + \beta_1 T(x)
$$

对 $x$ 做一个已知的变换 $T(\cdot)$，然后线性回归。常见选择：

| 变换 | 公式 | 适用情况 |
|------|------|---------|
| 对数 | $T(x) = \ln(x)$ | 效应递减、正偏分布 |
| 保零对数 | $T(x) = \ln(1 + x)$ | $x$ 可以取 0 |
| 幂变换 | $T(x) = x^c$（$c$ 已知） | 如 $\sqrt{x}$（$c = 0.5$） |

**关键点**：$\beta_1$ 的解释**随变换而变**！

例如 $T(x) = \ln(x)$：$\beta_1$ 表示 "$x$ 增加 1%（而非增加 1 个单位）时 $Y$ 的期望变化约为 $\beta_1 / 100$"。

> **仍然是线性模型！** 因为对参数 $\beta$ 仍然是线性的：$E(Y) = \beta_0 + \beta_1 \cdot [\text{新变量}]$。
> 只是"新变量"不再是 $x$，而是 $T(x)$。

---

### 方式 ③ 多项式回归（Polynomiale Regression）

$$
E(Y) = \beta_0 + \beta_1 x + \beta_2 x^2 + \beta_3 x^3 + \cdots + \beta_k x^k
$$

用 $x$ 的幂次来拟合曲线。

**核心问题**：$k$ 选多大？

**解决思路**：逐步检验——用**序贯平方和（Sequentielle Quadratsummen）**

从高次项开始检验 $H_0: \beta_l = \cdots = \beta_k = 0$，如果不显著就降低次数。

**优点**：简单直观
**缺点**：

- 高次多项式在数据边缘容易剧烈振荡（Runge 现象）
- $x, x^2, x^3, \ldots$ 之间高度共线
- 实践中一般不超过 3 次或 4 次

---

### 方式 ④ 分段常数函数（Stückweise konstante Funktion）

$$
E(Y) = \begin{cases} \beta_0 & x \leq x_0 \\ \beta_1 & x_0 < x < x_1 \\ \vdots \\ \beta_h & x > x_{h-1} \end{cases}
$$

本质上就是把连续变量**分组/分箱（Kategorisierung）**，变成类别变量。

**用途**：当 $x$ 只有少数几个取值时，可以用来**检验线性假设是否成立**。

做法：拟合分段常数模型 vs 简单线性模型，比较拟合效果。如果两者差不多，说明线性假设合理。

**缺点**：

- 信息损失（把连续变量变成离散）
- 切点（断点）的选择是主观的
- 函数不连续，跳跃式

---

### 方式 ⑤ 分段线性函数（Stückweise lineare Funktion）

$$
E(Y) = \beta_0 + \beta_1 x + \beta_2 (x - g_1)_+ + \beta_3 (x - g_2)_+ + \cdots + \beta_h (x - g_k)_+
$$

其中截断函数 $t_+ = \max(t, 0)$，$g_1, g_2, \ldots, g_k$ 是**已知的断点（Knoten / Knots）**。

**工作原理**：

- 当 $x < g_1$ 时：$E(Y) = \beta_0 + \beta_1 x$（一条直线）
- 当 $g_1 \leq x < g_2$ 时：$E(Y) = \beta_0 + \beta_1 x + \beta_2(x - g_1)$，斜率变为 $\beta_1 + \beta_2$
- 每过一个断点，斜率增加 $\beta_j$

→ 整体是一条**连续的折线**，在每个断点处斜率可以改变。

**优点**：连续，简单，可解释
**缺点**：在断点处**不可微**（有"尖角"）

---

### 方式 ⑥ 回归样条（Regressionsspline）

$$
E(Y) = \beta_0 + \beta_1 x + \beta_2 x^2 + \beta_3 x^3 + \beta_4 (x - g_1)^3_+ + \beta_5 (x - g_2)^3_+
$$

这是一个**三次样条（kubischer Spline）**：

- 基础部分 $\beta_0 + \beta_1 x + \beta_2 x^2 + \beta_3 x^3$ 是全局三次多项式
- 每个断点处加一个 $(x - g_j)^3_+$ 项来局部调整形状

**为什么用三次？**

因为 $x^3$ 在 0 点处**二阶连续可微**（$f(0) = 0, f'(0) = 0, f''(0) = 0$ 都连续）。

所以 $(x - g_j)^3_+$ 保证了**整个函数在断点处二阶连续可微**——曲线光滑，没有尖角！

| 方法 | 连续性 | 一阶可微 | 二阶可微 |
|------|--------|---------|---------|
| 分段常数 | ✗ | ✗ | ✗ |
| 分段线性 | ✓ | ✗ | ✗ |
| 三次样条 | ✓ | ✓ | ✓ |

**断点的选择**：仍然需要预先指定。实践中常用数据的分位数。

---

### 方式 ⑦ 分数多项式（Fraktionale Polynome）

例：Grad (-2, 2)

$$
E(Y) = \beta_0 + \beta_1 x^{-2} + \beta_2 x^{-1} + \beta_3 x^{-0.5} + \beta_4 \ln(x) + \beta_5 x^{0.5} + \beta_6 x + \beta_7 x^2
$$

**想法**：不再局限于整数幂 $x, x^2, x^3, \ldots$，而是允许**分数幂和对数**。

标准的幂次集合为 $\{-2, -1, -0.5, 0, 0.5, 1, 2, 3\}$，其中 $x^0$ 约定为 $\ln(x)$。

**优点**：

- 非常灵活，能拟合各种 U 形、J 形、饱和型曲线
- 有系统的选择策略（Sauerbrei et al., Uni Freiburg）
- 不需要选择断点

**缺点**：

- 需要 $x > 0$
- 解释不如简单线性直观

---

### 方式 ⑧ 三角多项式（Trigonometrische Polynome）

用于建模**周期性/季节性**效应：

$$
E(Y) = \beta_0 + \beta_1 \sin\!\left(\frac{2\pi}{T} x\right) + \beta_2 \cos\!\left(\frac{2\pi}{T} x\right) + \beta_3 \sin\!\left(\frac{2\pi}{T} \cdot 2x\right) + \beta_4 \cos\!\left(\frac{2\pi}{T} \cdot 2x\right)
$$

其中 $T$ = 周期长度，$x$ = 时间。

**数学事实**：

$$
A_1 \cos(x) + A_2 \sin(x) = A_3 \sin(x + \phi)
$$

即 sin + cos 的线性组合仍然是一个正弦波（振幅 $A_3 = \sqrt{A_1^2 + A_2^2}$，相移 $\phi$）。

**替代方案**：季节性虚拟变量（Saison-Dummy），但虚拟变量需要更多参数。

---

## 3 通用基函数框架（Allgemeiner Ansatz mit Basisfunktionen）

**统一视角**：上述所有方法都可以写成

$$
E(Y) = \beta_0 + \beta_1 B_1(x) + \beta_2 B_2(x) + \beta_3 B_3(x) + \cdots
$$

其中 $B_1, B_2, B_3, \ldots$ 是**已知的基函数（Basisfunktionen）**。

| 方法 | 基函数 $B_j(x)$ |
|------|----------------|
| 简单线性 | $B_1(x) = x$ |
| 多项式 | $B_j(x) = x^j$ |
| 分段线性 | $B_1(x) = x,\; B_{j+1}(x) = (x - g_j)_+$ |
| 三次样条 | $x, x^2, x^3, (x-g_1)^3_+, (x-g_2)^3_+, \ldots$ |
| 三角多项式 | $\sin(\cdot), \cos(\cdot)$ |

> **关键领悟**：只要基函数是已知的，模型对**参数 $\beta$** 仍然是线性的！
> → 所有标准线性模型的理论（OLS、F-检验、置信区间……）都直接适用。

**P-Splines（惩罚 B 样条）**：选择大量基函数 + 用惩罚项防止过拟合。这属于更高级的方法（广义线性模型课程内容）。

---

## 4 示例：巴登-符腾堡州狐狸种群趋势

$Y$：射杀狐狸数量（种群大小指标），$t$：时间

**模型 1**：二次多项式

$$
\ln(Y) = \beta_0 + \beta_1 t + \beta_2 t^2
$$

**模型 2**：三次样条

$$
\ln(Y) = \beta_0 + \beta_1 t + \beta_2 t^2 + \beta_3 t^3 + \beta_4 (t - 70)^3_+ + \beta_5 (t - 85)^3_+
$$

注意：因变量取了对数 $\ln(Y)$，这是对**因变量的变换**，使得原始尺度上的效应是乘性的。

---

## 5 多重检验与同时置信区间（Multiples Testen）

### 5.1 问题

在多元线性模型中，我们经常同时做**多个检验**或构建**多个置信区间**。

如果每个检验单独的显著性水平是 $\alpha = 0.05$，做 $k$ 个检验时：

$$
P(\text{至少一个错误拒绝}) \gg 0.05
$$

例如 $k = 20$ 个独立检验：$P(\text{至少一个}) = 1 - (1-0.05)^{20} \approx 0.64$。

### 5.2 三种水平的区分

| 概念 | 定义 |
|------|------|
| **局部水平（lokales Niveau）** | 单个检验的显著性水平 |
| **全局水平（globales Niveau）** | 当**所有** $H_0$ 都成立时，至少拒绝一个的概率 |
| **多重水平（multiples Niveau）** | 至少**错误地**拒绝一个 $H_0$ 的概率 $< \alpha$（FWER） |

### 5.3 四种常用策略

#### ① Bonferroni 校正

做 $k$ 个检验时，每个检验用水平 $\alpha_l = \alpha / k$。

**原理**：Bonferroni 不等式

$$
P\!\left(\bigcup E_i\right) \leq \sum P(E_i)
$$

所以 $\sum_{i=1}^k \frac{\alpha}{k} = \alpha$，保证了多重水平 $\leq \alpha$。

| 优点 | 缺点 |
|------|------|
| 极其通用，对任何检验都适用 | 当 $k$ 很大时过于保守 |
| 实施简单 | 检验力（Power）低 |
| 严格控制 FWER | — |

#### ② Scheffé 方法

适用于线性模型中**任意多个**线性组合 $\gamma = \mathbf{a}'\boldsymbol{\beta}$ 的同时置信区间。

对参数 $\beta_j$ 的 Scheffé 同时置信区间为：

$$
\hat{\beta}_j \pm \sqrt{p' \cdot F_{1-\alpha}(p', n - p')} \;\cdot\; \hat{\sigma}_{\hat{\beta}_j}
$$

对任意线性组合 $\gamma = \mathbf{a}'\boldsymbol{\beta}$：

$$
\hat{\gamma} \pm \sqrt{p' \cdot F_{1-\alpha}(p', n - p')} \;\cdot\; \hat{\sigma}_{\hat{\gamma}}
$$

其中 $p'$ 是参数个数，$n$ 是样本量。

**注意**：和普通 $t$-区间的区别仅在于用 $\sqrt{p' \cdot F_{1-\alpha}(p', n-p')}$ 代替了 $t_{1-\alpha/2}(n-p')$。当只看一个参数时，Scheffé 区间比 $t$-区间更宽（更保守），但它同时对**所有可能的线性组合**有效。

#### ③ Tukey / Hothorn et al. 方法（最大统计量）

给定同时假设 $\mathbf{A}\boldsymbol{\beta} = \mathbf{0}$，对每一行构造学生化统计量：

$$
T_i = \frac{\mathbf{d}_i'\boldsymbol{\beta}}{\hat{\sigma}_{\mathbf{d}_i'\boldsymbol{\beta}}}, \quad i = 1, \ldots, a
$$

取最大值：

$$
Q = \max_i |T_i|
$$

设 $Q(n - p', \mathbf{R})$ 是在相关矩阵 $\mathbf{R}$ 下 $t$-分布最大值的分布。找临界值 $q$ 使得 $P(Q \leq q) = 1 - \alpha$，则：

$$
P\!\left(\forall\, i:\; \mathbf{d}_i'\boldsymbol{\beta} \in \left[-q \cdot \hat{\sigma}_{\mathbf{d}_i'\boldsymbol{\beta}},\; +q \cdot \hat{\sigma}_{\mathbf{d}_i'\boldsymbol{\beta}}\right]\right) = 1 - \alpha
$$

即得到**同时置信区间**。

**优点**：

- 比 Bonferroni 更精确（利用了统计量之间的相关结构）
- 非常通用，可用于任何（渐近）正态的检验统计量
- R 实现：`multcomp` 包（Hothorn, Bretz, Westfall, 2008；被引用超过 9000 次）

#### ④ 封闭检验过程（Abschlussprozeduren）

对于检验系统（如 3 组比较），通常比其他方法**更有检验力**。

**应避免的方法**：Duncan multiple range、Fisher LSD、Newman-Keuls——这些方法**不控制多重水平**。

---

## 6 置信椭球体（Konfidenzellipsoide）

### 6.1 参数向量的联合置信域

在正态线性模型 $\mathbf{Y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon}$，$\boldsymbol{\varepsilon} \sim N(\mathbf{0}, \sigma^2 \mathbf{I})$ 下：

$$
\left\{\boldsymbol{\beta} \;\middle|\; (\boldsymbol{\beta} - \hat{\boldsymbol{\beta}})' (\mathbf{X}'\mathbf{X}) (\boldsymbol{\beta} - \hat{\boldsymbol{\beta}}) \leq p' \hat{\sigma}^2 F_{1-\alpha}(p', n - p')\right\}
$$

是 $\boldsymbol{\beta}$ 的 $1-\alpha$ 置信域。

**几何含义**：这是以 $\hat{\boldsymbol{\beta}}$ 为中心、以 $(\mathbf{X}'\mathbf{X})^{-1}$ 决定形状的**椭球体**。

### 6.2 线性变换的置信域

对于 $\boldsymbol{\gamma} = \mathbf{A}\boldsymbol{\beta}$：

$$
\left\{\boldsymbol{\gamma} \;\middle|\; (\boldsymbol{\gamma} - \hat{\boldsymbol{\gamma}})' \widehat{V(\hat{\boldsymbol{\gamma}})}^{-1} (\boldsymbol{\gamma} - \hat{\boldsymbol{\gamma}}) < \dim(\boldsymbol{\gamma}) \cdot F_{1-\alpha}(\dim(\boldsymbol{\gamma}), n - p')\right\}
$$

是 $\boldsymbol{\gamma}$ 的 $1-\alpha$ 置信域。

> **联系**：Scheffé 同时置信区间就是从这个椭球体投影到各个坐标轴上得到的。
> 椭球体内的**每一个点**对应的所有线性组合都被同时覆盖。

---

## 7 各方法对比总结

| 方法 | 灵活度 | 参数数量 | 光滑度 | 需要选择 | 仍是线性模型？ |
|------|--------|---------|--------|---------|-------------|
| ① 简单线性 | ★☆☆☆ | 2 | ∞ | 无 | ✓ |
| ② 变量变换 | ★★☆☆ | 2 | 取决于 $T$ | 变换形式 | ✓ |
| ③ 多项式 | ★★★☆ | $k+1$ | $C^\infty$ | 阶数 $k$ | ✓ |
| ④ 分段常数 | ★★☆☆ | $h+1$ | 不连续 | 断点 | ✓ |
| ⑤ 分段线性 | ★★★☆ | $k+2$ | $C^0$ | 断点 | ✓ |
| ⑥ 三次样条 | ★★★★ | $k+4$ | $C^2$ | 断点 | ✓ |
| ⑦ 分数多项式 | ★★★★ | 灵活 | 取决于幂次 | 幂次集合 | ✓ |
| ⑧ 三角多项式 | ★★★☆ | $2m+1$ | $C^\infty$ | 周期 $T$、阶数 | ✓ |

> **最重要的结论**：所有这些方法都是**对参数线性的**，因为基函数 $B_j(x)$ 是已知的。
> 所以 OLS 估计、F-检验、$t$-检验、置信区间——**全部照常使用**！