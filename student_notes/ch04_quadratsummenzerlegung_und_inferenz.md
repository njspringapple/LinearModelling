
---

## 🎯 贯穿全章的例子：Lesetest

**研究问题**：什么因素影响小学生的阅读错误数？

**数据**：180 个 3-4 年级学生，变量：

- Y = 错误数
- GE = 性别（1=男，0=女）
- JG = 年级（1=4年级，0=3年级）
- LZ = 学校阅读时间（分钟）
- WOL = 课外阅读频率（0-4）
- WOG = 打游戏机频率（0-4）
- WOTV = 看电视频率（0-4）

**完整模型**：

$$
Y = \beta_0 + \beta_1 \text{GE} + \beta_2 \text{JG} + \beta_3 \text{LZ} + \beta_4 \text{WOL} + \beta_5 \text{WOG} + \beta_6 \text{WOTV} + \varepsilon
$$

所以 $p=6$ 个斜率，$p'=p+1=7$ 个参数，$n=180$。

---

## 1. 平方和分解（Quadratsummenzerlegung）

### 1.1 核心恒等式

$$
\underbrace{(Y - \bar{Y})'(Y - \bar{Y})}_{SST} = \underbrace{(Y - \hat{Y})'(Y - \hat{Y})}_{SSE} + \underbrace{(\hat{Y} - \bar{Y})'(\hat{Y} - \bar{Y})}_{SSM}
$$

**直觉翻译**：

|符号|中文含义|回答的问题|
|---|---|---|
|**SST**|总变异|"如果我啥模型都没有，只能猜平均值 $\bar{Y}$，我会错多少？"|
|**SSE**|残差平方和|"用了模型以后，还剩多少错？"|
|**SSM**|模型平方和|"模型把我从瞎猜拉到多准？"|

### 1.2 几何直觉（PPT 最没讲清楚的！）

$Y$、$\hat{Y}$、$\bar{Y}$ 构成一个**直角三角形**：

```
         Y  (真实值)
        /|
       / |
      /  | SSE (残差，正交于模型空间)
     /   |
    /    |
   /_____Ŷ (预测值)
   Ȳ    SSM
```

**毕达哥拉斯定理** → $SST = SSM + SSE$。

这就是为什么投影矩阵 $P$（投到模型空间）和 $Q = I-P$（投到残差空间）**正交**：$PQ=0$。它们把 $\mathbb{R}^n$ 分成两个互相垂直的子空间。

### 1.3 Lesetest 例子中的数字

```r
fit_full <- lm(Fehler ~ GE + JG + LZ + WOL + WOG + WOTV, data = lese)
anova(fit_full)

# 假设结果:
# SST = 8500 (总变异)
# SSM = 2800 (模型解释的)
# SSE = 5700 (没解释的)
```

**R²**：  
$$
R^2 = \frac{SSM}{SST} = \frac{2800}{8500} \approx 0.33
$$

解释：模型解释了 33% 的错误数变异。剩下 67% 由别的因素或随机噪声驱动。

### 1.4 不带截距的版本（3.3）

$$
\underbrace{Y'Y}_{SST^*} = \underbrace{(Y-\hat{Y})'(Y-\hat{Y})}_{SSE} + \underbrace{\hat{Y}'\hat{Y}}_{SSM^*}
$$

区别：不扣平均值。几乎不用，除非你的模型真的没有截距（比如过原点回归）。

---

## 2. 自由度与平均平方

### 2.1 为什么要除自由度？

SSE 越大模型越差？**不对！** 样本量大了，SSE 自然大。要除以自由度，才能公平比较。

$$
MSE = \frac{SSE}{n-p'}
$$

**自由度的直觉**（和上一题 $(n-1)S^2/\sigma^2 \sim \chi^2(n-1)$ 一个道理）：

- 你用了 $p'=p+1$ 个参数去拟合数据
- 每估计一个参数，数据就少一个"自由"的维度
- 所以残差有 $n-p'$ 个自由度

### 2.2 各平均平方

|量|公式|自由度|含义|
|---|---|---|---|
|MST|$SST/(n-1)$|$n-1$|数据整体变异（= 样本方差）|
|MSE|$SSE/(n-p')$|$n-p'$|**$\sigma^2$ 的无偏估计**|
|MSM|$SSM/p$|$p$|模型单位参数的贡献|

**Lesetest 例子**：

```
MSE = 5700 / (180 - 7) = 32.95
σ̂ = √32.95 ≈ 5.74  → "模型预测通常差约 5.74 个错误"
```

---

## 3. 必备分布工具箱

### 3.1 多元正态分布

$Z \sim N_n(\mu, \Sigma)$ 的两个你必须记住的性质：

**性质 1：线性变换还是正态**  
$$
Z \sim N_n(\mu, \Sigma) \Rightarrow AZ \sim N_m(A\mu, A\Sigma A')
$$

这是为什么 $\hat{\beta} = (X'X)^{-1}X'Y$ 也是正态——它只是 $Y$ 的线性变换。

### 3.2 χ² 分布

定义：若 $Z \sim N_n(\mu, I)$，则 $Z'Z \sim \chi^2(n, \delta)$，其中 $\delta = \mu'\mu$。

**关键事实（3.21，Cochran 定理核心）**：

如果 $A$ 是秩 $r$ 的**投影矩阵**（$A^2=A, A'=A$），$Z \sim N_n(\mu, I)$，那么：

$$
Z'AZ \sim \chi^2(r, \mu'A\mu)
$$

**为什么重要**？ 因为 SSE、SSM 都是 $Y'(\text{投影矩阵})Y$ 的形式，所以它们都是 $\chi^2$ 分布的！

### 3.3 F 分布

$$
X = \frac{W_1/n_1}{W_2/n_2}, \quad W_1 \sim \chi^2(n_1, \delta), W_2 \sim \chi^2(n_2), \text{ 独立}
$$

→ $X \sim F(n_1, n_2, \delta)$

**直觉**：F = 信号/噪声 的比值。所有的 F-Test 都是这个形式。

---

## 4. KQ 估计量的分布（大结论）

正态假设下的 5 个关键结论（公式 3.24-3.29）：

$$
\hat{\beta} \sim N(\beta, \sigma^2(X'X)^{-1}) \tag{3.24}
$$

$$
(n-p')\frac{\hat{\sigma}^2}{\sigma^2} \sim \chi^2(n-p') \tag{3.26}
$$

$$
\hat{\sigma} \text{ 与 } \hat{\beta} \text{ 独立} \tag{3.27}
$$

$$
\frac{\hat{\beta}_k - \beta_k}{\hat{\sigma}_{\hat{\beta}_k}} \sim t(n-p') \tag{3.29}
$$

**Lesetest 例子里这意味着什么**？

```r
summary(fit_full)
# 假设输出:
# Coefficients:
#             Estimate  Std.Error  t value  Pr(>|t|)
# (Intercept)  15.20      2.10      7.24   <0.001 ***
# GE            1.80      0.88      2.05    0.042  *   
# JG           -3.50      0.90     -3.89   <0.001 ***
# LZ           -0.12      0.04     -3.00    0.003  **
# WOL          -1.40      0.45     -3.11    0.002  **
# WOG           0.80      0.40      2.00    0.047  *
# WOTV          0.30      0.42      0.71    0.478
```

- **Estimate** = $\hat{\beta}_k$
- **Std.Error** = $\hat{\sigma}_{\hat{\beta}_k} = \sqrt{c_{kk}} \cdot \hat{\sigma}$
- **t value** = $\hat{\beta}_k / \hat{\sigma}_{\hat{\beta}_k}$ ← 就是公式 (3.29)
- **Pr(>|t|)** = 基于 $t(173)$ 分布的双侧 p 值

**解读**：

- `JG` 的 t = −3.89 → 四年级比三年级**少错** 3.5 个（显著）
- `LZ` 的 t = −3.00 → 每多 1 分钟阅读时间，**少错** 0.12 个（显著）
- `WOTV` 的 t = 0.71 → 看电视频率**没显著影响**（p = 0.478）

---

## 5. Overall-Test：模型整体有没有用？

### 5.1 统计量

$$
F_O = \frac{MSM}{MSE} \sim F(p, n-p')
$$

### 5.2 假设

$$
H_0^O: \beta_1 = \beta_2 = \dots = \beta_p = 0
$$

"**所有斜率都是 0**" ↔ "模型就是个平均值，变量全是废的"

### 5.3 直觉

$$
F_O = \frac{\text{模型号称自己有多强}}{\text{纯噪声有多大}}
$$

- F 大 → 信号远大于噪声 → **拒绝 $H_0$**，模型有用
- F ≈ 1 → 信号和噪声差不多 → 模型没啥用

### 5.4 Lesetest 例子

```
F = MSM/MSE = (2800/6) / (5700/173) 
  = 466.7 / 32.95 
  = 14.16

p-value < 0.001  → 拒绝 H₀
```

**结论**：至少有一个变量对错误数有影响。但这**不告诉你是哪一个**——要继续做个别检验。

---

## 6. 一般线性假设 = 所有检验的"万能枪"

### 6.1 一般形式

$$
H_0: A\beta = c
$$

其中 $A$ 是 $a \times p'$ 矩阵，$c$ 是 $a$ 维向量，$\operatorname{rg}(A)=a$。

### 6.2 三个实用例子（Lesetest 版）

**例 1：WOTV 真的没用吗？**（单个系数）

$$
H_0: \beta_6 = 0
$$

$$
A = (0, 0, 0, 0, 0, 0, 1), \quad c = 0
$$

**例 2：媒体使用（WOL, WOG, WOTV）整体有用吗？**（变量组）

$$
H_0: \beta_4 = \beta_5 = \beta_6 = 0
$$

$$
A = \begin{pmatrix} 0 & 0 & 0 & 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 & 0 & 0 & 1 \end{pmatrix}, \quad c = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
$$

**例 3：打游戏和看电视的效应相同吗？**（系数相等）

$$
H_0: \beta_5 = \beta_6 \Leftrightarrow \beta_5 - \beta_6 = 0
$$

$$
A = (0, 0, 0, 0, 0, 1, -1), \quad c = 0
$$

### 6.3 Wald-Test 统计量

$$
SSH = (A\hat{\beta} - c)' [A(X'X)^{-1}A']^{-1} (A\hat{\beta} - c)
$$

$$
TF = \frac{MSH}{MSE} = \frac{SSH/a}{SSE/(n-p')} \sim F(a, n-p')
$$

### 6.4 **更实用的等价形式**（公式 3.40）

$$
TF = \frac{[SSE(H_0) - SSE] / a}{SSE / (n-p')}
$$

**翻译**：

1. 拟合**完整模型**，得 SSE
2. 拟合 **$H_0$ 下的受限模型**（强制满足约束），得 $SSE(H_0)$
3. 比较它们

**核心直觉**：

> **如果约束是合理的，$SSE(H_0)$ 不会比 $SSE$ 大多少。**  
> 如果差很多，说明约束是错的 → 拒绝 $H_0$。

### 6.5 Lesetest 例子：媒体使用变量组

```r
# 完整模型
fit_full <- lm(Fehler ~ GE + JG + LZ + WOL + WOG + WOTV, data = lese)
SSE_full <- sum(resid(fit_full)^2)  # 假设 = 5700, df = 173

# 受限模型（去掉 3 个媒体变量）
fit_restricted <- lm(Fehler ~ GE + JG + LZ, data = lese)
SSE_restricted <- sum(resid(fit_restricted)^2)  # 假设 = 6300, df = 176

# F 统计量
F_stat <- ((SSE_restricted - SSE_full) / 3) / (SSE_full / 173)
# = (600/3) / (5700/173) = 200 / 32.95 = 6.07

# p 值
pf(F_stat, df1 = 3, df2 = 173, lower.tail = FALSE)
# ≈ 0.0006

# 或者直接用 R 的 anova:
anova(fit_restricted, fit_full)
```

**结论**：$F = 6.07$，$p < 0.001$ → 媒体使用变量组**整体显著**，即使 WOTV 单独不显著。

---

## 7. 所有检验本质上是一个！

PPT 把 Overall-Test、t-Test、Wald-Test、LQ-Test 当成不同东西讲，让人迷惑。其实：

|检验|约束矩阵 A|本质|
|---|---|---|
|**双侧 t-Test**（$\beta_k = \beta_k^0$）|$A = (0,\dots,1,\dots,0)$|Wald-Test 特例|
|**Overall-Test**|$A$ = 选所有斜率|Wald-Test 特例|
|**变量组 F-Test**|$A$ = 选一组|就是 Wald-Test|
|**Likelihood-Quotienten-Test**|任意约束|数学上等价于 Wald-Test（公式 3.45）|

**LQ 与 Wald 的精确关系**：

$$
\tau_{LQ} = \left[\tau_W \cdot \frac{a}{n-p'} + 1\right]^{-n/2}
$$

g 严格单调 → 两种检验产生**同样的拒绝域、同样的 p 值**。

### 7.1 Overall-Test 是 Wald-Test 的特例

若 $H_0: \beta_1 = \dots = \beta_p = 0$，则受限模型 $Y = \beta_0 + \varepsilon$ 的最优估计是 $\hat{\hat{\beta}}_0 = \bar{Y}$：

$$
SSE(H_0) = \sum(Y_i - \bar{Y})^2 = SST
$$

$$
SSH = SSE(H_0) - SSE = SST - SSE = SSM
$$

所以 Overall-Test 确实只是特殊情况。

---

## 8. Sequential vs. Partial 平方和（理解 R 输出的关键）

### 8.1 问题

看 Lesetest 模型，`anova(fit)` 和 `summary(fit)` 给的 p 值**不一样**！ 为什么？

### 8.2 两种计算方式

**Sequential SS (Type I)**：按顺序加变量，每个变量的贡献 = 加入时 SSE 的下降。

```
M0: Y = β0                             SSE = 8500
M1: + GE                               SSE = 8200   → GE 贡献 300
M2: + JG                               SSE = 7500   → JG 贡献 700
M3: + LZ                               SSE = 6800   → LZ 贡献 700
M4: + WOL                              SSE = 6200   → WOL 贡献 600
M5: + WOG                              SSE = 5900   → WOG 贡献 300
M6: + WOTV                             SSE = 5700   → WOTV 贡献 200
```

加起来 = $2800 = SSM$ ✓ （完美分解 SST）

**但问题**：如果把 JG 和 LZ 调换顺序，分出来的"贡献"会变！

**Partial SS (Type III)**：每个变量的贡献 = 假设其他变量**全部在**，再加这个变量带来的 SSE 下降。

```
R(β_WOTV | 其他所有变量) = SSE(无 WOTV) - SSE(完整) = 小
R(β_JG | 其他所有变量)    = SSE(无 JG)   - SSE(完整) = 大
```

### 8.3 **关键区别**

- **Sequential SS** 相加 = $SSM$ ✓
- **Partial SS** 相加 ≠ $SSM$ ❌（公式 3.48）

当变量之间**相关**时，两者会不同。

### 8.4 类比

- **Sequential**：公司按入职顺序给功劳
- **Partial**：假设另外两个人已经在了，评估你**边际**贡献

### 8.5 在 R 里

```r
anova(fit_full)        # Type I (Sequential) - 顺序敏感！
car::Anova(fit_full)   # Type III (Partial) - summary() 里的 p 值
```

**实务建议**：做假设检验看 **Type III (Partial)**；做方差分解（比如 ANOVA）看 Type I。

---

## 9. 变量的"净效应"：Bereinigung（公式 3.54）

### 9.1 核心结论（Frisch-Waugh-Lovell 定理）

把回归 $Y = X_1\beta_1 + X_2\beta_2 + \varepsilon$ 拆成两步：

1. 把 $Y$ 和 $X_1$ **都对 $X_2$ 回归**，取残差：$Y^*=Q_2Y$，$X_1^*=Q_2X_1$
2. 把 $Y^*$ 对 $X_1^*$ 回归 → **得到的系数就是原模型的 $\hat{\beta}_1$**！

$$
\hat{\beta}_1 = (X_1^{*\prime} X_1^*)^{-1} X_1^{*\prime} Y^*
$$

### 9.2 翻译成人话

**"多元回归中 $\beta_1$ 的意思"** = "把 $X_1$ 和 $Y$ 中被其他变量能解释的部分**都扣掉**，剩下的纯净关系"。

### 9.3 Lesetest 例子

`LZ` 的系数 $\hat{\beta}_3 = -0.12$ 怎么解释？

**错误**：每多阅读 1 分钟，错误减少 0.12。

**正确**：**在性别、年级、课外阅读、游戏、电视都保持不变的前提下**，每多阅读 1 分钟，错误减少 0.12。

怎么可视化？**Partial Leverage Plot**：

1. 把 `LZ` 对其他 5 个变量回归 → 取残差 `LZ*`
2. 把 `Fehler` 对其他 5 个变量回归 → 取残差 `Fehler*`
3. 画 `Fehler*` vs `LZ*` 的散点图 + 回归线

这条线的斜率就是 $\hat{\beta}_3 = -0.12$。

```r
library(car)
avPlots(fit_full)   # 一次画所有 partial leverage plots
```

---

## 10. 预测：置信区间 vs. 预测区间

### 10.1 场景

一个新学生：男生、4 年级、阅读时间 30 分钟、WOL=2、WOG=1、WOTV=3。  
他的错误数是多少？

**点预测**：  
$$
\hat{Y}_{n+1} = x'_{n+1}\hat{\beta}
$$

假设算出 $\hat{Y}_{n+1} = 10$。

### 10.2 两种区间（常被混淆！）

**置信区间**（for $\mu_{n+1}$）：  
$$
\hat{Y}_{n+1} \pm t_{1-\alpha/2}(n-p') \cdot \hat{\sigma} \sqrt{x'_{n+1}(X'X)^{-1} x_{n+1}}
$$

**回答**："**这类学生的平均错误数**在哪个范围？"

**预测区间**（for $Y_{n+1}$）：  
$$
\hat{Y}_{n+1} \pm t_{1-\alpha/2}(n-p') \cdot \hat{\sigma} \sqrt{1 + x'_{n+1}(X'X)^{-1} x_{n+1}}
$$

**回答**："**这个具体学生**的错误数在哪个范围？"

### 10.3 为什么预测区间更宽？

预测区间多了那个 **+1**：不仅要考虑 $\hat{\beta}$ 的不确定性，还要加上**新个体自己的噪声** $\varepsilon_{n+1}$。

**Lesetest 数字**：

```
95% CI for μ:  [8.5, 11.5]    宽度 3     "这类学生平均错 8.5-11.5 个"
95% PI for Y:  [1.2, 18.8]    宽度 17.6  "他本人会错 1.2-18.8 个"
```

个体波动远大于平均值波动，所以 PI 宽得多。

### 10.4 实务建议

- **展示模型时** → 画 CI（体现"你的估计有多准"）
- **对个体做预测时** → 报 PI（体现"实际会多不准"）

```r
predict(fit_full, newdata, interval = "confidence")  # CI
predict(fit_full, newdata, interval = "prediction")  # PI
```

---

## 11. Gauss-Markov 定理：KQ 估计量好在哪？

### 11.1 陈述

在假设 $E(\varepsilon)=0$、$V(\varepsilon)=\sigma^2 I$ 下（**不需要正态性**！）：

> KQ 估计量 $\hat{\beta}$ 是**BLUE**（Best Linear Unbiased Estimator）——在所有无偏线性估计量中方差最小。

### 11.2 翻译

- **L**inear：$\hat{\beta}=CY$ 形式
- **U**nbiased：$E(\hat{\beta})=\beta$
- **B**est：方差最小（在所有 LU 估计量中）

### 11.3 注意事项

- 如果你愿意接受**有偏**估计，可能方差更小（例如 Ridge Regression）
- 如果误差**不是同方差**，KQ 就不 BLUE 了，要用 WLS 或 GLS

---

## 12. 大样本渐近性质

### 12.1 一致性（公式 3.60）

$$
\lim_{n \to \infty}(X_n' X_n)^{-1} = 0
$$

$$
\hat{\beta}^{(n)} \xrightarrow{P} \beta
$$

**人话**：样本量越大，估计越接近真值。

### 12.2 渐近正态性（公式 3.62）

$$
\hat{\beta}^{(n)} \stackrel{a}{\sim} N(\beta, \hat{\sigma}_n^2(X_n' X_n)^{-1})
$$

**实用意义**：即使误差**不是正态**的，只要 n 够大，所有 t 检验、F 检验的结果**近似**仍然成立。

---

## 🎁 总结：一张图看完整章

```
┌──────────────────────────────────────────────────────────┐
│                                                          │
│  数据 Y ──[拟合模型]──> 预测 Ŷ                           │
│                                                          │
│  SST = SSM + SSE   (毕达哥拉斯：直角三角形)              │
│   ↓      ↓     ↓                                         │
│  总变异 模型贡献 残差                                    │
│                                                          │
│  自由度:  n-1    p    n-p'                               │
│                                                          │
│  ┌─────────────────────────────────────┐                 │
│  │ 所有假设检验的统一框架:             │                 │
│  │                                     │                 │
│  │   F = [SSE(H₀) - SSE] / a           │                 │
│  │       ─────────────────────         │                 │
│  │       SSE / (n-p')                  │                 │
│  │                                     │                 │
│  │   "约束后 SSE 变大很多吗?"          │                 │
│  └─────────────────────────────────────┘                 │
│                                                          │
│  特例:                                                   │
│  • t-Test      → 1 个系数                                │
│  • Overall-Test → 所有斜率 = 0                           │
│  • 变量组检验   → 一组系数 = 0                           │
│  • LQ-Test     → 数学上等价                              │
│                                                          │
│  预测:                                                   │
│  • CI (均值):   ±t·σ̂·√(x'(X'X)⁻¹x)                       │
│  • PI (个体):   ±t·σ̂·√(1+x'(X'X)⁻¹x)   (宽得多)          │
│                                                          │
│  估计量的好性质:                                         │
│  • BLUE (Gauss-Markov, 有限样本)                         │
│  • 一致 + 渐近正态 (大样本)                              │
│                                                          │
└──────────────────────────────────────────────────────────┘
```

---

## 📝 一页小抄

|我想做...|用什么|R 代码|
|---|---|---|
|整个模型有没有用？|Overall-Test (F)|`summary(fit)` 最后一行|
|单个 $\beta_k$ 是 0 吗？|t-Test|`summary(fit)`|
|一组 $\beta$ 都是 0 吗？|F-Test / Wald|`anova(fit_0, fit_1)`|
|两个系数相等吗？|线性假设|`car::linearHypothesis()`|
|$\beta_k$ 的"净效应"什么样？|Partial Leverage|`car::avPlots(fit)`|
|预测新观测的均值|置信区间|`predict(..., "confidence")`|
|预测新个体的实际值|预测区间|`predict(..., "prediction")`|
|$\sigma^2$ 的估计|MSE|`summary(fit)$sigma^2`|
|$R^2$|$SSM/SST$|`summary(fit)$r.squared`|

---

**有什么想深入的？** 告诉我你最卡的是哪一块（几何直觉？Cochran 定理证明？partial vs. sequential？），我可以再展开。
