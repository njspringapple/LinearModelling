# Linear Modelling PPT 与三本教材章节映射笔记

依据：根目录下的 `LinearModellingPPT.pdf`、三本教材 PDF，以及已经拆分出的 `chapter_pdfs/LinearModellingPPT/` 目录。

## 1. 哪本是主教材

主教材应当是 Fahrmeir, Kneib, Lang, Marx 的 `Regression: Models, Methods and Applications`，也就是根目录下的 `Regression.pdf`。

判断理由：

- PPT 的 `Literatur` 页把 Fahrmeir et al. 的 `Regression` 放在第一项，并注明可在线获取；两本 Faraway 书列在后面。
- PPT 正文大量引用 `vgl. Fahrmeir et al.; S. ...`。按文本检索，`Fahrmeir` 在 PPT 中出现 40 次，`Faraway` 只出现在参考书列表中 1 次。
- PPT 的课程顺序与 `Regression.pdf` 的目录非常贴合：先是回归导论与经典线性模型，再到诊断、广义/加权最小二乘、混合模型、logistic 回归、变量选择和非参数/样条扩展。

所以：`Regression.pdf` 是课程主线；`Linear Models With R` 和 `Extending The Linear Model With R` 更像 R 实操和扩展阅读。

## 2. 三本教材的区别

### `Regression.pdf`

角色：理论主教材。

特点：覆盖面最像本课讲义，数学推导更集中，章节从线性回归的矩阵表达、假设检验、模型选择，延伸到 GLS、正则化、GLM、混合模型、非参数回归和附录中的矩阵/概率推断基础。PPT 中的页码引用基本都指向这本书。

适合用途：复习证明、理解公式来源、对照 PPT 的主线内容。

### `Linear Models With R (Julian James Faraway).pdf`

角色：经典线性模型的 R 实操辅助书。

特点：聚焦普通线性模型、诊断、误差/预测变量问题、变量选择、ANOVA、ANCOVA、实验设计，并配合 R 做数据分析。它对 PPT 前半部分和变量选择很有帮助，但不系统覆盖 logistic 回归和混合模型。

适合用途：用 R 跑线性回归、诊断图、ANOVA/ANCOVA、变量选择；把 PPT 里的抽象公式落到代码和例子上。

### `Extending The Linear Model With R (Second Edition) (Julian James Faraway).pdf`

角色：线性模型之后的扩展模型 R 辅助书。

特点：可以看作 Faraway `Linear Models with R` 的后续：重点是 GLM、binary/binomial/count response、logistic 变体、随机效应、重复测量/纵向数据、GLMM、非参数回归、additive models、trees 和 neural networks。

适合用途：补 PPT 后半部分的 logistic 回归、mixed model、样条/非参数模型；也适合后续 `Generalisierte Regression` 类课程的预习。

## 3. PPT 主题到教材章节的映射

### 01 `Einführung und Beispiele`

教学内容总结：课程导入，解释 regression 的概念和应用场景；用身高回归、环境政策、教育、医学/流行病学等例子说明为什么要建模变量之间的关系，并给出课程整体路线。

主教材 `Regression.pdf`：Chapter 1 `Introduction`；Chapter 2 `Regression Models`，尤其 2.1-2.3 的模型概览。

`Linear Models With R`：Chapter 1 `Introduction`，尤其初步数据分析和何时使用回归。

`Extending The Linear Model With R`：Chapter 1 `Introduction` 可作为扩展模型的预告，不是本主题的主读物。

### 02 `Das einfache lineare Regressionsmodell`

教学内容总结：简单线性回归模型、最小二乘思想、参数解释、估计量性质、残差、误差方差估计、置信区间、平方和分解、`R^2`、预测与线性近似。

主教材 `Regression.pdf`：Chapter 2.2.1 `Simple Linear Regression Model`；Chapter 3.1 `Model Definition`；Chapter 3.2 `Parameter Estimation`；Chapter 3.3.2 `Confidence Regions and Prediction Intervals`。

`Linear Models With R`：Chapter 2 `Estimation`；Chapter 3 `Inference` 的预测和置信区间部分。

`Extending The Linear Model With R`：Chapter 1 是复习性质；不是本主题的主要教材。

### 03 `Das multiple lineare Regressionsmodell`

教学内容总结：多元线性模型的矩阵表达 `Y = X beta + epsilon`、设计矩阵、模型假设、最小二乘估计、hat matrix、residual matrix、投影几何、残差性质和误差方差估计。

主教材 `Regression.pdf`：Chapter 2.2.2 `Multiple Linear Regression`；Chapter 3.1 `Model Definition`；Chapter 3.2 `Parameter Estimation`；Appendix A `Matrix Algebra`。

`Linear Models With R`：Chapter 2 `Estimation`，尤其矩阵表示、最小二乘、Gauss-Markov theorem 和 goodness of fit。

`Extending The Linear Model With R`：Chapter 1 只作复习。

### 04 `Quadratsummenzerlegung und statistische Inferenz im multiplen linearen Regressionsmodell`

教学内容总结：平方和分解、带/不带截距的分解、均方、正态假设下的估计量分布、卡方/t/F 分布、Cochran 定理、overall test、一般线性假设、Wald test、似然比检验、受限模型、顺序/部分平方和、partial leverage plot、预测区间、Gauss-Markov、一致性和渐近正态性。

主教材 `Regression.pdf`：Chapter 3.3 `Hypothesis Testing and Confidence Intervals`；Chapter 3.4 `Model Choice and Variable Selection`；Appendix B `Probability Calculus and Statistical Inference`；必要时配合 Appendix A 的二次型和矩阵结果。

`Linear Models With R`：Chapter 3 `Inference`；Chapter 8 `Variable Selection` 中的模型比较思想。

`Extending The Linear Model With R`：Appendix A `Likelihood Theory` 对 ML/LRT 的扩展背景有帮助，但线性模型推断主线仍看 `Regression.pdf`。

### 05 `Diskrete Einflussgrößen: Dummy- und Effektkodierung, Mehrfaktorielle Varianzanalyse`

教学内容总结：离散解释变量、dummy coding、effect coding、组均值模型、单因素/多因素方差分析、因素交互、离散变量与连续变量结合的协方差分析、不同类型平方和、前后测比较。

主教材 `Regression.pdf`：Chapter 3.1.3 `Modeling the Effects of Covariates`；Chapter 3.3 的线性假设检验；Chapter 3.4 的模型比较。

`Linear Models With R`：Chapter 13 `Analysis of Covariance`；Chapter 14 `One-Way Analysis of Variance`；Chapter 15 `Factorial Designs`；Chapter 16 `Block Designs`。这部分 Faraway 反而非常好用，因为它把分类解释变量和 ANOVA/ANCOVA 展开得更实操。

`Extending The Linear Model With R`：不是本主题的主要读物；若把分组/区组作为随机效应，可预读 Chapter 10 `Random Effects`。

### 06 `Metrische Einflussgrößen: Interaktionen, Polynome, Splines, Transformationen`

教学内容总结：连续解释变量的交互、非线性趋势的线性模型表示、多项式回归、三角多项式、回归样条、基函数思想、变量变换、多重检验和同时置信区间。

主教材 `Regression.pdf`：Chapter 3.1.3 `Modeling the Effects of Covariates`；Chapter 8.1 `Univariate Smoothing`，尤其 8.1.1 `Polynomial Splines` 和 8.1.2 `Penalized Splines`；Chapter 9.3 `Models with Interactions`。

`Linear Models With R`：Chapter 7 `Transformation`；Chapter 13 `Analysis of Covariance` 可辅助理解连续变量与分组变量的交互。

`Extending The Linear Model With R`：Chapter 14 `Nonparametric Regression`；Chapter 15 `Additive Models`，适合补样条和更灵活的非参数建模。

### 07 `Modelldiagnose`

教学内容总结：残差类型、标准化/学生化残差、交叉验证残差、递归残差、正态性问题、异方差、自相关和 Durbin-Watson test、离群点、高杠杆点、Cook's distance、模型设定错误、partial leverage plot、共线性、VIF/条件数、测量误差。

主教材 `Regression.pdf`：Chapter 3.4.4 `Model Diagnosis`；Chapter 4.1.3 `Heteroscedastic Errors`；Chapter 4.1.4 `Autocorrelated Errors`；Chapter 4.2 `Regularization Techniques` 可作为共线性/正则化补充。

`Linear Models With R`：Chapter 4 `Diagnostics`；Chapter 5 `Problems with the Predictors`；Chapter 6 `Problems with the Error`。这本书是本主题的强辅助材料。

`Extending The Linear Model With R`：Chapter 8.4 `GLM Diagnostics` 适合之后类比到 GLM；普通线性模型诊断仍优先看前两本。

### 08 `Das allgemeine lineare Modell: Gewichtete KQ-Methode, autokorrelierte und heteroskedastische Störterme`

教学内容总结：一般线性模型、非独立同方差误差结构、模型变换、weighted least squares、generalized least squares、广义 Gauss-Markov theorem、方差结构、AR(1) 时间序列误差、方差参数估计、ML/REML 和相关推断。

主教材 `Regression.pdf`：Chapter 4.1 `The General Linear Model`，尤其 4.1.2 `Weighted Least Squares`、4.1.3 `Heteroscedastic Errors`、4.1.4 `Autocorrelated Errors`。

`Linear Models With R`：Chapter 6.1 `Generalized Least Squares`；Chapter 6.2 `Weighted Least Squares`。

`Extending The Linear Model With R`：没有完全对应的普通线性模型章节；GLM 的 sandwich/robust ideas 可参考 Chapter 8.5-8.6，但这不是本主题主线。

### 09 `Das gemischte lineare Regressionsmodell (Linear mixed Model)`

教学内容总结：随机截距/方差分量模型、边际模型、纵向数据和重复测量、固定效应与随机效应、混合模型作为一般线性模型的表示、模型推断和结果解释。

主教材 `Regression.pdf`：Chapter 7 `Mixed Models`，尤其 7.1 `Linear Mixed Models for Longitudinal and Clustered Data`、7.3 `Likelihood Inference in LMMs`、7.7 `Practical Application of Mixed Models`。

`Linear Models With R`：没有系统的 mixed model 章节；Chapter 16 `Block Designs` 与区组/随机效应思想相邻，但不足以覆盖本主题。

`Extending The Linear Model With R`：Chapter 10 `Random Effects`；Chapter 11 `Repeated Measures and Longitudinal Data`；如果进入非正态响应的混合模型，再看 Chapter 13 `Mixed Effect Models for Nonnormal Responses`。

### 10 `Das logistische Regressionsmodell`

教学内容总结：二元响应变量、logistic regression、logit/odds/log-odds 解释、分类问题、最大似然估计、MLE 渐近性质、存在唯一性问题、Wald test、likelihood-ratio test、拟合优度、ROC/AUC、病例-对照研究和 propensity score 相关计算。

主教材 `Regression.pdf`：Chapter 2.3 `Regression with Binary Response Variables: The Logit Model`；Chapter 5.1 `Binary Regression`；Appendix B.4 `Likelihood Inference`。

`Linear Models With R`：基本不覆盖本主题，只能提供线性模型背景。

`Extending The Linear Model With R`：Chapter 2 `Binary Response`；Chapter 4 `Variations on Logistic Regression`；Chapter 8 `Generalized Linear Models` 可补 GLM 总框架。

### 11 `Variablenselektion`

教学内容总结：建模目标、解释型建模与预测型建模、过拟合和偏差-方差权衡、模型优度指标、`R^2`/调整 `R^2`、AIC/BIC、交叉验证、前向/后向/逐步选择、变量选择后的推断风险、模型不确定性，以及 LASSO/boosting 等后续方向。

主教材 `Regression.pdf`：Chapter 3.4 `Model Choice and Variable Selection`；Chapter 4.2 `Regularization Techniques`，尤其 ridge 和 LASSO；Chapter 4.3 `Boosting Linear Regression Models` 可作为拓展。

`Linear Models With R`：Chapter 8 `Variable Selection`；Chapter 9 `Shrinkage Methods`；Chapter 10 `Statistical Strategy and Model Uncertainty`。

`Extending The Linear Model With R`：Chapter 2.5 `Model Selection` 可看 logistic 场景；Appendix A.3 `Model Selection` 可看 likelihood 视角；更高级的树和神经网络在 Chapter 16-17，但已经超出本讲义主线。

## 4. 学习顺序建议

1. 先按 PPT 主题读 `Regression.pdf` 的对应章节，因为它是课程主线。
2. 遇到需要上机或 R 代码的普通线性模型、诊断、ANOVA、变量选择，再补 `Linear Models With R`。
3. 遇到 logistic regression、mixed model、splines/additive models，再补 `Extending The Linear Model With R`。
4. PPT 还列了 Thomas Nagler 的线性代数讲义；它不是三本教材之一，但对 hat matrix、projection、quadratic forms 等证明很有用。
