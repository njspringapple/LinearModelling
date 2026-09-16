data(mtcars)
# r语言的公式语法，右侧是把一些参数加入自变量，不是数学的加法
# 意思是 mpg = beta0 + beta1*wt + beta2*hp + beta3*cyl
fit <- lm(mpg ~ wt + hp + cyl, data = mtcars)
fit

# 提取模型的公式，给人看
# y = Xbeta + e
X <- model.matrix(fit)
y <- mtcars[["mpg"]]

dim(X)
head(X)
head(y)

summary(fit) 

#Call:
#  lm(formula = mpg ~ wt + hp + cyl, data = mtcars)
#Residuals:
#  Min      1Q  Median      3Q     Max 
#-3.9290 -1.5598 -0.5311  1.1850  5.8986 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 38.75179    1.78686  21.687  < 2e-16 ***
#  wt          -3.16697    0.74058  -4.276 0.000199 ***
#  hp          -0.01804    0.01188  -1.519 0.140015    
#  cyl         -0.94162    0.55092  -1.709 0.098480 .  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 2.512 on 28 degrees of freedom
#Multiple R-squared:  0.8431,	Adjusted R-squared:  0.8263 
#F-statistic: 50.17 on 3 and 28 DF,  p-value: 2.184e-11
