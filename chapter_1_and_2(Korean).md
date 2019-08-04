## Preface

- 책은 총 5개의 theme으로 구성됨
  1. introduces the reader to basic concepts and terminology (1~3)
  2. the estimation of summary survival statistics (4~6)
  3. hypothesis testing (7)
  4. regression analysis for censored and truncated data (★, 8~12)
  5. multivariate models for survival data (13)



## 1. Examples of Survival Data

#### 1.2 Remission Duration from a Clinical Trial for Acute Leukemia

- n = 42 (급성 백혈병 판정 후 완전/부분 회복중인 환자들)
- 완전/부분 회복 여부에 따라 2인 1개조로 편성 후 투약(6-MP) / 위약으로 구분함
- 재발까지 걸린 시간(월단뒤)을 측정하였고, 일부 left-censoring이 존재함
- 사용 예제
  - chapter 4 - product-limit estimator / Nelson-Aalen estimator / mean survival time and standard error
  - chapter 6.4 - estimating survival function using Bayesian approaches
  - chapter 7.5 - stratified log rank test
  - chapter 9.3 - stratified proportional hazards model

#### 1.3 Bone Marrow Transplantation for Leukemia

- n = 137 (99 AML, 38 ALL) - 급성 백혈병 종류임
- 42명 재발, 41명 회복 중 사망, 26명 급성 GVHD 발병, 17명 혈소판 회복과정 없이 바로 재발 또는 사망
- 몰라 존나 복잡함...  Appendix D 참조
- 사용 예제
  - chapter 4 - product-limit estimator / Nelson-Aalen estimator / pointwise confidence intervals and confidence bands / competing risks (relapse and death)
  - chapter 6.2 - construction of estimates of the hazard rate
  - chapter 7 - tests for the equality of $K$ survival curves (stratified / unstratified)
  - chapter 8 - tests for the equality of $K$ hazard rates adjusted for possible fixed-time confounders
  - chapter 9 - include covariates on the model in chapter 8, whose values change over time
  - chapter 11 - regression diagnostics

#### 1.4 Times to Infection of Kedney Dialysis Patients

- n = 119 (신부전증으로 인한 투석 치료 환자)
- 43 - 수술로 카테터 삽관 / 76 - 경피 카테터 삽관
- 삽관부위 감염증상이 나타날 때까지 걸린 시간을 측정 (월단위), 일부 censoring 존재
- 사용 예제
  - chapter 7 - inference about the equality of two survival curves (two-sample weighted, log-rank test) / Cramor-von Mises test for censored data
  - chapter 8 - methods for constructing partial likelihoods and the subsequent testing of equality of the survival curve
  - chapter 9 - testing for proportional hazards (proportional hazards assumption for this data is not correct)

#### 1.5 Times fo Death for a Breast-Cancer Trial

- n=45 (겨드랑이 부위 림프절에서 SLM 방법에 의해 음성 판정을 받은 환자들)
- immunohistochemical(IH) 방법으로 재검사한 결과 9명은 양성, 36명은 음성 판정을 받음
- 각 그룹에 대해 생존 시간을 측정 (월단위), 일부 censoring 존재
- 사용 예제
  - chapter 8 - construction of partial likelihood functions / two-sample test based on proprotional hazards with no ties with right-censored data
  - chapter 10 - least-squares estimation in additive hazards model

#### 1.6 Times to Infection for Burn Patients



## 2. Basic Quantities and Models

#### 2.1 Introduction

**target random variable**

- $X$ : time until some specified event
  - nonnegative random variable from a homogeneous population

**Characteristic functions of** $X$

- $S(x)$ (*survival function*) : probability of an individual surviving to time $x$
- $\mathcal{h}(x)$ (*hazard rate, risk*) : chance an individual of age $x$ experiences the event in the next instant in time
- $f(x)$ (*pdf*) : unconditional probability of the event's occuring at time $x$
- $\text{mrl}(x)$ (*mean residual life*) : the mean time to the event of interest, given the event has not occurred a $x$
- If we know any one of these four, then the other three can be "uniquely" determined.

**Competing risk context**

- $h_i(t)$ (*cause-specific hazard rate*) : the rate at which subjects who have yet to experience any of the competing risks are experiencing the $i$th competing cause of failure



#### 2.2 The Survival Function

- $S(x) = Pr(X>x)$ : $x$ 시점까지 이벤트가 발생하지 않을(생존할) 확률
  - $\displaystyle S(x) = \int_x^{\infty}{f(t)dt}$ ($X$가 연속형일 때)
  - $\displaystyle S(x) = \sum_{x_j>x}p(x_j)$ ($X$가 이산형일 때)
  - 조건: monotonic decreasing function (x=0일 때 1이고, x가 무한일 때 0)



#### 2.3 The Hazard Function

**definitions**

- **hazard function** : $\displaystyle h(x) = \lim_{\Delta x \rightarrow 0}{\cfrac{P[x\le X < x+\Delta x | X \ge x]}{\Delta x}}$ 

  - $x$시점까지 생존했을 때 instant하게 이벤트가 발생할(죽을) 조건부 가능성

  - $\displaystyle h(x) = f(x)/S(x) = -d\ln[S(x)]/dx$
  - 유일한 조건은 모든 $x$에 대해 non-negative일 것

- **cumulative hazard function** : $\displaystyle H(x) = \int_0^x{h(u)du} = -\ln[S(x)]$ 

  - $S(x) = \exp[-H(x)]$

**hazard function의 종류들**

- increasing hazard : 시간이 지날수록 hazard가 (monotonic) 증가
  - natural aging or wear를 잘 설명함
- decreasing hazard : 시간이 지날수록 hazard가 (monotonic) 감소
  - 일반적이지는 않지만 특정 가전제품 수명 등을 설명하는데에 사용
- bathtub-shaped hazard : 초기와 후기의 hazard가 높고 그 중간이 낮음
  - 가장 많이 사용, 사람의 생애 hazard를 표현
- constant hazard : 시간 변화에 복립적인 hazard
  - 공산품의 수명을 설명할 수 있음
- hump-shaped hazard : 초기 hazard가 높고 종국엔 hazard가 0으로 수렴
  - 수술 후 경과나 감염증 치료 후 생존률 모델링에 사용

#### Notes

- hazard function이 event 발생에 대한 정보를 더 많이 담고있기 때문에 주된 관심의 대상이 됨

- 연속형 확률변수 $X$의 경우

  - **IFR** (increasing failure-rate) property:

    $h(x) : \text{monotonic increasing }\forall x\ge0$

  - **IFRA** (increasing failure-rate on the average):

    $H(x)/x\text{ : monotonic increasing }\forall x >0$

  - **DFR** (decreasing failure-rate) property:

    $h(x)\text{ : monotonic decreasing }\forall x \ge 0$
    
  - 위 내용들은 stochastic process의 "Reliability theory"에 관한 내용임



#### 2.4 The Mean Residual Life Function and Median Life

**definition**

- **mean residual life** : $\text{mrl}(x)=E(X-x|X>x)$

  - measurement of expected "remaining" lifetime given age $x$

  - $\text{mrl}(x) = \cfrac{\int_x^{\infty}{(t-x)f(t)dt}}{S(x)}=\cfrac{\int_x^{\infty}{S(t)dt}}{S(x)}$ (연속형인 경우)

**Relationships**

- $\displaystyle E(X) =\text{mrl}(0)= \int _0^{\infty}tf(t)dt = \int_0^{\infty}S(t)dt$
- $\displaystyle Var(X)=2\int_0^{\infty}tS(t)dt - \left[ \int_0^{\infty}S(t)st \right]^2$
- $x_p\ (p\text{th quantile of }X) = \inf\left\{ t:S(t)\le 1-p \right\}$
  - can be found by solving $S(x_p)=1-p$

**Notes**

- **median residual life** : $\text{mdrl}(x) = \inf\left\{ t:S(t)\le 0.5 \ |\ t\ge x \right\}$



#### 2.5 Common Parametric Models for survival Data

(density, survival, hazard function은 교재 37페이지 참고)

(여기엔 특징만 간단히 언급)

- **Exponential Distribution**
  - hazard function이 상수함수임 (constant hazard rate)
  - "no-aging" property(memoryless) : $P(X\ge x+z\ |\ X\ge x) = P(X\ge z)$
  - constant hazard rate appears too restrictive in both health and industrial applications
  - special case of both the **Weibull & gamma** distribution
- **Weibull Distribution**
  - describing the life length of materials
  - fully explained with scale parameter($\lambda$) and shape parameter($\alpha$)
  - shape parameter $\alpha$ can describe increasing($\alpha>1$), decreasing($\alpha<1$), constant hazard($\alpha=1$) rate
  - popular because of simple form of hazard rate and rich explanation
  - useful to work with the logarithm of the lifetime ($Y=\log X$)
    - general linear model format : $Y=\mu+\sigma E$
    - $\mu = (-\ln \lambda)/\alpha$
    - $\sigma={\alpha}^{-1}$
    - $E \sim \exp(w-e^w)$ : standard extreme value distribution
      ($-\infty<w<\infty$)
- **Log Normal Distribution**
  - $Y=\ln X \sim \text{normal distribution}$
  - popularizaed because of its relationship to the normal distribution
  - fully explained with scale($\mu $) and shape($\sigma $) parameter
  - $\displaystyle S(x) = 1-\Phi\left[ \cfrac{\ln x - \mu}{\sigma} \right]$ : can be expressed with standard normal
  - **hump-shaped hazard rate**
  - may fit cases where large values of $x$ are not of interest
- **Log Logistic Distribution**
  - $Y = \ln X \sim \text{logistic distribution}$
  - survival function is mathematically more tractable than normal
  - general linear model format : $Y = \mu + \sigma W$
    - $W$ : standardized logistic distribution with $\mu = 0,\ \sigma =1$
  - hazard function $h(x) = \cfrac{\alpha\lambda x^{\alpha -1}}{1+\lambda x^{\alpha}}$
    - numerator : Weibull hazard rate
    - denominator : make $h(x)$ as monotonic function
  - hazard is similar with log normal, but has heavy tail than that (?)
- **Gamma distribution**
  - similar to the Weibull distribution, but not mathematically tractable
  - fully explained with scale ($\lambda $) and shape ($\beta$) parameters
  - $\beta$: make $h(x)$ monotonic function

**Notes**

- (exponential) plotting $H(x)$ vs $x$ : empirical check for an exponential fit
  result should be a straignt line through the origin, with slope $\lambda$
- (Weibull) plotting $\ln\left[ H(x) \right]$ vs $\ln (x)$ : empirical chech for an Weibull fit
  result should be a straigt line with slphe $\alpha$ and intercept $\ln (\lambda)$
- $\phi$ : guarantee time parameter (최소 $\phi$ 동안은 이벤트가 발생하지 않음)
  - ex) $S(x) = \exp\left[ -\lambda (x-\phi )^{\alpha} \right]$ : Weibull guarantee model



#### 2.6 Regression Models for Survival Data

idea: homogeneous population $\rightarrow $ adjusting the survival function with covariates

**General Ideas**

- $X$ : failure time (as discussed before)
- $\mathbf{Z}^t = (Z_1,\cdots,Z_p)$ : explanatory variables associated with $X$
  - qualitative or quantitative
  - time dependent $\mathbf{Z}^t(x) = \left[ Z_1(x),\cdots,Z_p(x) \right]$
  - time dependent covariates can express **intermediate event** has occured
- **Main interest: ascertain the relationship** b/w $X$ and $\mathbf{Z}$

**Ⅰ.  accelerated failure-time model (AFT) approach (classical linear regression)**

- target: $Y = \ln(X)$
- model: $Y = \mu + \boldsymbol{\gamma}^t\mathbf{Z}+\sigma W$
  - $\boldsymbol{\gamma}$ : vector of regression coefficients
  - $W$ : error distribution
  - (log normal) $W$ : standard normal distribution
  - (Weibull) $W$ : extreme value distribution ($\exp(w-e^w)$)
  - (log logistic) $W$ : logistic distribution
- estimation: via ML methods (in chapter 12)
- description
  - $\begin{align} Pr\left[ X > x |\mathbf{Z} \right] & = Pr\left[ Y > \ln x |\mathbf{Z} \right]\\&= Pr\left[ \mu + \sigma W > \ln x - \boldsymbol{\gamma}^t\mathbf{Z} | \mathbf{Z} \right] \\ &= Pr \left[ e^{\mu + \sigma W }>x \exp{(-\boldsymbol{\gamma}^t\mathbf{Z})} \right] \\& =S_0\left[ x\exp{(-\boldsymbol{\gamma}^t\mathbf{Z})} \right] \end{align}$
  - survival function is accelerated by $\exp{(-\boldsymbol{\gamma}^t\mathbf{Z})}$ (change time scale)
  - $h_0(x|\mathbf{Z}) = h_0\left[ x\exp{(-\boldsymbol{\gamma}^t\mathbf{Z})} \right]\exp{(-\boldsymbol{\gamma}^t\mathbf{Z})}$ (accelerated hazard function)
- charateristics
  - provides a direct extension of the classical linear model's construction
  - use is restricted by the error distribution $W$

**Ⅱ. Hazard rate model approach**

idea: <u>model the conditional hazard rate as a function of the covariates</u>

**Ⅱ.1 multiplicative hazard rate model**

- target: $h(x|\mathbf{z})$

- model: $h(x|\mathbf{z}) = h_0(x)c(\boldsymbol{\beta}^t\mathbf{z})$

  - $h_0(x)$ : specified parametric form (or arbitrary nonnegative function)
  - $c(\cdot)$: link function
  - $\mathbf{z}$ : covariate vector
  - $\boldsymbol{\beta}$ : coefficient vector of $\mathbf{z}$
  - common choice of $c(\boldsymbol{\beta},\mathbf{z})=\exp{(\boldsymbol{\beta}^t\mathbf{z})}$ : **Cox's model**

- characteristics

  - with fixed time 0, hazard rates are proportional. i.e.

    $\cfrac{h(x|\mathbf{z}_1)}{h(x|\mathbf{z}_2)}=\cfrac{h_0(x)c(\boldsymbol{\beta}^t\mathbf{z}_1)}{h_0(x)c(\boldsymbol{\beta}^t\mathbf{z}_2)}=\cfrac{c(\boldsymbol{\beta}^t\mathbf{z}_1)}{c(\boldsymbol{\beta}^t\mathbf{z}_2)}$ (constant independent of time)

  - $S(x|\mathbf{z})$ (baseline survival function) $= S_0(x)^{c(\boldsymbol{\beta}^t\mathbf{z})}$

  - "Lehmann Alternative" : relationship in nonparametric statistics above

- Usage

  - modeling relative survival (6.3장)
  - form the basis for modeling covariate effects (8,9장)

- Note

  - Weibull is the <u>only continuous distribution</u> with has the property of being <u>both an AFT and multiplicative hazards</u> model (그래서 자주 쓰이는 듯)

**Ⅱ.2 additive hazard rate model**

- target: $h(x|\mathbf{z})$
- model: $\displaystyle h(x|\mathbf{z}) = h_0(x) + \sum_{j=1}^{p}{z_j(x)\beta_j(x)}$
  - $z_j,\ \beta_j$ : **time-dependent**
- estimation: via weighted least-squares method
- Usage
  - modeling excess mortality (6.3장)
  - modeling regression effects (10장)



#### 2.7 Models for Competing Risks

**definitions**

- $T = \min(X_1,\cdots,X_p)$ : failure time from any cause ($p$ : number of possible causes)
- $\delta = i$ if $T = X_i$
- $\displaystyle h_i = \lim_{\Delta t\rightarrow 0}{\cfrac{P\left[ t\le T < t+\Delta t, \delta = i \ | \ T \ge t \right]}{\Delta t}}$ : cause-specific hazard rate
  - rate at which subjects who have yet to experience any of the competing risks are experiencing the $i$th competing cause of failure
- $\displaystyle h_T(t) = \sum_{i=1}^K{h_i(t)}$ : overall hazard rate of the time-to-failure
- $S(t_1,\cdots,t_k)=Pr\left[ X_1 > t_1,\cdots,X_K > t_K \right]$ : joint survival function
  - $h_i = \cfrac{-\partial S(t_1,\cdots,t_K)/\partial t_i |_{t_1=\cdots=t_K=t}}{S(t,\cdots,t)}$

**Assumptions on dependency**

- independent competing risks: marginal hazard = cause-specific hazard

  $h(t) = h_i(t)\ \ \forall i$

- dependent competing risks: need to make some assumptions about the dependence structure b/w potential failure times

- identifiability dilemma: 우리가 관측 가능한 것은 사망 시간과 원인 뿐이며, 이를 통해 dependence structure에 대한 assumption이 합당한지 competing risk data만으로는 testing이 불가능함

**likelihood approach**

- competing risk data에서는 (cause-specific) hazard rate와 더불어 특정 원인의 발생 가능성 (likelihood)이 관심 대상이다
- crude probability
  - $F_i(t) = Pr\left[ T \le t, \delta = i \right]$ : 모든 cause가 발생 가능할 때 특정 cause에 의해 사망할 확률
  - cumulative incidence function이라고도 함
  - 예시: 어떤 남자가 50세 이전에 심장마비로 사망할 확률
  - $\displaystyle F_i(t)=\int_0^t{h_i(u)\exp{\left\{ -H_T(u) \right\}}du}$
  - depending all the competing risks (not simply on specific cause)
- net probability
  - $S_i(t)$ : 오직 한 가지 원인에 의해서만 사망한다는 가정 하에서의 생존함수
  - marginal probability와 같음
  - 예시: 다른 원인이 발생하지 않는다는 가정 하에서, 어떤 남자가 50세 이전에 심장마비로 사망할 확률
  - independent / dependent 여부에 따라 구하는 방법이 달라짐
- partial crude probability
  - $\lambda_i^J(t)$: 일부 사망 cause가 이론적으로 제거된 상태에서의 crude prob.
    - $J$ : indicator of remaining set of causes

**Notes**

- crude probability $\ne$ net probability

