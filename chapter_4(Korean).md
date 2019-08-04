## 4. Nonparametric Estimation of Basic Quantities for Right-Censored and Left-Truncated Data

#### 4.1 Introduction

**Basic quantities**

- $d_i$ : number of event at time $t_i$
- $Y_i$ : number of observations <u>at risk</u> at time $t_i $
- $\cfrac{d_i}{Y_i}$ : estimate of the conditional probability $Pr[T= t_i\ |\ T>t_i]$
  - basic quantity of construction of estimators of the $S(t)$ and $H(t)$

#### 4.2 Estimators of the Survival and Cumulative Hazard Functions for Right-Censored Data

- **Product-Limit estimator** (discretized Kaplan-Meier estimator)
  - $\hat{S}(t) = \begin{cases}1, & \text{if }t<t_1 \\ \prod_{t_i\le t}[1-{d_i\over Y_i}], & \text{if }t_1 \le t\end{cases}$
  - $\displaystyle \hat{V}[\hat{S}(t)] = \hat{S}(t)^2\sum_{t_i\le t}\cfrac{d_i}{Y_i(Y_i-d_i)}=\hat{S}(t)^2\sigma_S^2(t)\ \text{(Greenwood's formula)}$
  - $d_i$ : number of events(death) at time $t_i$
  - $Y_i$ : number at risk at time $t_i$
    - censoring이 발생할 경우 $Y_i$에서 제외해야함!
  - $\hat{H}(t) = -\ln[\hat{S}(t)]$
- **Nelson-Aalen estimator** (using modern counting process techniques)
  - $\tilde{H}(t) = \begin{cases}0, & \text{if }t\le t_1 \\ \sum_{t_i \le t}{d_i\over Y_i}, & \text{if }t_1 \le t\end{cases}$
  - $\displaystyle \hat{\sigma}_H^2(t)=\hat{V}[\tilde{H}(t)]=\sum_{t_i \le t}\cfrac{d_i}{Y_i^2}$
  - $\tilde{S}(t) = \exp[-\tilde{H}(t)]$
  - **primary uses in analyzing data**
    - selection among parametric models for the time-to-event data
    - providing crude estimates of the hazard rate $h(t)$ (slope of $\tilde{H}(t) $)
- **Notes**
  - assumption: knowledge of a censoring time for and individual provides <u>no further information about this person's likelihood</u> of survival at a future time had the individual continued on the study
  - Kaplan-Meier estimator는 관측된 시간의 최대값 이전 사건들에 대해서 잘 정의되지만, 그 이후의 사건에 대해서는 성능이 좋지 않음



#### 4.3 Pointwise Confidence Intervals

**Confidence Intervals for survival function** $S(t) $

- Basic statistic: **product-limit estimator** ($\hat{S}(t)$) and its standard error

- linear confidence interval
  - $\text{CI} = \hat{S}(t_o) \pm z_{1-\alpha/2}\sigma_S(t_o)\hat{S}(t_o)$
- 1st transformed CI: **log-transformation** confidence interval
  - $\text{CI} = [\hat{S}(t_o)^{1/\theta},\ \hat{S}(t_o)^{\theta}],\ \text{where }\theta=\exp\left\{ \cfrac{Z_{1-\alpha/2}\sigma_S(t_o)}{\ln[\hat{S}(t_o)]} \right\}$
- 2nd transformed CI: **arcsine-square root transformation** confidence interval
  - $\text{LB: } \sin^2\left\{ \max\left[ 0,\arcsin(\hat{S}(t_o)^{1/2})-0.5Z_{1-\alpha/2}\sigma_S(t_o)\left( \cfrac{\hat{S}(t_o)}{1-\hat{S}(t_o)} \right)^{1/2} \right] \right\}$
  - $\text{UB: } \sin^2\left\{ \min\left[ \cfrac{\pi}{2},\arcsin(\hat{S}(t_o)^{1/2}) + 0.5Z_{1-\alpha/2}\sigma_S(t_o)\left( \cfrac{\hat{S}(t_o)}{1-\hat{S}(t_o)} \right)^{1/2} \right] \right\}$

**Confidence Intervals for cumulative hazard function** $H(t) $

- Basic statistic: **Nelson-Aalen estimator** ($\tilde{H}(t)$) and its standard error
- $1-\alpha$ linear confidence interval
  - $\text{CI} = \tilde{H}(t_o) \pm Z_{1-\alpha/2}\sigma_H(t_o)$
- log-transformation confidence interval
  - $\text{CI}=[\tilde{H}(t_o)/\phi,\phi\tilde{H}(t_o)],\ \text{where }\phi=\exp\left[\cfrac{Z_{1-\alpha/2}\sigma_H(t_o)}{\tilde{H}(t_o)}\right]$
- arcsine-square root transformation confidence interval
  - $\text{LB}$ : $-2\ln\left\{ \sin\left[ \min\left( \cfrac{\pi}{2}, \arcsin\left[\exp\left(-\cfrac{\tilde{H}(t_o)}{2}\right)\right]+0.5Z_{1-\alpha/2}\sigma_H(t_o)\left( \exp \tilde{H}(t_o)-1 \right)^{-1/2} \right) \right] \right\}$
  - $\text{UB}$ : $-2\ln\left\{ \sin\left[ \max\left( 0, \arcsin\left[\exp\left(-\cfrac{\tilde{H}(t_o)}{2}\right)\right]-0.5Z_{1-\alpha/2}\sigma_H(t_o)\left( \exp \tilde{H}(t_o)-1 \right)^{-1/2} \right) \right] \right\}$

**Notes**

- log-transformed & arcsine-square root transformed CI **perform better** than linear (Borgan, Liestol (1990)) with smaller sample size and large censoring individuals
- log-transformed & arcsine-square root transformed CI are **not symmetric**, thus appropriate for small samples where the <u>point estimators are biased</u> and the <u>distribution of the estimators is skewed</u>
- suggested CIs above are not confidence band (not for all values $t$)



#### 4.4 Confidence Bands for the Survival Function

**Basic idea**

- 특정 시점 $t_o$가 아닌 추정하고자 하는 시간 구간 $[t_L,t_U]$에 대한 신뢰대역(confidence band)를 세우고자 함
- confidence interval에 비례하는 방법과 비례하지 않는 방법의 2가지 접근법이 있음
- 두 접근법 모두 세 방법(linear, log-trans, arcsine-square root trans)에 의해 CB를 계산할 수 있음

**1) EP(equal probability) band**

- implementation: 모든 시점 $t\in [t_L,t_U]$에 대해, pointwise CI 계산식의 critical value $Z_{1-\alpha/2}$를 $c_{\alpha}(a_L,a_U)$로 바꾸어 계산하면 됨
  - pointwise CI에 비례함
  - 주어진 세 가지 방법(linear, log-trans, arcsine-square root trans)에 모두 동일하게 적용 가능
- $c_{\alpha}(a_L,a_U)$ : confidence coefficient, 미리 계산된 값들 (table이 존재함)
  - $a_L = \cfrac{n\sigma_S^2(t_L)}{1+n\sigma_S^2(t_L)},\ \ a_U = \cfrac{n\sigma_S^2(t_U)}{1+n\sigma_S^2(t_U)}$
  - $n$ : given sample size
  - $t_L$ : 관측된 event-time의 최소값보다 작아서는 안됨 + 0보다 커야함
    $t_U$ : 관측된 event-time의 최대값보다 커서는 안됨

**2) Hall-Wellner confidence band**

- implementation: 모든 시점 $t \in [t_L, t_U]$에 대해, pointwise CI 계산식의 $Z_{1-\alpha/2}\sigma_S(t)$를 $k_{\alpha}(a_L,a_U)\cfrac{1+n\sigma_S^2(t)}{n^{1/2}}$로 대체하여 계산함
  - pointwise CI에 비례하지 않음
  - 주어진 세 가지 방법(linear, log-trans, arcsine-square root trans)에 모두 동일하게 적용 가능
- $k_{\alpha}(a_L,a_U)$ : confidence coefficient, 미리 계산된 값들 (table이 존재함)
  - $a_L, a_U$는 위에서 주어진 식과 동일함
  - $t_L = 0$로 설정할 수 있음

**Notes**

- Pointwise CI보다 더 넓은 구간을 추정하는 경향이 있음
- Hall-Wellner CB는 작은 $t$값에서는 더 넓은, 큰 $t$값에서는 더 좁은 band를 도출
- Cumulative hazard function의 CB도 위와 같은 방법으로 계산할 수 있음
  ($Z_{1-\alpha/2}$ 또는 $Z_{1-\alpha/2}\sigma_H(t)$ 대체하여 계산)
- linear CB는 sample size가 작을 때 $(<200)$ 성능이 좋지 않다고 함, 또한 arcsine transformation이 성능이 약간 더 좋다고 함



#### 4.5 Point and Interval Estimates of the Mean and Median Survival Time

**Ideas**

- MRL이나 MDRL, $E(X)$, $Var(X)$ 등 quantity의 nonparametric estimator는 $\hat{S}(t)$ 또는 $\tilde{H}(t)$를 대입하여 바로 도출할 수 있음

**Mean survival time** ($E(X)$) **estimator**

- $\displaystyle \hat{\mu}=\int_0^{\infty}\hat{S}(t)dt,\ \ \text{where}\ \hat{S}(t)\text{ : Product-Limit estimator}$
- Corrected estimator for $\mu$ : 마지막 관찰시점이 censored value일 경우
  - basic form: $\displaystyle \hat{\mu}_{\tau}=\int_0^{\tau}\hat{S}(t)dt$ (**interval이 restricted 됨**)
  - Efron's tail correction: 마지막 관찰시점을 largest observed death time $t_{\max}$으로 대체하여 계산하는 방법
  - preassigned interval: 마지막 관찰시점을 임의의 시점 $\tau$로 설정하여 계산하는 법 (모든 개체가 생존할 수 있는 시점으로 설정)
  - $\displaystyle \hat{V}[\hat{\mu}_{\tau}] = \sum_{i=1}^D\left[ \int_{t_i}^{\tau}\hat{S}(t)dt \right]^2 \cfrac{d_i}{Y_i(Y_i - d_i)}$
  - $1-\alpha$ CI for $\mu$ : $\hat{\mu}_{\tau} \pm Z_{1-\alpha/2}\sqrt{\hat{V}[\hat{\mu}_{\tau}]}$

$p$**th quantile of survival time** ($x_p$) **estimator**

- $\hat{x}_p=\inf\left\{ t:\hat{S}(t)\le 1-p \right\}$
- Confidence intervals for $x_p $
  - 위에서와 마찬가지로 linear, log-trans, arcsine-square root trans의 3가지 방법으로 구할 수 있음
  - 아래 조건들을 만족하는 모든 $t$의 집합들이 $x_p $의 CI가 됨
  - **linear**: $-Z_{1-\alpha/2}\le \cfrac{\hat{S}(t)-(1-p)}{\hat{V}^{1/2}[ \hat{S}(t) ]} \le Z_{1-\alpha/2}$
  - **log-trans**: $-Z_{1-\alpha/2}\le \cfrac{ [\ln(-\ln\hat{S}(t)) - \ln(-\ln(1-p))]*\hat{S}(t)*\ln\hat{S}(t) }{\hat{V}^{1/2}[ \hat{S}(t) ]} \le Z_{1-\alpha/2}$
  - **arcsine-square root trans**: $-Z_{1-\alpha/2}\le \cfrac{ 2\left[ \arcsin\left(\sqrt{\hat{S}(t)}\right) - \arcsin\left(\sqrt{1-p}\right)\right]*\{\hat{S}(t)(1-\hat{S}(t))\}^{1/2} }{\hat{V}^{1/2}[ \hat{S}(t) ]} \le Z_{1-\alpha/2}$



#### 4.6 Estimators of the Survival Function for Left-Truncated and Right-Censored Data

**Basic Ideas**

- right-censored에 더해 "<u>left-truncation</u>"이 발생한 data를 handling하는 기법을 다룸

- $j$ : $j$th individual

- $L_j$ : study entering time of individual $j$

- $T_j$ : death or censoring time of individual $j$

- $t_1 < t_2 < \cdots <t_i<\cdots < t_D$ : ordered death times (not censoring)

- $d_i$ : number of death at time $t_i$

- $Y_i$ : remaining individuals at risk at time $t_i$

  - **(left-truncated data)** $Y_i = \text{number of individuals s.t. }L_j < t_i \le T_j$

  - 이 지점이 left-truncation이 발생했을 경우의 차이점

- redefined $Y_i$를 사용하여 앞서 논의한 추정량들을 거의 바로 적용할 수 있음

**Product-Limit estimator & Nelson-Aalen estimator**

- Product-Limit estimator: $\hat{S}(t)$의 해석이 $Pr[X>t \ | \ X\ge L] = S(t)/S(L) \ \ (L\text{ : the smallest of the entry times})$의 조건부 확률 추정량으로 바뀜
- Nelson-Aalen estimator: $\tilde{H}(t)$의 해석이 $\displaystyle \int_L^th(t)dt$의 추정량으로 바뀜 (단, 조건부 확률의 추정량은 아님)

**Notes**

- independent truncation 가정이 필요함
- left-truncation이 발생한 경우 Product-Limit estimator는 작은 $t$에 대해 큰 분산을 가지는 경향이 있음



#### (4.7 생략)