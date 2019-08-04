## 3. Censoring and Truncation

#### 3.1 Introduction

- censoring: 모든 구간에서 관측이 가능하지만 그 값을 정확히 알 수 없고, 특정 구간에서 사건이 발생했다는 사실만 알 수 있는 경우
  - 부분적인 정보를 포함하고 있기 때문에 censoring이 발생한 개체에 대해 <u>partial likelihood construction</u>이 가능함
- truncation: 특정 구간에서의 관측 자체가 불가능한 경우, 또는 특정 구간을 의도적으로 조사하지 않은 경우 (30세 이상 생존한 개체만을 대상으로 조사하는 경우 등)
  - 따라서 likelihood construction에 <u>conditional probability</u>가 사용됨 (특정 구간을 제외하는 조건)
- censoring과 truncation은 design(sampling)에 따라 다른 likelihood function을 가짐
- 서로 다른 likelihood structure를 대상으로 한 common approach가 존재함
- likelihood construction: parametric modeling의 기반
- partial likelihood construction: semi-parametric modeling의 기반
- counting process
  - censored / truncated data에 대한 nonparametric techniques의 기반
  - survival analysis에 대한 parametric / nonparametric 방법론의 중요한 성질



#### 3.2 Right Censoring

**Terminology**

- right censored data : $(T,\delta)$
- $T = \min(X,C_r)$
  - $X \overset{\text{iid}}{\sim} f(x),\ S(x)$ : random survival lifetime
  - $C_r$ : right censored lifetime of the individuals
- $\delta = \begin{cases}1, & T=X \\ 0,&T=C_r\end{cases}$

**Type 1 censoring**

- Type 1 censoring: event is observed only if it occurs prior to some prespecified ("fixed") time
  - $T$ is **random** if $T=X$
- progressive type 1 censoring : individuals have ①different ②fixed censoring time
  - advantage is that the censored individuals give information on the natural history of nonlethal event
- generalized type 1 censoring : individuals ①enter the study at different times and has ②different ③fixed censoring time
  - shift each individual's starting time to 0
  - or using Lexis diagram

**Type 2 censoring**

- Type 2 censoring: continue study until the failure of the first $r$ individuals
  - often used in testing of equipment life
  - may save time and money because it could take a very long time for all items to fail.
  - statistical treatment of type 2 censored data is simpler, because theory of order statistics is directly applicable to determining the likelihood
  - number of censored observations $r$ is **fixed**
  - observer censoring time $T_{(r)}$ is **random**
- progressive type 2 censoring (generalized type 2 censoring)
  - repeated type 2 censoring, removing some remaining individuals
  - each observed censoring time $T_{(i)}$ are **random**

**Competing risks censoring**

- type of random censoring (복수의 알 수 없는 원인에 의해 무작위로 censoring이 발생함을 가정)
- 특정 cause의 marginal distribution에 대한 추정이 필요하지만, 몇몇의 개체가 competing risks를 보유하고 있는 경우 필요한 가정임
- 해당하는 개체가 사망한 경우, 관심있는 cause가 아닌 "돌발상황"에 의해 실험에서 이탈한 상황으로 간주하여 **random right censored** 처리함 (Censoring time이 random이 됨)
  $\rightarrow $ 따라서 $X$(marginal event)와 censoring time이 서로 독립임을 가정해야 함
- quantity estimation
  1. marginal probability (net probability)
     - other competing risks are considered as random observations except for the target event
     - need and assumption of independence b/w target event vs competing events
  2. crude probability
     - each competing risk is modeled by an cumulative incidence curve
     - no independence assumption is needed
  3. partial crude probabilities
     - competing risks to be eliminated: treated as random obs.
     - competing risks to be remained: modeled by a cumulative incidence curve



#### 3.3 Left or Interval censoring

**Left censoring**

- left censored data : $(T,\epsilon)$
- $T = \max(X,C_l)$
  - $X \overset{\text{iid}}{\sim} f(x),\ S(x)$ : random survival lifetime
  - $C_l$ : left censored lifetime of the individuals
- $\epsilon = \begin{cases}1, & T=X \\ 0,&T=C_l\end{cases}$

**Doubly censoring**

- doubly censored data : $(T,\delta)$
- $T = \max[\min(X,C_r),\ C_l]$
  - $C_r$ : time after some individuals experience the event
  - $C_l$ : time before some individuals experience the event
  - $C_l \le x \le C_r$
- $\delta = \begin{cases} 1,&T=X \\ 0,& T = C_r \\ -1. & T = C_l \end{cases}$

**Interval censoring**

- general expression of censored data
- only known to the event time fall in an interval $(L_i,R_i]$
- periodic inspection for proper functioning of items are appropriate



#### 3.4 Truncation

**General ideas**

- Truncation: only the events within $(Y_L,Y_R)$ are observed, and others are not
- not observed individuals have no information to the investigator
- main contrast to censoring, which is at least partial information available (events may occur in certain interval)
- inference is restricted to <u>conditional estimation</u>

**Left truncation**

- left truncation: $Y_R \rightarrow \infty$
- often called delayed entry time
- the event of interest prior to the truncation time are not observed
- truncated individuals are never considered for inclusion into the study

**Right truncation**

- right truncation: $Y_L \rightarrow 0$
- particularly relevant to studies of AIDS



#### 3.5 Likelihood Construction for Censored and Truncated Data

**basic ideas**

- available informations from observation

  $\begin{matrix}\text{exact lifetimes} & - & f(x) \\ \text{right-censored obs.} & - & S(C_r) \\ \text{left-censored obs.} & - & 1-S(C_l) \\ \text{interval-censored obs.} & - & [S(L) - S(R)] \\ \text{left-truncated obs.} & - & f(x)/S(Y_L) \\ \text{right-truncated obs.} & - & f(x)/[1-S(Y_R)] \\ \text{interval-truncated obs.} & - & f(x)/[S(Y_L) - S(Y_R)]\end{matrix}$

  - $X$ : approximate density function of $X$
  - Censored obs. : info. on survival function
  - truncated obs. : info. on conditional density of $X$

- basic likelihood function

  - $\displaystyle L \propto \prod_{i \in D}f(x_i)\prod_{i \in R}S(C_r)\prod_{i \in L}(1-S(C_l))\prod_{i \in I}[S(L_i) - S(R_i)]$
  - $D$ : the set of (observed) death times
  - $R$ : the set of right-censored obs.
  - $L$ : the set of left-censored obs.
  - $I$ : the set of interval-censored obs.

- left-truncated data

  - $f(x)\ \rightarrow \  \cfrac{f(x)}{[S(Y_{Li}) - S(Y_{Ri})]}$
  - $S(C_i) \ \rightarrow \ \cfrac{S(C_i)}{[S(Y_{Li}) - S(Y_{Ri})]}$
  - $(Y_{Li}, Y_{Ri})$ : truncation interval

- right-truncated data

  - $\displaystyle L \propto \prod_{i}f(Y_i)/[1-S(Y_i)]$
  - $Y_i$ : observed death time

- different failure distribution

  - $\displaystyle L = \prod_{i \in D}f_i(x_i)\prod_{i \in R}S_i(C_r)\prod_{i \in L}(1-S_i(C_l))\prod_{i \in I}[S_i(L_i) - S_i(R_i)]$

**Likelihood construction for right-censored data**

- **type 1 censoring (ordinary right-censoring)** 
  - $\displaystyle L_Ⅰ=\prod_{i=1}^n{Pr[t_i,\delta_i]} = \prod_{i=1}^n[f(t_i)]^{\delta_i}[S(t_i)]^{1-\delta_i}$
  - $t_i$ : observed death time or censored time
  - $\displaystyle L_Ⅰ= \prod_{i=1}^n[h(t_i)]^{\delta_i}\exp{[-H(t_i)]}$
- **type 2 censoring (ordered observations)**
  - simple case (data consist of the $r$th smallest lifetimes)
    - $\displaystyle L_{Ⅱ,1} = \cfrac{n!}{(n-r)!}\left[ \prod_{i=1}^n f(x_{(i)}) \right][S(x_{(r)})]^{n-r}$
  - progressive type 2 censoring (2 repeated study)
    - $\displaystyle L_{Ⅱ,2} \propto \prod_{i=1}^{r_1}f(x_{(i)})[S(x_{(r_i)})]^{n_1}\prod_{i=1}^{r_2}f(x_{(i)}^*)[X(x_{(r_2)}^*)]^{n-r_1-n_1-r_2}$
    - $x^*$ : obs. from second study
  - likelihood can be constructed by the theory of order statistics
- **type 3 censoring (random censoring)**
  - in the type 3 censoring, censoring time have to be considered as random variable (since cause of event cannot be revealed)
  - dealing with joint density function $f(x,c_r)$
  - marginal density of $C_r \sim g(c_r),\ G(c_r)$
  - independent $X$ and $C_r$
    - $\displaystyle L_{Ⅲ,1}=\left\{ \prod_{i=1}^nG(t_i)^{\delta_i}g(t_i)^{1-\delta_i} \right\}\left\{ \prod_{i=1}^n f(t_i)^{\delta_i}S(t_i)^{1-\delta_i} \right\}$
  - dependent $X$ and $C_r$
    - $\displaystyle L_{Ⅲ,2} \propto \prod_{i=1}^n\left\{ \left[ \cfrac{-\partial S(x,t_i)}{\partial x} \right]_{x=t_i} \right\}^{\delta_i}\left\{ \left[ \cfrac{-\partial S(t_i,c)}{\partial c} \right]_{c=t_i} \right\}^{1-\delta_i}$



#### 3.6 Counting Processes

**Definitions**

- $N(t)\ (t\ge 0)$ : counting process
  - $N(0) = 0, \ N(t)<\infty \ \ \forall t$
- $\displaystyle N(t) =\sum_{i=1}^n{N_i(t)}=\sum_{t_i\le t}{\delta_i}$ (for right-censored data)
  - $\delta_i$ : indicator for the death of $i$th individual
  - $N_i(t)=I[T_i \le t,\ \delta_i = 1]$
  - simply counts the number of deaths in the sample at (or prior to) $t$
- $\mathbf{F}_t$ : history of the counting process at time $t$
  - accumulated knowledge about what has happened to patients up to time $t$
  - $\mathbf{F}_{t-}$ : an instant just prior to time $t$

**Basic counting (stochastic) processes and quantities**

- $dN(t)$ : the change in the process $N(t)$ over a short time interval $[t,t+dt)$
  - $dN(t) = \begin{cases}1, & \text{death occurred at }t\\ 0, & \text{o/w}\end{cases}$
- $Y(t)$ : the number of individuals with a study time $T_i \ge t$
  - provides the number of individuals **at risk** at a given time
- $\lambda(t)=Y(t)h(t)$ : intensity process (of the counting process)
  - $E[dN(t)|\mathbf{F}_{t-}] = \lambda(t)dt$
- $\Lambda(t) = \int_{0}^t\lambda(s)ds$ : cumulative intensity process
  - $E[N(t)|\mathbf{F}_{t-}] = \Lambda(t)$
- $M(t)=N(t) - \Lambda(t)$ : counting process martingale
  - $E[dM(t)|\mathbf{F}_{t-}] = 0$
  - $Var[dM(t)|\mathbf{F}_{t-}]=d\left\langle M \right\rangle(t)$
    - $\left\langle M \right\rangle(t)$ : predictable variation process
  - interpretation
    - $N(t)$ : target nondecreasing step function
    - $\Lambda(t)$ : compensator of the $N(t)$ (smoothing process)
    - martingale can be considered as mean zero noise

**Nelson-Aalen estimator**

- $\hat{H}(t)\text{ (Nelson-Aalen estimator)} = H^*(t)+W(t)$
  - basic relationship: $dN(t) = Y(t)h(t)dt + dM(t)$
    - $t$시점에서의 사망자 수는 ($t$시점에서의 생존자 수)*($t$시점에서의 hazard) + random error
    - $\mathbf{F_{t-}}$에서의 $Y(t)$는 상수 취급
  - **formal definition**: $\displaystyle \underset{\hat{H}(t)}{\int_0^tJ(u)\cfrac{dN(u)}{Y(u)}} = \underset{H^*(t)}{\int_0^tJ(u)h(u)du} + \underset{W(t)}{\int_0^t\cfrac{J(u)}{Y(u)}dM(u)}$
    - $J(t)$ : indicator of whether $Y(t)$ is positive value
    - 위의 basic relationship으로부터 도출된 관계
    - 정의에 의해 $\hat{H}(t)$는 **cumulative hazard function의 nonparametric estimator**가 됨
  - **limiting distribution**: $\displaystyle \sqrt{n}[\hat{H}(t)-H^*(t)] \underset{\cdot}{\overset{\cdot}{\sim}}N(0,\sigma^2),\ \ \sigma^2 = \int_0^t\cfrac{h(u)du}{y(u)}$
    - this result can be justified by **martingale central limit theorem**
    - employed to evaluate confidence interval or band of $\hat{H}(t)$
    - $\displaystyle \hat{\sigma}^2=n\int_0^t\cfrac{dN(u)}{Y(u)^2}$

**Kaplan-Meier estimator**

- $\displaystyle \hat{S}(t)\text{ (Kaplan-Meier estimator)} = \prod_{s=0}^t[1-d\hat{H}(t)]=\prod_{s=0}^t\left[ 1-\cfrac{dN(s)}{Y(s)} \right]$
  - 정의에 의해 $\hat{S}(t)$는 **survival function의 nonparametric estimator가 됨**
  - discrete case expression
  - confidence interval and band can be found using limiting distribution suggested above

**Likelihood construction using counting process**

- $N_j(t)$ : individual counting process ($0$ or $1$)
- $dN_j(t) \underset{\cdot}{\overset{\cdot}{\sim}}\text{Bernoulli}(\lambda_j(t))$ (특정 시점에 사망할 intensity가 확률이 됨)
- **full likelihood**: $\displaystyle L=\left[ \prod_{j=1}^n\lambda_j(t)^{dN_j(t)} \right]\exp\left[ -\sum_{j=1}^n\int_0^{\tau}\lambda_j(u)du \right]$
  - $\tau$ : upper limit of time range
- **right-censored likelihood**: $\displaystyle L \propto \left[ \prod_{j=1}^n h(t_j)^{\delta_j} \right]\exp\left( -\sum_{j=1}^nH(t_j) \right)$
  - using relationship $\lambda_j(t) = Y_j(t)h(t)$
  - $Y_j(t) = \begin{cases}1,& t\le t_j \\ 0, & t > t_j\end{cases}$
