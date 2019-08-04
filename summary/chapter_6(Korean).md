## 6. Topics in Univariate Estimation

#### 6.1 Introduction

- KM, NA estimator는 한정된 정보만을 제공함 → NA estimator의 slope가 제공하는 <u>crude estimator of hazard function의 해석</u>을 풍부하게 하는 방법을 알아보자 (kernel-smoothing method)
- excess mortality를 표현하는 두 모형을 알아보자 (additive & multiplicative)
- Bayesian perspective에서 right-censored data를 handling하는 방법을 알아보자 (analytic solution / Gibbs sampler)



#### 6.2 Estimating the Hazard Function

**Ideas**

- NA estimator는 $H(t)$에 대한 추정량을 제공하지만, 대부분 $H(t)$보다는 $h(t)$에 대해 관심이 많음.
- $h(t)$는 $H(t)$의 미분값이므로 일종의 기울기로 해석할 수 있음. 그러므로 $H(t)$ 추정량의 미분값을 계산할 수 있으면 $h(t)$의 crude estimator를 얻을 수 있음
- 따라서 discretized estimator (step function)인 $\tilde{H}(t)$를 **smoothing시켜서 미분 가능하게 만들자**는 것이 핵심 (교재에서는 kernel-smoothing 방법 사용)

**Basic Quantities**

- $\tilde{H}(t),\ \hat{V}(\tilde{H}(t))$ : NA estimator and its sample variance
- $0=t_0<t_1<\cdots<t_D$ : event times (jumping points)
- $\Delta\tilde{H}(t_i)=\tilde{H}(t_t)-\tilde{H}(t_{i-1})$
  - magnitude of the jumps in $\tilde{H}(t_i)$
- $\Delta\hat{V}(\tilde{H}(t_i))=\hat{V}(\tilde{H}(t_i))-\hat{V}(\tilde{H}(t_{i-1}))$
  - magnitude of the jumping variance in $\tilde{H}(t_i)$
- $b$ : bandwidth of kernel-smoothed estimator
  - range to determine smoothness of the estimator
- $K()$: kernel function
  - function to determine weight of closeness to time $t$
  - defined on the interval $[-1, 1]$
  - uniform kernel: $K(t) = 1/2$
  - Epanechnikov kernel: $K(t) = 0.75(1-t^2)$
  - biweight kernel : $K(t) = 15/16(1-t^2)^2$

**Building Kernel-smoothed estimator of** $h(t)$

- $(b \le t \le t_D-b)$ 인 경우

  - $\displaystyle \hat{h}(t)=\cfrac{1}{b}\sum_{i=1}^D K\left(\cfrac{t-t_i}{b}\right)\Delta\tilde{H}(t_i)$

  - $\displaystyle \hat{V}(\hat{h}(t))=\cfrac{1}{b^2}\sum_{i=1}^DK\left(\cfrac{t-t_i}{b}\right)^2\Delta\hat{V}(\tilde{H}(t_i))$

- $(t<b)$ 인 경우

  - 위 식의 kernel function을 $K_q$ (modified kernel function, Gasser and Muller(1979))으로 대체

  - estimator와 variance에 모두 적용 가능

  - $q=t/b$

  - modified uniform kernel

    $K_q(t) = \cfrac{4(1+q^3)}{(1+q)^4}+\cfrac{6(1-q)}{(1+q)^3}t$

  - modified Epanechnikov kernel

    $K_q(t) = K(t)(\alpha_E+\beta_Et)$

    $(\alpha_E=\cfrac{64(2-4q+6q^2-3q^3)}{(1+q)^4(19-18q+3q^2)}, \ \beta_E=\cfrac{240(1-q)^2}{(1+q)^4(19-18q+3q^2)})$

  - modified biweight kernel

    $K_q(t) = K(t)(\alpha_{BW}+\beta_{BW}t)$

    $\begin{align}(\alpha_E=\cfrac{64(8-24q+48q^2-45q^3+15q^4)}{(1+q)^5(81-168q+126q^2-40q^3+5q^4)},\\ \beta_E=\cfrac{1120(1-q)^3}{(1+q)^5(81-168q+126q^2-40q^3+5q^4)})\end{align}$

- $(t>t_D-b)$ 인 경우 (upper tail correction)

  - modified kernel function의 $t$를 $-t$로 대체하여 계산하면 됨
  - $q=(t_D-t)/b$

- Confidence Interval & Confidence band

  - CI와 CB는 4장의 방법을 그대로 사용하면 됨

- Notes

  - $\hat{h}(t)$는 **smoothed function의 추정량**이기 때문에 (그래서 crude estimator) 해석에 주의해야 함 (true hazard rate의 CI / CB가 아님)
  - biweight > Epanechnikov > uniform 순으로 smooth한 경향이 있음
    (uniform이 제일 jagged)
  - bandwidth가 클수록 smooth한 경향이 있으나, <u>bias가 함께 증가하여 optimization이 필요</u>함

**Finding optimal bandwidth** $b$ **with cross-validation (CV)**

- cross-validation technique을 적용하기 위해
  ① mean integrated square error(MISE)를 정의하고
  ② 그 추정량(목적함수) $g(b)$을 계산한 다음
  ③ $g(b)$을 최소화하는 $b$ 값을 최적의 bandwidth로 채택하여 사용한다
- MISE
  - $(\tau_L,\tau_U)$ : MISE를 정의하는 시구간
  - 
    $\begin{align}MISE(b) = &E\int_{\tau_L}^{\tau_U}[\hat{h}(u) - h(u)]^2du\text{ (kind of squared error)} \\ =& E\int_{\tau_L}^{\tau_U}\hat{h}^2(u)du - 2E\int_{\tau_L}^{\tau_U}\hat{h}(u)h(u)du + E\int_{\tau_L}^{\tau_U}h^2(u)du \end{align}$
  - 첫 번째 항과 두 번째 항이 $b$에 대한 식이므로, 이 두 항을 추정하여 최소화한다
- 목적함수 $g(b)$
  - $\tau_L=u_1 < \cdots<u_M = \tau_U$ : grid points (CV의 정밀함, $M$개에 대해 탐색)
  - $\displaystyle g(b)=\sum_{i=1}^{M-1}\left( \cfrac{u_{i+1}-u_i}{2} \right)[\hat{h}^2(u_i)+\hat{h}^2(u_{i+1})] - 2\cfrac{1}{b}\sum_{i \ne j}K\left( \cfrac{t_i - t_j}{b} \right)\Delta\tilde{H}(t_i)\Delta\tilde{H}(t_j)$
- optimal value $b^*$
  - $b^* = \underset{b}{\arg\min}\ {g(b)}$



#### 6.3 Estimation of Excess Mortality

**Ideas**

- 특정 집단의 사망 가능성이 전체 모집단과 다를 수 있으니, 그 차이의 크기를 비교해보고자 함
- 크게 모집단 사망가능성에 비례한다는 가정**(relative excess mortality model)**과 모집단 사망가능성에 추가로 가능성이 더해진다는 가정**(additive excess mortality model)**으로 나눌 수 있음

**Basic assumptions**

- $n$ : 연구 대상 집단의 표본수

- $\theta_j(t)$ : $j$번째 개체에 대한 $t$시점 모집단 hazard rate (reference rate)

  이미 알려진 값을 사용하며(전체 인구추계 등), 개별 개체 특성(성별, 연령)에 따라 다르게 적용 가능함

- $h_j(t)$ : $j$번째 개체의 실제 hazard rate

  <u>relative 또는 excess model에 의해 보정되어야 하는 hazard</u>임

- $0 \le t_1 < \cdots < t_D$ : event times of data

- $d_i$ : $t_i$에 발생한 number of events

- $Y_j(t)$ : $t$시점의 $j$번째 개체 생존 여부 (살았으면 1, 죽었으면 0)

- $Y(t)$ : $t$시점의 at risk 수 ($=\sum_{j=1}^nY_j(t)$)

- right-censoring + left-truncation data에 적용 가능함

**Relative mortality model**

- **model**: $h_j(t)=\beta(t)\theta_j(t)$

  - $\beta(t)$ : 연구 대상 집단의 실제 relative mortality rate

    이 <u>집단 특성에 비례</u>해 개체의 실제 사망률이 결정된다는 가정

  - $B(t) = \int_0^t\beta(u)du$ : cumulative relative mortality

    데이터를 통한 추정의 대상임

- $B(t)$와 $\beta(t)$의 추정

  - $\displaystyle \hat{B}(t)=\sum_{t_i \le t}d_i/Q(t_i),\ \ Q(t)=\sum_{j=1}^n\theta_j(t)Y_j(t)$

    $\displaystyle \hat{V}(\hat{B}(t))=\sum_{t_i \le t}d_i/Q^2(t_i)$

  - large-sample normal approximation이 성립하므로, NA estimator를 활용한 CI, CB construction method를 그대로 적용할 수 있음

  - $\hat{\beta}(t)$ : kernel-smoothing method등을 활용해 얻을 수 있는 $\hat{B}(t)$의 미분값으로 추정할 수 있음

**Excess mortality model**

- **model**: $h_j(t) = \alpha(t)+\theta_j(t)$

  - $\alpha(t)$ : 연구 대상 집단의 실제 excess mortality rate

  - $A(t)=\int_0^t\alpha(u)du$ : cumulative excess mortality

    마찬가지로 데이터를 통한 추정의 대상임

- $A(t)$와 $\alpha(t)$의 추정

  - $\displaystyle \hat{A}(t)=\sum_{t_i \le t}\cfrac{d_i}{Y(t)}-\Theta(t), \ \ \Theta(t) = \sum_{j=1}^n \int_0^t\theta_j(u)\cfrac{Y_j(u)}{Y(u)}du\text{ (expected }H(t))$

    $\displaystyle \hat{V}(\hat{A}(t))=\sum_{t_i \le t}\cfrac{d_i}{Y(t)^2}$

  - 마찬가지로 large-sample normal approx.와 kernel-smoothing method를 적용하여 CI, CB, $\hat{\alpha}(t)$를 도출할 수 있음

- **corrected survival function**

  - $\Theta(t)$를 활용하여 KM estimator $\hat{S}(t)$를 보정할 수 있다는 idea

  - $S^C(t)\text{ (corrected survival function)} = \hat{S}(t)/S^*(t)$

    $S^*(t) = \exp{[-\Theta(t)]}$

  - 실제 모집단의 cumulative hazard가 주어져있으니, 그 정보를 활용한다는 의미

**Notes**

- Mixed excess mortality model도 활용할 수 있음 ($h_j=\alpha(t) + \beta(t)\theta_j(t)$)

  이 model은 additive regression formulation (chapter 10)으로 fitting 가능함



#### 6.4 Bayesian Nonparametric Methods

**Ideas**

- survival function $S(t)$는 Bayesian nonparametric method을 사용하여 추정할수도 있음
- survival function에 대한 prior knowlegde를 parameter로 표현하고, 관측된 데이터로 이를 보정한다는 컨셉임
- 크게 두 가지 접근법이 있는데,
  ① conjugate prior 사용하여 analytic solution을 도출하는 방법과
  ② Gibbs sampler를 사용하여 근사하는 방법이 있음
- conjugate prior 접근법은 다시 두 가지로 나뉘어, survival function에 부여하는 방향과 cumulative hazard에 부여하여 간접추론하는 방향이 있음

**Conjugate prior - on survival function (via Dirichlet process)**

- assumptions

  - $S(t) \sim \text{Dir}(\alpha) $

  - $\alpha([t,\infty))=cS_0(t)$

    - $S_0(t)$ : prior guess on survival function
    - $c$ : measure of weight to put on prior guess

  - $E[S(t)]=\cfrac{\alpha(t,\infty)}{\alpha(0,\infty)}=S_0(t)$

    $V[S(t)]=\cfrac{[\alpha(0,\infty)-\alpha(t,\infty)]\alpha(t,\infty)}{[\alpha^2(0,\infty)+\alpha^3(0,\infty)]}=\cfrac{S_0(t)(1-S_o(t))}{c+1}$

  - $0=t_0<t_1<\cdots<t_M<t_{M+1}=\infty$ : $M$ distinct event times
    (including censoring)

  - $d_i$ : number of death at time $t_i$

  - $Y_i$ : number at risk at time $t_i$

  - $\lambda_i$ : number of censored obs. at time $t_i$

  - $\delta_j$ : indicator of death of $j$th individual

  - $\Delta_i=\begin{cases}1,& d_i>0 \\ 0,& d_i=0\end{cases}$

- posterior estimation

  - $\displaystyle \alpha^*((a,b)) = \alpha((a,b)) + \sum_{j=1}^n I[\delta_j>0,\ a<T_j <b]$
    $\text{(posterior parameter of Dirichlet distribution)}$
  - $\displaystyle \tilde{S}_D(t)=\cfrac{\alpha(t,\infty)+Y_{i+1}}{\alpha(0,\infty)+n}\prod_{k=1}^n\cfrac{\alpha(t_k,\infty)+Y_{k+1}+\lambda_k}{\alpha(t_k,\infty)+Y_{k+1}}$
    $\text{for }t_i \le t < t_{i+1},\ \ i=0,1,\cdots,M$
  - $n \rightarrow \infty$ 이면 위의 베이즈 추정량은 KM 추정량과 같게 된다

**Conjugate prior - on cumulative hazard function (via beta process)**

- assumptions

  - $A_i=[a_{i-1},a_i),\ i=1,2,\cdots,k\text{ (non-overlapping intervals)}$
    $0=a_0<a_1<\cdots<a_k$

  - $H(a_i)-H(a_{i-1})=W_i\overset{iid}{\sim}\text{Beta}(p_i,q_i)$

    - $p_i=c\left(\cfrac{a_i+a_{i-1}}{2}\right)*[H_0(a_i)-H_0(a_{i-1})]$

      $q_i=c\left(\cfrac{a_i+a_{i-1}}{2}\right)*[1-(H_0(a_i)-H_0(a_{i-1}))]$

    - $H_0(t)$ : prior guess at $t$ of the $H(t) $

    - $c(t)$ : measure of weight to put on prior guess

  - $E(W_i)=H_0(a_i) - H_0(a_{i-1})$

  - $W_i$는 $a_i$의 간격이 매우 작을 때, $dH(s)(=h(s))$가 beta dist.를 따른다는 가정을 표현하는 식과 같다.

  - 나머지 assumption은 Dirichlet prior와 같다

- posterior estimation

  - $\displaystyle \tilde{S}_B(t)=\exp\left\{ -\sum_{k=1}^i \int_{t_{k-1}}^{t_k} \cfrac{c(u)h_0(u)}{c(u)+Y_k} -\int_{t_{i}}^t \cfrac{c(u)h_0(u)}{c(u)+Y_{i+1}}du\right\}\times\prod_{k:t_k\le t}\left[ 1-\cfrac{c(t_k)h_0(t_k)+d_k}{c(t_k) + Y_k} \right]^{\Delta_k}$
    $\text{if }t_i \le t < t_{i+1}$
  - $c(t) \rightarrow 0$ 일 때 위 추정량은 KM과 같다

**Gibbs sampler method (via Dirichlet distribution)**

- Gibbs sampler method(MCMC)는 conjugate prior를 활용하는 것보다 유연한 추정이 가능하다 (**모든 censoring / truncation scheme에 대해 적용 가능하며 regression 문제에도 적용할 수 있음**)

- assumptions

  - $0<t_1<\cdots<t_M$ : $M$ timepoints

  - $d_j$ : number of deaths in the interval $[t_{j-1},t_j)$

  - $\lambda_j$ : number of right-censored obs. at $t_j$

  - $\theta_j = S(t_{j-1}) - S(t_j)$, $\theta_{M+1}=S(t_M)$

  - $\boldsymbol{\theta}\sim\text{Dir}(\boldsymbol{\alpha})$ (multivariate Dirichlet distribution)

    - $\alpha_j=C[S_0(t_{j-1})-S_0(t_j)]\ \text{ for }j=1,\cdots,M+1$

    - $S_0(t)$ : prior guess on survival function

      $S_0(t_{M+1})=0$

    - $C()$ : weight

- sampling for $i$th iteration

  0. initiate $\theta_j$s with random sample from $W_j\overset{iid}{\sim}\text{Gamma}(\alpha_j,1)$ and evaluate $\theta_j=W_j/\sum W_k$

  1. generate $Z_{j+1,j},\cdots,Z_{M+1,j}$ if $\lambda_j>0$ from $\text{multi}(\lambda_j,\rho_{j+1},\cdots,\rho_{M+1})$

     where $\rho_k=\cfrac{\theta_k^i}{\sum_{h=j+1}^{M+1}\theta_h^i}$

  2. update $R_{h}^{i+1}=\alpha_h+d_h+\sum_{j=1}^MZ_{h,j}$ (posterior Dirichlet parameter)

  3. sample $\boldsymbol{\theta}^{i+1}=(\theta_1^{i+1},\cdots,\theta_{M+1}^{i+1})$ from $\text{Dir}(\mathbf{R}^{i+1})$

     where $\mathbf{R}^{i+1}=(R_1^{i+1},\cdots,R_{M+1}^{i+1})$

  4. repeat 1~3 for $i$ times ($i$ : relatively small order such as 10 or 20) and store $\mathbf{R} $

- posterior estimation

  - $S$ : number of generated Gibbs sample ($\mathbf{R}_{s}$,  $s=1,2,\cdots,S$)
  - $\displaystyle \tilde{\theta}_h=\cfrac{1}{S}\sum_{s=1}^S\cfrac{R_{hs}^i}{\sum_{k=1}^{M+1}R_{ks}^i}$ (posterior estimator of magnitude (jump) of survival function)
  - can employ $\alpha, 1-\alpha$th quantile as $1-\alpha$ CI for $\tilde{\theta} $

