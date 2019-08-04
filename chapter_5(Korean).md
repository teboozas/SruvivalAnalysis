## 5. Estimation of Basic Quantities for Other Sampling Schemes

#### 5.1 Introduction

- "right-censored + left-truncated" 상황 외의 다른 sampling schemes에 대한 추정 방법을 간단히 살펴보자

#### 5.2 Estimation of the Survival Function for Left, Double, and Interval Censoring

**Ideas**

- $\hat{S}(t)$를 left, double, interval censored data에 적용해보자
- 각각의 sampling scheme에 따라 survival function의 construction이 달라짐

**Left-censoring scheme**

- time scale을 reverse시켜서 해결함 (Ware and Demets, 1976)

- 추정의 대상이 $P[\tau-X>t] = P[X<\tau-t]$가 됨

  (그래서 left-censoring의 경우 likelihood가 $1-S(C_l)$로 fitting 됨, quick note 3.5)

**Double-censoring scheme (both right and left)**

- modified Product-Limit estimator를 사용 (Turnbull, 1974)

- 소위 **Turnbull algorithm**이라는 알고리즘에 의해 근사됨 (closed-form이 없음)

- 기본 요소

  - $0=t_0<t_1<t_2<\cdots<t_m$ : grid of timepoints
  - $d_i$ : number of death at time $t_i$, $d_i=0$일 수 있음!
  - $r_i$ : number of right-censored individuals at time $t_i$
  - $c_i$ : number of left-censored individuals at time $t_i$
    - 이 때 실제 이벤트 발생 시간을 $t_j\ (\le t_i)$라 가정함

- **Turnbull's algorithm**

  1. left-censoring이 발생한 모든 $t_j$에 대해 초기 추정치 $S_o(t_j)$를 추정한다. 이 때 권장 추정량은 left-censored obs.를 무시하고 계산한 Product-Limit est.이다.

  2. 계산된 $\hat{S}_K$를 활용해 $\hat{p}_{ij}=\cfrac{S_K(t_{j-1})-S_K(t_j)}{1-S_K(t_j)}$를 계산한다.

     ($p_{ij}=P[t_{j-1}< X \le t_i\ | \ X \le t_i]$)

  3. 계산된 $\hat{p}_{ij}$를 활용해 $\hat{d}_{i}=d_i+\sum_{i=j}^m c_ip_{ij}$를 업데이트 한다.

  4. 계산된 $\hat{d}_i$를 활용해 일반적인 Product-Limit estimator를 계산한다. 이 값이 이전에 계산한 Product-Limit estimator와 모든 $t_i$에서 거의 같다면 업데이트를 중단하고 이 값을 사용한다. 만약 그렇지 않다면 2-4를 반복한다.

**Interval-censoring scheme**

- modified Turnbull's algorithm을 사용 (Turnbull, 1976)
- 기본 요소
  - $(L_i,R_i]$ : censoring interval of $i$th individual
  - $0=\tau_0<\tau_1 < \cdots < \tau_m$ : grid of timepoints including all $L_i,R_i$
  - $\alpha_{ij}=\begin{cases} 1,& (\tau_{j-1},\tau_j] \sub (L_i,R_i] \\ 0, & \text{o/w} \end{cases}$
- **modified Turnbull's algorithm**
  1. initial guess $S(\tau_j)$를 계산한다
  2. 계산된 $S$를 활용해 $p_j = S(\tau_{j-1}) - S(\tau_j)$를 계산한다
  3. 계산된 $p_j$를 활용해 $\displaystyle d_i=\sum_{i=1}^n \cfrac{\alpha_{ij}p_j}{\sum_k \alpha_{ik}p_k}$를 계산한다
  4. 계산된 $d_i$를 활용해 $t_i = \sum_{k=j}^m d_k$를 계산한다
  5. 계산된 pseudo data $d_i, Y_i$를 활용해 일반적인 Product-Limit estimator를 계산하고 그 값을 비교하여 업데이트 한다. 만약 값이 비슷하지 않으면 2-5를 반복한다.



#### 5.3 Estimation of the Survival Function for Right-Truncated data

**Ideas**

- censoring과 더불어 right-truncation이 발생했을 때 handling 방법을 알아보자
- 감염질환(infectious diseases) data에서 자주 발생하는 truncation 유형

**Right-truncation scheme**

- 시간축을 reverse시켜서 추정한다
- 기본 요소
  - $T_i$ : 감염 시간
  - $X_i$ : (감염 ~ 발병) 사이의 시간
  - $(0,\tau)$ : 환자가 sampling 되는 기간 → $\tau$시점 전까지 이미 감염된 사람만 포함된 data (right-truncation)
  - $R_i = \tau - X_i$ : left-truncation time (왜냐면 $T_i \le R_i$인 개체들만 포함됨)
- 방법: left-truncation method (조건부 Product-Limit est. 또는 limited Nelson-Aalen est.)를 사용한다



#### 5.4 Estimation of Survival in the Cohort Life Table

**Ideas**

- cohort: 특정 시점으로부터 동시에 출발하여 event time의 관찰이 진행되는 동일 집단
- 인구 사망 추계 데이터등이 대표적이며, 연구 특성상 censoring이 발생하기 쉬움
- 코호트 연구를 기반으로 한 생명표를 **cohort life table**이라 함

**Basic construction of Cohort life table**

- cohort life table은 기본적으로 **11개의 column**으로 구성됨
- **1st column**
  - $I_j = (a_{j-1},a_j],\ j=1,\cdots,k+1\ \text{(interval of study)}$
    - $a_j$ : interval을 나누는 경계값, $a_0 = 0,\ a_{k+1} = \infty$
  - event 또는 censoring은 반드시 저 구간들 중 하나에만 포함됨
    (미포함 또는 중복포함 없음)
- **2nd column**
  - $Y_j^{'}\ \text{(who have not experienced the event, entering }j\text{th interval)}$
- **3rd column**
  - $W_j \ \text{(number of censored(lost-to-follow-up) individuals in }j\text{th interval)}$
- **4th column**
  - $Y_j = Y_j^{'}-W_j/2\ \text{(estimator of true }Y_j)$
- **5th column**
  - $d_j\ \text{(number of events in }j\text{th interval)}$
- **6th column**
  - $\hat{S}(a_{j-1})\ \text{(estimate of survival function)}$
  - $\hat{S}(a_j) = \hat{S}(a_{j-1})[1-d_j/Y_j]$
  - $\hat{S}(a_0)=1$
- **7th column**
  - $\hat{f}(a_{mj}) \ \text{(estimate of pdf at }a_{mj})$
  - $a_{mj} = (a_j + a_{j-1})/2$ : midpoint of $j$th interval
  - $\hat{f}(a_{mj})=[\hat{S}(a_{j-1})-\hat{S}(a_{j})]/(a_j-a_{j-1})$
- **8th column**
  - $\hat{h}(a_{mj})\ \text{(estimate of hazard rate)}$
  - $\hat{h}(a_{mj}) = \hat{f}(a_{mj})/\hat{S}(a_{mj}) = \cfrac{2\hat{f}(a_{mj})}{[\hat{S}(a_j) + \hat{S}(a_{j-1})]}$ (첫 번째 방법)
  - $\hat{h}(a_{mj}) =  \cfrac{d_j}{(a_j - a_{j-1})(Y_j - d_j/2)}$ (두 번째 방법)
  - $\hat{h}(a_{mj}) = \cfrac{2\hat{q}_j}{(a_j - a_{j-1})(1+\hat{p}_j)},\ \left( \hat{q}_j = d_j/Y_j,\ \ \  \hat{p}_j=1-\hat{q}_j \right)$ (세 번째 방법)
- **9th column**
  - $s.e.(\hat{S}(a_{j-1}))\ \text{(standard error of survival function)}$
  - $\displaystyle s.e.(\hat{S}(a_{j-1})) = \hat{S}(a_{j-1})\sqrt{\sum_{i=1}^{j-1}\cfrac{d_i}{Y_i(Y_i - d_i)}}$
  - $s.e.(\hat{S}(a_{0})) = 0$ (constant)
- **10th column**
  - $s.e.(\hat{f}(a_{mj}))\ \text{(standard error of pdf)}$
  - $\displaystyle s.e.(\hat{f}(a_{mj})) = \cfrac{\hat{S}(a_{j-1})\hat{q}_j}{(a_j - a_{j-1})} \sqrt{\sum_{i=1}^{j-1}[\hat{q}_i/(Y_i\hat{p_i})] + [\hat{p}_j/(y_j\hat{p}_j)]}$
- **11th column**
  - $s.e.(\hat{h}(a_{mj}))\ \text{(standard error of hazard rate)}$
  - $s.e.(\hat{h}(a_{mj})) = \left\{\cfrac{1-[\hat{h}(a_{mj})(a_j - a_{j-1})/2]}{Y_jq_j}\right\}^{1/2}\cdot\hat{h}(a_{mj})$
- **summary**
  - 1~2 : 구간정의 및 생존자 숫자 counting
  - 3~5 : censoring 여부 체크 + 보정된 생존자수와 event 발생 수 counting
  - 6~8 : 구간별 $S, f, h$ 추정
  - 9~11 : 구간별 $S,f,h$의 표준오차 추정

**Cohort life table의 활용**

- **median life time 추정**
  - $\hat{x}_{0.5} = a_{j-1} + [\hat{S}(a_{j-1})-0.5]/\hat{f}(a_{mj})$
  - censoring이 발생한 경우 그 유형에 따라 $\hat{S}(t)$를 추정하는 방법이 달라서 mean보다는 median 추정을 선호함
- **median residual life time 추정**
  - $\hat{\text{mdrl}}(a_{i-1})=(a_{j-1} - a_{i-1}) + \cfrac{[\hat{S}(a_{j-1})-\hat{S}(a_{i-1})/2](a_j-a_{j-1})}{\hat{S}(a_{j-1}) - \hat{S}(a_j)}$
  - $a_{i-1}$ : 특정 개체들의 현재까지의 생존 시간

