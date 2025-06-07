# 拟Newton法

> 学完这里后，需要能够意识到，很多在跟黑塞矩阵近似的迭代式满足$$H_kd_k+\nabla f(x_k)=0 $$这个迭代式很有可能就是类似拟牛顿法的方法。需要直觉上考察其与黑塞矩阵的全局或局部的接近程度。

## 思想

拟Newton法希望保持快速收敛，但是不想要每次计算黑塞矩阵还要求它正定。需要构造一个对称正定矩阵近似黑塞矩阵。

$$B_k\approx \nabla^2 f(x_k),\text{or} \quad H_k\approx \nabla^2 f(x_k)^{-1}$$

> 一般$B_k$用在理论分析多，$H_k$用在实际计算多。

## 计算方法

令$$B_{k+1}s_k=y_k$$，其中$s_k=x_{k+1}-x_k,\quad y_k=\nabla f(x_{k+1})-\nabla f(x_k) $

这称为拟Newton方程，或称割线方程。

对于$H_{k+1}$，等价地写$$H_{k+1}y_k=s_k$$

算出来$B_{k+1},H_{k+1} $后代替Newton法中的黑塞矩阵即可。

## 收敛性与Dennis-Moré条件

拟Newton法能保持**超线性收敛**。

**定理**

$f$二次连续可微，对于拟牛顿法生成的$d_k$，若满足$\nabla f(x^*)=0,\quad \nabla^2 f(x^*) $正定(仅在收敛点上要求正定了)，古典迭代$$x_{k+1}=x_k+d_k $$超线性收敛于$x^*$当且仅当满足$$\lim_{k\to\infty}\frac{\Vert (B_k-\nabla^2 f(x^*))d_k\Vert}{\Vert d_k\Vert}=0 $$，又可写作$$ \lim_{k\to\infty}\frac{\Vert (B_k-\nabla^2 f(x^*))s_k\Vert}{\Vert s_k\Vert}=0$$

上面的条件称为Dennis-Moré条件。

该定理与之前的[超线性收敛定理](10-收敛速度.md/#定理2(超线性收敛定理))是一致的。


> 证明：
>
> 充分性：
>
> 首先回忆[收敛速度的定理2](10-收敛速度.md#定理2(超线性收敛定理))，可知 $$\lim_{k\to\infty}\frac{\Vert \nabla f(x_k)+\nabla^2f(x_k)d_k\Vert}{\Vert d_k\Vert}=0 $$满足时就超线性收敛。这个条件也被称为Dennis-Moré条件。
>
> 然后只需要证明$\Vert \nabla f(x_k)+\nabla^2f(x_k)d_k\Vert\leq \Vert (B_k-\nabla^2 f(x^*))d_k\Vert+o(\Vert d_k\Vert )$
> 必要性：
>
> 利用超线性收敛的定义，利用三角不等式，推出$$\lim_{k\to \infty}\frac{\Vert x_k-x^*\Vert}{\Vert d_k\Vert}=1 $$
>
> **这个不等式经常用到**
>
> 然后再依照条件推出$$ (B_k-\nabla^2 f(x^*))d_k\leq o(\Vert d_k\Vert) $$即可

# 三类重要的修正公式

不用计算黑塞矩阵一定程度上是很好的。但是求解方程组仍然是一个计算量较大的问题。

希望能够利用低秩对称矩阵来构造$B_k$的递推公式，来避免解方程组。$$B_{k+1}=B_k+\Delta_k $$

## SR1修正公式

取$\Delta_k$的秩为1，那么就有$$\Delta_k=\beta_ku_ku_k^T ,\beta_k\in R$$，将其代入方程整理，有$$\beta_ku_ku_k^Ts_k=y_k-B_ks_k $$

发现后面的$u_k^Ts_k$为标量，就可以知道$u_k\Vert y_k-B_ks_k$，也就是共线。

不妨直接取$u_k=y_k-B_ks_k$，那么也就能直接得到$\beta_k$的表达式。

最终就得到秩1修正公式$$B_{k+1}=B_k+\frac{(y_k-B_ks_k)(y_k-B_ks_k)^T}{(y_k-B_ks_k)^Ts_k} $$

其实，将$y_k$与$s_k$交换，$B$改为$H$就有$H$的表达式。

1. 这个方法不需要搜索，可以直接取$\alpha_k=1$，并且具有二次终止性。
2. 遗传性质，保持$B_ks_i=y_i,\forall i<k $，也就是自对偶性。
3. 但是它们无法保证正定性。

## BFGS修正公式

现在取其秩为2，那么$\Delta_k=a_ku_ku_k^T+b_kv_kv_k^T$

代入有$$B_{k+1}=B_k+a_ku_ku_k^T+b_kv_kv_k^T $$

就有$$ (a_ku_ku_k^T+b_kv_kv_k^T)s_k=y_k-B_ks_k $$

上式没有唯一的取法，但是可以简单地一对一取，也就是令$u_k=y_k,\quad v_k=B_ks_k$

这个时候就有$$B_{k+1}=B_k-\frac{B_ks_ks_k^TB_k}{s_k^TB_ks_k}+\frac{y_ky_k^T}{y_k^Ts_k} $$

BFGS公式是迄今最好的拟Newton修正公式。
1. 精确搜索下具有二次终止性
2. 遗传性质，保持$B_ks_i=y_i,\forall i<k $，也就是自对偶性。
3. 正定性保证，容易成立的$y_k^Ts_k>0$条件下就能保证$B_{k+1}$对称正定。

### 如何保证$y_k^Ts_k>0$

只要下列条件之一成立，则$y_k^Ts_k>0$成立。

1. 算法采用Wolfe-Powell搜索或精确搜索
2. 函数$f$二次连续可微，$\forall x\in R^n,\nabla^2 f(x) $正定

> 为了改进其在Armijo搜索与非凸函数上缺乏保证的正定性，还有一些改进形式被提出。

> $B_0$如何选取？可以选单位矩阵，也可以尝试在尺度上接近黑塞矩阵。无论如何，简单而接近即可。
