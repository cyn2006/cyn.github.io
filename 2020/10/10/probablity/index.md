[$\color{purple}{\text{back to top}}$](https://cyn2006.github.io)

萌新第一次学概率期望。

如无特殊说明，以下题目中取一样本的概率相等。

<div>
    <font size="6" style="font-family:'Trebuchet MS','Lucida Sans Unicode','Lucida Grande','Lucida Sans',Arial,sans-serif">
        Problem-Set
    </font>
	</br>
	<font size="5" style="font-family:'Trebuchet MS','Lucida Sans Unicode','Lucida Grande','Lucida Sans',Arial,sans-serif">
        例 1
    </font>
</div>

有一个硬币，每次可以投正或反面，当投到连续的三个面为“正反正”时停止投币，求期望的投币次数。

<div>
	<font size="5" style="font-family:'Trebuchet MS','Lucida Sans Unicode','Lucida Grande','Lucida Sans',Arial,sans-serif">
        Solution
    </font>
</div>

我们设目前投币状态为 $f(0\ldots3)$：

- $f(0)$：当前满足终止条件（下简称条件）的投币集合为空集。
- $f(1\ldots 3)$：当前满足条件的投币集合大小分别为 $1\ldots 3$。

<div>
    <font size="4" style="font-family:'Trebuchet MS','Lucida Sans Unicode','Lucida Grande','Lucida Sans',Arial,sans-serif">f(3)</font>
</div>

现在我们考虑对于一个最终状态 $f(3)$，它的投币期望次数为 $0$。

<div>
    <font size="4" style="font-family:'Trebuchet MS','Lucida Sans Unicode','Lucida Grande','Lucida Sans',Arial,sans-serif">f(2)</font>
</div>

然后对于 $f(2)$，它投一次币后：

- 如果投的是正面，那么可以到“正反正”；
- 反之可以到空集。

考虑到每种状态的转移是等概率的，然后逆向地推：$f(2)=\dfrac{1}{2} f(0)+\dfrac{1}{2} f(3)+1$（$+1$ 是因为期望次数由之前状态转以后应增加）。

<div>
    <font size="4" style="font-family:'Trebuchet MS','Lucida Sans Unicode','Lucida Grande','Lucida Sans',Arial,sans-serif">f(1)</font>
</div>

然后对于 $f(1)$，它投一次币后：

- 如果投的是正面，仍然是“正”；
- 反之可以到状态“正反”。

所以 $f(1)=\dfrac{1}{2} f(2)+\dfrac{1}{2} f(1)+1$。

<div>
    <font size="4" style="font-family:'Trebuchet MS','Lucida Sans Unicode','Lucida Grande','Lucida Sans',Arial,sans-serif">f(0)</font>
</div>

然后对于 $f(0)$，它投一次币后：

- 如果投的是正面，可以到状态“正”；
- 反之仍然是空集。

所以 $f(0)=\dfrac{1}{2} f(0)+\dfrac{1}{2}f(1)+1$。

<div>
    <font size="4" style="font-family:'Trebuchet MS','Lucida Sans Unicode','Lucida Grande','Lucida Sans',Arial,sans-serif">Solving</font>
</div>

开始~~毒奶zzc初赛98分~~解方程了：
$$
\begin{cases} f(0) = \dfrac{1}{2} f(0) + \dfrac{1}{2} f(1) + 1 \\ f(1) = \dfrac{1}{2} f(1) + \dfrac{1}{2} f(2) + 1 \\ f(2) = \dfrac{1}{2} f(0) + \dfrac{1}{2} f(3) + 1 \\ f(3)=0 \end{cases} \Longrightarrow f(2) =5\Longrightarrow f(0)=10
$$
所以从状态空集： $0$ 开始的期望终止次数为 $10$。

<div>
    <font size="5" style="font-family:'Trebuchet MS','Lucida Sans Unicode','Lucida Grande','Lucida Sans',Arial,sans-serif">
        例 2
    </font>
</div>

假设现在有无数个红球、蓝球和黄球，现在开始取球，累计取到一个红球或者两个黄球则停止，问期望取球的数量。

<div>
	<font size="5" style="font-family:'Trebuchet MS','Lucida Sans Unicode','Lucida Grande','Lucida Sans',Arial,sans-serif">
        Solution
    </font>
</div>

我们设：

- $f(0)$ 表示目标状态只取了蓝球的期望数量；
- $f(1)$ 表示目标状态取了一个黄球和若干蓝球的期望数量；
- $f(2)$ 表示目标状态终止了的期望数量。

对于 $f(2)$：因为取 $1$ 个红球和 $2$ 个黄球等价，故归到一类，共三类。

那么
$$
\begin{cases}
f(0) = \dfrac{1}{3} f(2) + \dfrac{1}{3} f(1) + \dfrac{1}{3} f(0) + 1\\
f(1) = \dfrac{2}{3} f(2) + \dfrac{1}{3} f(1) + 1\\
f(2) = 0
\end{cases}
$$
解得：$f(0)=\dfrac{9}{4}$。