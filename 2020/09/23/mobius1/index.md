[$\color{purple}{\text{back to top}}$](https://cyn2006.github.io)

# 莫反学习笔记 $1$

我太菜了。

我是菜鸡。

持续更新。

以下为紫题集合

参考 $\text{OI-Wiki}$。

## $\mathcal {Lemma}$

**$1.1$**, $\forall a,b,c\in \mathbb Z, \left\lfloor \dfrac{\left\lfloor \dfrac{a}{b}\right\rfloor}{c}\right\rfloor=\left\lfloor\dfrac{a}{bc}\right\rfloor$.

**$2.1$**, $\epsilon(n)=[n=1]$.

**$2.2$**, $\operatorname{id}_{k}(n)=n^k \operatorname{id}_{1}(n)$，简记作 $\operatorname{id}(n)$.

**$2.3$**, $1(n)=1$.

**$2.4$**, $\varphi(n)=\sum\limits_{i=1}^n[\gcd(i,n)=1]$

**$2.5$**, $\mu(n)=\begin{cases} 1\ (n=1)\\0\ (\exists d>1:d^2|n)\\(-1)^{\omega(n)}\ (otherwise)\end{cases}$

**$3.1$**, $\sum\limits_{i=1}^n\sum\limits_{j=1}^m [\gcd(i,j)=1]=\sum\limits_{i=1}^n\sum\limits_{j=1}^m \sum\limits_{d|i,d|j} \mu(d)$，特别地，对于 $m=i$，等价于 $\sum\limits_{i=1}^n \varphi(i)$。

---------

## P3455 [POI2007]ZAP-Queries

> Problem b 的弱化版。

考虑直接前缀即可。

## P2522 [HAOI2011]Problem b

形式化地，差分之后需要求出 $\sum\limits_{i=1}^n\sum\limits_{j=1}^n [\gcd(i,j)=k]$。

我们令 $S=\sum\limits_{i=1}^n\sum\limits_{j=1}^n [\gcd(i,j)=k]$，那么

$$A=\sum_{i=1}^{\lfloor \frac{n}{k} \rfloor}\sum_{j=1}^{\lfloor \frac{m}{k} \rfloor} [\gcd(i,j)=1]$$

$\because \varepsilon(n)=\sum_{d|n} \mu(d)$

$$\therefore A=\sum_{i=1}^{\lfloor \frac{n}{k} \rfloor}\sum_{j=1}^{\lfloor \frac{m}k{} \rfloor} \sum_{d|\gcd(i,j)} \mu(d)$$

$$A=\sum_{d=1}^{\lfloor \frac{n}{k} \rfloor} \mu(d) \sum_{i=1}^{\lfloor \frac{n}{k} \rfloor} [d|i] \sum_{i=1}^{\lfloor \frac{m}{k} \rfloor} [d|j]$$

$$A=\sum_{d=1}^{\lfloor \frac{n}{k} \rfloor} \mu(d) \lfloor\frac{n}{kd}\rfloor \lfloor \frac{m}{kd} \rfloor$$

于是整除分块即可。

$Code:$

```cpp
#include<bits/stdc++.h>

#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;++i)
#define ll long long
#define N 100005

int mu[N],pri[N],tot,n;
bool f[N];
inline void init(int n){
	mu[1]=1;
	rep(i,2,n){
		if(!f[i]) pri[++tot]=i,mu[i]=-1;
		for(int j=1;j<=tot&&i*pri[j]<=n;++j){
			f[i*pri[j]]=1;
			if(i%pri[j]==0){mu[i*pri[j]]=0; break;}
			mu[i*pri[j]]=-mu[i];
		}
	}
	rep(i,1,n) mu[i]+=mu[i-1];
}

inline ll solve(int x,int y,ll ans=0){
	for(int l=1,r;l<=std::min(x,y);l=r+1){
		r=std::min(x/(x/l),y/(y/l));
		ans+=1ll*(mu[r]-mu[l-1])*(x/l)*(y/l);
	}
	return ans;
}

int main(){
	std::ios::sync_with_stdio(0);
	init(1e5);
	std::cin>>n;
	while(n--){
		int a,b,c,d,k; std::cin>>a>>b>>c>>d>>k;
		ll ans=solve(b/k,d/k)-solve(b/k,(c-1)/k)-solve((a-1)/k,d/k)+solve((a-1)/k,(c-1)/k);
		std::cout<<ans<<'\n';
	}
	return 0;
}
```

----------------

## P2257 YY的GCD

$$
\large S(n,m)=\sum_{p\in \mathcal{Prime}} \sum_{i=1}^{\lfloor \frac{n}{p}\rfloor} \sum_{j=1}^{\lfloor \frac{m}{p} \rfloor} [\gcd(i,j)=1]
$$

$$
\large =\sum_{p\in \mathcal{Prime}} \sum_{i=1}^{\lfloor \frac{n}{p}\rfloor} \sum_{j=1}^{\lfloor \frac{m}{p} \rfloor}  \sum_{d|\gcd(i,j)} \mu(d)
$$

$$
\large =\sum_{p\in \mathcal{Prime}} \sum_{d=1}^{\lfloor \frac{n}{p} \rfloor} \mu(d) \lfloor \frac{n}{pd}\rfloor \lfloor \frac{m}{pd}\rfloor
$$

然后我们先枚举 $pd$ 再枚举 $p$，并且交换整除分块与 $\mu$ 的求和顺序，即得

$$
\large S(n,m)=\sum_{pd=t=1}^n \lfloor \frac{n}{t}\rfloor \lfloor \frac{m}{t}\rfloor \sum_{p|t,p\in \mathcal{Prime}} \mu(\frac{t}{p})
$$

$Code:$

```cpp
#include<bits/stdc++.h>

#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;++i)
#define ll long long
#define N 100005

int mu[N],pri[N],tot,n;
ll sum[N];
bool f[N];
inline void init(int n){
	mu[1]=1;
	rep(i,2,n){
		if(!f[i]) pri[++tot]=i,mu[i]=-1;
		for(int j=1;j<=tot&&i*pri[j]<=n;++j){
			f[i*pri[j]]=1;
			if(i%pri[j]==0){mu[i*pri[j]]=0; break;}
			mu[i*pri[j]]=-mu[i];
		}
	}
	rep(i,1,tot) for(int j=pri[i];j<=n;j+=pri[i])
		sum[j]+=mu[j/pri[i]];
	rep(i,1,n) sum[i]+=sum[i-1];
}

inline ll solve(int x,int y,ll ans=0){
	for(int l=1,r;l<=std::min(x,y);l=r+1){
		r=std::min(x/(x/l),y/(y/l));
		ans+=1ll*(sum[r]-sum[l-1])*(x/l)*(y/l);
	}
	return ans;
}

signed main(){
	std::ios::sync_with_stdio(0);
	init(5e4);
	std::cin>>n;
	while(n--){
		int a,b; std::cin>>a>>b;
		std::cout<<solve(a,b)<<'\n';
	}
	return 0;
}
```

-----------------------

## P3312 [SDOI2014]数表

> 这里我们认为 $\sigma_1(x)$ 为 $x$ 的约数的和，即 $\sigma_1(x)=\sum\limits_{d|x} d$。

首先考虑没有~~duliu~~ $a$ 的限制，我们令 $S(n,m)=$
$$
S(n,m)=\sum_{i=1}^n \sum_{j=1}^m \sigma_1 (\gcd(i,j))
$$

$$
=\sum_{d=1}^n\sum_{i=1}^{n}\sum_{j=1}^{m} \sigma_1(d)[\gcd(i,j)=d]\quad (*)
$$

$$
=\sum_{d=1}^n\sum_{i=1}^{\lfloor \frac{n}{d}\rfloor}\sum_{j=1}^{\lfloor \frac{m}{d} \rfloor} \sigma_1(d)[\gcd(i,j)=1]
$$

$$
=\sum_{d=1}^n\sigma_1 (d) \sum_{x=1}^{\lfloor \frac{n}{d}\rfloor} \sum_{i=1}^{\lfloor \frac{n}{d}\rfloor}[x|i]\sum_{j=1}^{\lfloor \frac{n}{d}\rfloor}[x|j] \times \mu(x)
$$

$$
=\sum_{d=1}^n\sigma_1 (d) \sum_{x=1}^{\lfloor \frac{n}{d}\rfloor} \mu(x)\lfloor \frac{n}{dx} \rfloor \lfloor \frac{m}{dx} \rfloor
$$

> 对于 $(*)$：考虑枚举一个 $gcd$ 为 $d$，转化成常见的 $\mathcal{Mobius}$ 反演的形式。

令 $t=dx$，则有
$$
S(n,m)=\sum_{t=1}^n \lfloor \frac{n}{t} \rfloor \lfloor \frac{m}{t} \rfloor \sum_{d|t} \sigma_1(d) \mu(\lfloor \frac{t}{d} \rfloor)
$$
考虑到前面的柿子可以整除分块，则对于后面的柿子需要进行 $\leqslant a$ 的特殊处理。

我们令后面的和式为 $F$。

考虑到调和级数为 $n\ln n$，于是对于每一个 $a$，我们将它离线下来，按从小到大枚举，然后对于每一个**约数和**在 $(a_{lst},a_{now}]$ 的数，我们可以直接枚举 $kd$（即为 $d$ 的 $k$ 倍），在树状数组上单点更新对于 $F$ 的贡献，然后对于每一个新的 $(d,\lfloor \frac{t}{d} \rfloor)$，贡献加入到 $t$ 的 $F$ 中去。

最后整除分块时候用树状数组查询区间和即可。

$Code:$
```cpp
#include<bits/stdc++.h>

#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;++i)
#define ll long long
#define N 100005
int mu[N],pri[N],tot;
ll sum[N],d[N];
bool f[N];
std::pair<ll,int> p[N];
inline void init(int n){
	mu[1]=1;
	rep(i,2,n){
		if(!f[i]) pri[++tot]=i,mu[i]=-1;
		for(int j=1;j<=tot&&i*pri[j]<=n;++j){
			f[i*pri[j]]=1;
			if(i%pri[j]==0){mu[i*pri[j]]=0; break;}
			mu[i*pri[j]]=-mu[i];
		}
	}
	rep(i,1,n) for(int j=i;j<=n;j+=i) d[j]+=i;
	rep(i,1,n) p[i]=std::make_pair(d[i],i);
	std::sort(p+1,p+n+1);
}
struct tree{
	ll c[N];
	inline void upd(int x,ll v){for(;x<=1e5;x+=x&-x) c[x]+=v;}
	inline ll qry(int x,ll ret=0){for(;x;x-=x&-x) ret+=c[x]; return ret;}
	inline ll ask(int x,int y,ll ret=0){return qry(y)-qry(x-1);}
} T;
int q;
struct data{
	int n,m,x,id;
} a[N];

inline ll solve(int x,int y,ll ret=0){
	for(int l=1,r;l<=std::min(x,y);l=r+1){
		r=std::min(x/(x/l),y/(y/l));
		ret+=T.ask(l,r)*(x/l)*(y/l);
	}
	return ret;
}

inline bool cmp(data x,data y){
	return x.x<y.x;
}
ll ans[N];
signed main(){
	std::ios::sync_with_stdio(0);
	init(1e5);
	std::cin>>q;
	rep(i,1,q) std::cin>>a[i].n>>a[i].m>>a[i].x,a[i].id=i;
	std::sort(a+1,a+q+1,cmp);
	int head=1;
	rep(i,1,q){
		while(head<=1e5&&p[head].first<=a[i].x){
			for(int v=p[head].second,j=1;j*v<=1e5;++j)
				T.upd(j*v,mu[j]*p[head].first);
			++head;
		}
//		std::cerr<<i<<' '<<head<<'\n';
		ans[a[i].id]=solve(a[i].n,a[i].m);
	}
	const ll mod=(1ll<<31)-1;
	rep(i,1,q) std::cout<<(ans[i]&mod)<<'\n';
	return 0;
}
```

----------------

## P1390 公约数的和

本题有三种解法。

设所求答案为 $S(n)$。

### 方法一

$$
S(n)=\sum_{i=1}^n \sum_{j=1}^m \gcd(i,j)
$$

$$
\sum_{i=1}^n \sum_{j=1}^m \sum_{d=1}^n d[\gcd(i,j)=d]
$$

$$
\sum_{d=1}^nd \sum_{i=1}^n\sum_{j=1}^m[\gcd(i,j)=d]
$$

$$
\sum_{d=1}^nd \sum_{i=1}^{\lfloor \frac{n}{d} \rfloor}\sum_{j=1}^{\lfloor \frac{m}{d} \rfloor} [\gcd(i,j)=1]
$$

$$
\sum_{d=1}^nd \sum_{i=1}^{\lfloor \frac{n}{d} \rfloor}\sum_{j=1}^{\lfloor \frac{m}{d} \rfloor} \sum_{x\mid i,x\mid j} \mu(x)
$$

$$
\sum_{d=1}^nd \sum_{x=1}^{\lfloor \frac{n}{d} \rfloor}\mu(x)\sum_{i=1}^{\lfloor \frac{n}{d} \rfloor}[x\mid i] \sum_{j=1}^{\lfloor \frac{m}{d} \rfloor}[x\mid j]
$$

$$
\sum_{d=1}^n d\sum_{x=1}^{\lfloor \frac{n}{d} \rfloor} \mu(x)\lfloor \frac{n}{dx} \rfloor \lfloor \frac{m}{dx} \rfloor
$$

直接枚举 $d$，然后整除分块即可。

时间复杂度：$O(n\sqrt{\ln n})$。

时间复杂度证明：

首先我们设渐进时间复杂度为 $\mathcal O$。

那么 $\mathcal O=\sum\limits_{d=1}^n \sqrt{\dfrac{n}{d}}\leqslant n \sqrt{\dfrac{\sum\limits_{i=1}^n \dfrac{n}{d}}{n}}=n\sqrt{\dfrac{n\ln n}{n}}=n\sqrt{\ln n}$。

证毕。

$Code:$

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define ll long long
#define N 2000005
int pri[N],cnt;
bool f[N];
ll phi[N];
inline void init(int n){
	phi[1]=1;
	rep(i,2,n){
		if(!f[i]) pri[++cnt]=i,phi[i]=i-1;
		for(int j=1;j<=cnt&&i*pri[j]<=n;++j){
			f[i*pri[j]]=1;
			if(i%pri[j]==0){phi[i*pri[j]]=phi[i]*pri[j]; break;}
			phi[i*pri[j]]=phi[i]*(pri[j]-1);
		}
	}
	rep(i,1,n) phi[i]+=phi[i-1];
}
inline ll solve(int x,ll ret=0){
	for(int l=1,r;l<=x;l=r+1){
		r=x/(x/l);
		ret+=1ll*(phi[r]-phi[l-1])*(x/l)*(x/l);
	}
	return ret;
}
int main(){
	std::ios::sync_with_stdio(0);
	int n; std::cin>>n;
	init(n);
	ll ret=solve(n);
	std::cout<<(ret-(n*(n+1ll)>>1)>>1)<<'\n';
	return 0;
}
```

### 方法二

然后我们将 $dx$ 提前（$\color{black}{\texttt{z}}\color{red}{\texttt{houakngyang}}$ 教我的方法），枚举 $t=dx$：

$$
\sum_{dx=t=1}^n \lfloor \frac{n}{t} \rfloor \lfloor \frac{m}{t} \rfloor \sum_{d\mid t} \mu(\frac{t}{d})\times d
$$

线筛 $\mu$ 这样就可以做到 $\mathcal{O} (n)$ 了（by zzc）。

### 方法三

考虑到 $\varphi(x)=\sum\limits_{i\mid x}\mu(\frac{x}{i})\times i$ 即 $\varphi =\mu * \operatorname{id}$，那么：

$$
\sum_{t=1}^n \lfloor \frac{n}{t} \rfloor \lfloor \frac{m}{t} \rfloor \varphi(t)
$$

线筛一遍即可。

线筛之后可以整除分块。

时间复杂度 $O(n+\sqrt{n})$（这样就可以处理多次询问了）。

$Code:$

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define ll long long
#define N 2000005
int pri[N],cnt;
bool f[N];
ll phi[N];
inline void init(int n){
	phi[1]=1;
	rep(i,2,n){
		if(!f[i]) pri[++cnt]=i,phi[i]=i-1;
		for(int j=1;j<=cnt&&i*pri[j]<=n;++j){
			f[i*pri[j]]=1;
			if(i%pri[j]==0){phi[i*pri[j]]=phi[i]*pri[j]; break;}
			phi[i*pri[j]]=phi[i]*(pri[j]-1);
		}
	}
	rep(i,1,n) phi[i]+=phi[i-1];
}
inline ll solve(int x,ll ret=0){
	for(int l=1,r;l<=x;l=r+1){
		r=x/(x/l);
		ret+=1ll*(phi[r]-phi[l-1])*(x/l)*(x/l);
	}
	return ret;
}
int main(){
	std::ios::sync_with_stdio(0);
	int n; std::cin>>n;
	init(n);
	ll ret=solve(n);
	std::cout<<(ret-(n*(n+1ll)>>1)>>1)<<'\n';
	return 0;
}
```

### $\text{More}$

考虑 $n$ 的范围到 $10^9$？

然后杜教筛即可。

---------

## P3327 [SDOI2015]约数个数和

$$
S(n,m)=\sum_{i=1}^n \sum_{j=1}^m d(i\times j)
$$

$$
=\sum_{i=1}^n \sum_{j=1}^m \sum_{x\mid i} \sum_{y\mid j} [\gcd(x,y)=1]\quad (\mathcal{Lemma})
$$

$$
=\sum_{x=1}^n\sum_{y=1}^m[\gcd(x,y)=1]\sum_{i=1}^n [x|i]\sum_{j=1}^m [y|j]
$$

$$
=\sum_{x=1}^n\sum_{j=1}^m[\gcd(x,y)=1]\lfloor \frac{n}{x} \rfloor \lfloor \frac{m}{y} \rfloor
$$

$$
=\sum_{x=1}^n\sum_{j=1}^m\lfloor \frac{n}{x} \rfloor \lfloor \frac{m}{y} \rfloor\sum_{d\mid x,d\mid y} \mu(d)
$$

$$
=\sum_{d=1}^n \mu(d)\sum_{x=1}^n [d|x] \lfloor \frac{n}{d} \rfloor \sum_{y=1}^m [d|y] \lfloor \frac{m}{d} \rfloor
$$

$$
=\sum_{d=1}^n \mu(d)\sum_{x=1}^{\lfloor\frac{n}{d}\rfloor}\lfloor \frac{n}{dx} \rfloor \sum_{y=1}^{\lfloor\frac{m}{d}\rfloor}\lfloor\frac{m}{dx}\rfloor
$$

令 $k=\lfloor \frac{n}{d} \rfloor$，则
$$
S(n,m)=\sum_{d=1}^n \mu(d) \sum_{x=1}^k\lfloor \frac{k}{x} \rfloor \sum_{y=1}^m \lfloor \frac{k}{y} \rfloor
$$
那么我们进行整除分块，预处理出对于每一个 $i$ 的 $f(i)=\sum\limits_{j=1}^i \lfloor \frac{i}{j} \rfloor$，那么有
$$
S(n,m)=\sum_{d=1}^n \mu(d)\times f(\lfloor \frac{n}{d} \rfloor)\times f(\lfloor \frac{m}{d} \rfloor)
$$
则对于每一次询问，整除分块即可。

时间复杂度 $O(n\sqrt{n}+T\sqrt{n})$。

$Code:$

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define ll long long
#define N 2000005
int mu[N],pri[N],cnt;
ll sum[N];
bool f[N];
inline ll solve(int x,ll ret=0){
	for(int l=1,r;l<=x;l=r+1){
		r=x/(x/l);
		ret+=(r-l+1ll)*(x/l);
	}
	return ret;
}
inline void init(int n){
	mu[1]=1;
	rep(i,2,n){
		if(!f[i]) pri[++cnt]=i,mu[i]=-1;
		for(int j=1;j<=cnt&&i*pri[j]<=n;++j){
			f[i*pri[j]]=1;
			if(i%pri[j]==0){mu[i*pri[j]]=0; break;}
			mu[i*pri[j]]=-mu[i];
		}
	}
	rep(i,1,n) mu[i]+=mu[i-1];
	rep(i,1,n) sum[i]=solve(i);
}
inline ll solve2(int x,int y,ll ret=0){
	for(int l=1,r;l<=std::min(x,y);l=r+1){
		r=std::min(x/(x/l),y/(y/l));
		ret+=1ll*(mu[r]-mu[l-1])*sum[x/l]*sum[y/l];
	}
	return ret;
}
int main(){
	std::ios::sync_with_stdio(0);
	init(5e4);
	int T; std::cin>>T;
	while(T--){
		int n,m; std::cin>>n>>m;
		std::cout<<solve2(n,m)<<'\n';
	}
	return 0;
}
```

----------------

## P1829 [国家集训队]Crash的数字表格
令
$$
S(n,m)=\sum_{i=1}^n \sum_{j=1}^m \operatorname{lcm}(i,j)
$$

则
$$
S(n,m)=\sum_{i=1}^n \sum_{j=1}^m \dfrac{i\times j}{\gcd(i,j)}
$$

$$
=\sum_{i=1}^n \sum_{j=1}^m \sum_{d\mid i,d\mid j} [\gcd(i,j)=d]\times \dfrac{i\times j}{d}
$$

$$
=\sum_{d=1}^n d\sum_{i=1}^{\lfloor \frac{n}{d} \rfloor} \sum_{j=1}^{\lfloor \frac{m}{d} \rfloor} [\gcd(i,j)=1]\times i\times j
$$

$$
=\sum_{d=1}^n d\sum_{i=1}^{\lfloor \frac{n}{d} \rfloor} \sum_{j=1}^{\lfloor \frac{m}{d} \rfloor} \sum_{x\mid i,x\mid j} \mu(x)\times i\times j
$$

$$
=\sum_{d=1}^n d\sum_{x=1}^{\lfloor \frac{n}{d} \rfloor} \mu(x) \sum_{i=1}^{\lfloor \frac{n}{d} \rfloor}[x\mid i]\times i\sum_{j=1}^{\lfloor \frac{m}{d} \rfloor} [x\mid j] \times j
$$

$$
=\sum_{d=1}^n d\sum_{x=1}^{\lfloor \frac{n}{d} \rfloor} \mu(x) \times x^2 \sum_{i=1}^{\lfloor \frac{n}{dx} \rfloor} i\sum_{j=1}^{\lfloor \frac{m}{dx} \rfloor} j
$$

记 $ g(x)=\sum\limits_{i=1}^xi=\dfrac{x(x+1)}{2}$，则有

$$
S(n,m)=\sum_{d=1}^n d\sum_{x=1}^{\lfloor \frac{n}{d} \rfloor} \mu(x) \times x^2\times g(\lfloor \frac{n}{dx} \rfloor)\times g(\lfloor \frac{m}{dx} \rfloor)
$$

然后整除分块套整除分块即可。

时间复杂度 $O(\sqrt{n}\times \sqrt{n})=O(n)$。

$Code:$

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define ll long long
#define N 10000005
const int mod=20101009;
ll mu[N];
int pri[N],cnt;
bool f[N];
inline ll C(int x,int y){return ((0ll+x+y)*(y-x+1ll)>>1)%mod;}
inline ll upd(ll x){return x+(x>>31&mod);}
inline void init(int n){
	mu[1]=1;
	rep(i,2,n){
		if(!f[i]) pri[++cnt]=i,mu[i]=-1;
		for(int j=1;j<=cnt&&i*pri[j]<=n;++j){
			f[i*pri[j]]=1;
			if(i%pri[j]==0){mu[i*pri[j]]=0; break;}
			mu[i*pri[j]]=-mu[i];
		}
	}
	rep(i,1,n) mu[i]=upd((mu[i-1]+mu[i]*i%mod*i)%mod);
}
inline ll solve2(int x,int y,ll ret=0){
	for(int l=1,r;l<=std::min(x,y);l=r+1){
		r=std::min(x/(x/l),y/(y/l));
		(ret+=upd(mu[r]-mu[l-1])*C(1,x/l)%mod*C(1,y/l))%=mod;
	}
	return ret;
}
inline ll solve1(int x,int y,ll ret=0){
	for(int l=1,r;l<=std::min(x,y);l=r+1){
		r=std::min(x/(x/l),y/(y/l));
		(ret+=solve2(x/l,y/l)*C(l,r))%=mod;
	}
	return ret;
}
int main(){
	int n,m;
	std::cin>>n>>m;
	init(std::max(n,m));
	std::cout<<solve1(n,m)<<'\n';
	return 0;
}
```

-----------------

## P4213 【模板】杜教筛（Sum）

杜教筛。

首先我们要计算出 $S(n)=\sum\limits_{i=1}^n S(i)$。

$g(x)$ 是我们需要构造出来的一个数论函数。

现在有一些性质（显然数论函数必须满足）：
$$
\sum_{i=1}^n \sum_{d\mid i} f(d)\times g(\frac{i}{d})=\sum_{i=1}^ng(i)\times S(\lfloor \frac{n}{i} \rfloor)
$$

$$
\Longleftrightarrow \sum_{i=1}^n (f*g)(i)=\sum_{i=1}^n g(i)\times S(\lfloor \frac{n}{i} \rfloor)
$$

$f(d)g(\frac{i}{d})$ 就是对于所有的 $i\leqslant n$ 所做的贡献，那么**变换枚举顺序** $d,\frac{i}{d}$ 得到：

$$
\sum_{i=1}^n \sum_{d\mid i} f(d)\times g(\frac{i}{d})
$$

$$
=\sum_{d=1}^n \sum_{x=1}^{\lfloor \frac{n}{d}\rfloor} g(d)\times f(x)
$$

$$
=\sum_{d=1}^ng(d)\sum_{x=1}^{\lfloor \frac{n}{d} \rfloor} f(x)
$$

$$
=\sum_{d=1}^ng(d)\times S(\lfloor \frac{n}{d} \rfloor)
$$

$$
=g(1)S(n)+\sum_{d=2}^n g(d)\times S(\lfloor \frac{n}{d} \rfloor)
$$
联系性质可以得到：
$$
g(1)S(n)=\sum_{i=1}^n(f*g)(i)-\sum_{d=2}^n g(d)\times S(\lfloor \frac{n}{d} \rfloor)
$$
变量统一一下得到：
$$
g(1)S(n)=\sum_{i=1}^n(f*g)(i)-\sum_{i=2}^n g(i)\times S(\lfloor \frac{n}{i} \rfloor)
$$

如果我们能对于 $RHS$ 的前面一个柿子快速求出答案，并对后面一个柿子在 $O(\sqrt{n})$ 内整除分块，这样就能较为快速地求出 $g(1)S(n)$ 了。

-----

好了正文开始。

对于 $\mathcal{Mobius}$ 函数的前缀和：
$$
\because \epsilon = \mu * 1
$$

$$
\therefore \epsilon = \sum_{d\mid n} \mu(d)
$$

$$
S(n)=\sum_{i=1}^n \epsilon (i) -\sum_{i=2}^n S(\lfloor \frac{n}{i} \rfloor)
$$

$$
=1-\sum_{i=2}^n S(\lfloor \frac{n}{i} \rfloor)
$$
对于 $\mathcal{\varphi}$ 函数的前缀和：

考虑用莫反：
$$
\sum_{i=1}^n \sum_{j=1}^n [\gcd(i,j)=1]
$$

$$
=\sum_{i=1}^n\sum_{j=1}^m \sum_{d\mid i,d\mid j} \mu(d)
$$

$$
=\sum_{d=1}^n \mu(d)\sum_{i=1}^n [d\mid i] \sum_{j=1}^n [d\mid j]
$$

$$
=\sum_{d=1}^n \mu(d)\lfloor \frac{n}{d} \rfloor^2
$$
整除分块即可。

考虑用杜教筛：

求 $S(i)=\sum_{i=1}^n\varphi(i)$。

由性质可以知道：$\varphi * 1 = \operatorname{id}$。
$$
\sum_{i=1}^n (\varphi * 1)(i)=\sum_{i=1}^n S(\lfloor \frac{n}{i} \rfloor)
$$

$$
\sum_{i=1}^n\operatorname{id}(i)=\sum_{i=1}^nS(\lfloor \frac{n}{i} \rfloor)
$$

$$
\dfrac{n(n+1)}{2}=\sum_{i=1}^nS(\lfloor \frac{n}{i} \rfloor)
$$

$$
S(n)=\dfrac{n(n+1)}{2}-\sum_{i=2}^nS(\lfloor \frac{n}{i} \rfloor)
$$

莫反：

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)

#define ll long long
#define N 5000005
int T,mx,mu[N],pri[N],cnt;
bool f[N];
std::unordered_map<ll,ll> _mu,_phi;
inline void init(int n){
	mu[1]=1;
	rep(i,2,n){
		if(!f[i]) pri[++cnt]=i,mu[i]=-1;
		for(int j=1;j<=cnt&&i*pri[j]<=n;++j){
			f[i*pri[j]]=1;
			if(i%pri[j]==0){
				mu[i*pri[j]]=0;
				break;
			}
			mu[i*pri[j]]=-mu[i];
		}
	}
	rep(i,1,n) mu[i]+=mu[i-1];
}
inline ll sum_mu(ll x){
	if(!x) return 0;
	if(x<=mx) return mu[x];
	if(_mu[x]) return _mu[x];
	ll ret=0;
	for(ll l=2,r;l<=x;l=r+1){
		r=x/(x/l);
		ret+=(r-l+1)*sum_mu(x/l);
	}
	return _mu[x]=1-ret;
}
inline ll sum_phi(ll x,ll ret=0){
	for(ll l=1,r;l<=x;l=r+1){
		r=x/(x/l);
		ret+=(sum_mu(r)-sum_mu(l-1))*(x/l)*(x/l);
	}
	return ret+1>>1;
}
int main(){
	std::ios::sync_with_stdio(0);
	std::cin>>T;
	mx=5e6;
	init(mx);
	while(T--){
		ll n; std::cin>>n;
		std::cout<<sum_phi(n)<<' '<<sum_mu(n)<<'\n';
	}
	return 0;
}
```

杜教筛：

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)

#define ll long long
#define N 5000005
int T,mx,mu[N],pri[N],cnt;
ll phi[N];
bool f[N];
std::unordered_map<ll,ll> _mu,_phi;
inline void init(int n){
	mu[1]=phi[1]=1;
	rep(i,2,n){
		if(!f[i]) pri[++cnt]=i,mu[i]=-1,phi[i]=i-1;
		for(int j=1;j<=cnt&&i*pri[j]<=n;++j){
			f[i*pri[j]]=1;
			if(i%pri[j]==0){
				mu[i*pri[j]]=0;
				phi[i*pri[j]]=phi[i]*pri[j];
				break;
			}
			mu[i*pri[j]]=-mu[i];
			phi[i*pri[j]]=phi[i]*(pri[j]-1);
		}
	}
	rep(i,1,n) mu[i]+=mu[i-1],phi[i]+=phi[i-1];
}
inline ll sum_mu(ll x){
	if(x<=mx) return mu[x];
	if(_mu[x]) return _mu[x];
	ll ret=0;
	for(ll l=2,r;l<=x;l=r+1){
		r=x/(x/l);
		ret+=(r-l+1)*sum_mu(x/l);
	}
	return _mu[x]=1-ret;
}
inline ll sum_phi(ll x){
	if(x<=mx) return phi[x];
	if(_phi[x]) return _phi[x];
	ll ret=0;
	for(ll l=2,r;l<=x;l=r+1){
		r=x/(x/l);
		ret+=(r-l+1)*sum_phi(x/l);
	}
	return _phi[x]=(x*(x+1)>>1)-ret;
}
int main(){
	std::ios::sync_with_stdio(0);
	std::cin>>T;
	mx=5e6;
	init(mx);
	while(T--){
		ll n; std::cin>>n;
		std::cout<<sum_phi(n)<<' '<<sum_mu(n)<<'\n';
	}
	return 0;
}
```

