[$\color{purple}{\text{back to top}}$](https://cyn2006.github.io)

# 莫反学习笔记 $2$

我太菜了。

我是菜鸡。

持续更新。

以下为紫（或黑？）题集合。

**申明：对于整除分块套整除分块的复杂度其实应该是 $O(n^{\frac{3}{4}})$**。

## $\mathcal {Lemma}$

$1.1$, 狄利克雷卷积：$(f*g)(i)=\sum\limits_{d\mid i} f(d)g(\frac{i}{d})$.

## P3768 简单的数学题

非常的。。。

我调疯了！！！

参考[Crash的数字表格](https://www.luogu.com.cn/blog/CYN/mobius1)，我们可以得到：

$$
ANS=\sum_{i=1}^n\sum_{j=1}^nij\gcd(i,j)
$$

$$
ANS=\sum_{d=1}^n d^3 \sum_{x=1}^{\lfloor \frac{n}{d} \rfloor} \mu(x)g(\lfloor \frac{n}{dx}\rfloor)^2x^2
$$

$$
=\sum_{t=1}^n g(\lfloor\frac{n}{t}\rfloor)^2 \sum_{d\mid t} d^3\mu(\frac{t}{d})\times (\frac{t}{d})^2
$$

$$
=\sum_{t=1}^n g(\lfloor\frac{n}{t}\rfloor)^2 t^2 \varphi(t)
$$

令 $S(n)=\sum\limits_{i=1}^n i^2\varphi(i)$，

对于 $f(t)=t^2\varphi(t)$，我们将它转化成 $\operatorname{id}_2*f=\operatorname{id}_3$，

则令 $g=\operatorname{id}_2$，则 
$$
g(1)S(n)=\sum_{i=1}^n (f*g)(i) - \sum_{i=2}^n g(i)S(\lfloor \frac{n}{i} \rfloor)
$$
而 $g(1)=1$，所以
$$
S(n)=\sum_{i=1}^n (f*g)(i)-\sum_{i=2}^n g(i)S(\lfloor \frac{n}{i} \rfloor)
$$

$$
=\sum_{i=1}^n i^3-\sum_{i=2}^n i^2\sum_{j=1}^{\lfloor\frac{n}{i}\rfloor}j^2\varphi(j)
$$
则整除分块套杜教筛即可。

应注意：$f$ 是积性函数，并且杜教筛时直接筛出 $i^2 \varphi(i)$ 即可。

$Code:$

```cpp
#include<bits/stdc++.h>
#define ll long long
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define N 10000005
int pri[N],cnt,mx;
bool f[N];
ll phi[N],mod,inv2,inv6;
std::unordered_map<ll,ll> _phi;
inline ll slow(ll x){return (x%mod+mod)%mod;}
inline ll sum(ll x){
	x=slow(x);
	return x*(x+1)%mod*inv2%mod;
}
inline ll g(ll x){return sum(x)*sum(x)%mod;}
inline ll zzc(ll x){
	x=slow(x);
	return x*(x+1)%mod*(2*x+1)%mod*inv6%mod;
}

inline ll qpow(ll x,ll y,ll ret=1){
	for(;y;y>>=1,x=x*x%mod) if(y&1)
		ret=ret*x%mod; return ret;
}
inline void init(int n){
	phi[1]=1;
	rep(i,2,n){
		if(!f[i]) pri[++cnt]=i,phi[i]=i-1;
		for(int j=1;j<=cnt&&i*pri[j]<=n;++j){
			f[i*pri[j]]=1;
			if(i%pri[j]==0){
				phi[i*pri[j]]=slow(phi[i]*pri[j]);
				break;
			}
			phi[i*pri[j]]=slow(phi[i]*(pri[j]-1));
		}
	}
	rep(i,1,n) phi[i]=slow(phi[i]%mod*i%mod*i%mod+phi[i-1]);
}

inline ll sum_phi(ll x){
	if(x<=mx) return phi[x];
	if(_phi[x]) return _phi[x];
	ll ret=g(x);
	for(ll l=2,r;l<=x;l=r+1){
		r=x/(x/l);
		ret=slow(ret-slow(zzc(r)-zzc(l-1))*sum_phi(x/l)%mod);
	}
	return _phi[x]=ret;
}

inline ll solve(ll x){
	ll ret=0;
	for(ll l=1,r;l<=x;l=r+1){
		r=x/(x/l);
		ret=slow(ret+sum(x/l)*sum(x/l)%mod*slow(sum_phi(r)-sum_phi(l-1))%mod);
	}
	return ret;
}

int main(){
	mx=1e7; ll n;
	std::cin>>mod>>n;
	inv2=qpow(2,mod-2),inv6=qpow(6,mod-2);
	init(mx);
	std::cout<<solve(n)<<'\n';
	return 0;
}
```

## P5221 Product

题目要求

$$
\prod_{i=1}^n\prod_{j=1}^n \dfrac{\operatorname{lcm}(i,j)}{\gcd(i,j)}
$$

$$
=\prod_{i=1}^n\prod_{j=1}^n\dfrac{i\times j}{\gcd(i,j)^2}
$$

首先考虑求出分母是什么。

我们考虑先枚举 $d=\gcd(i,j)$，即
$$
\prod_{d=1}^nd\ ^{\large\sum\limits_{i=1}^n\sum\limits_{j=1}^n[\gcd(i,j)=d]}
$$

$$
\prod_{d=1}^nd\ ^{\large \sum\limits_{i=1}^{\lfloor \frac{n}{d} \rfloor}\sum\limits_{j=1}^{\lfloor \frac{n}{d} \rfloor}[\gcd(i,j)=1]}
$$

$$
\prod_{d=1}^n d\ ^{\large\sum\limits_{x=1}^{\lfloor \frac{n}{d} \rfloor} \mu(x) \sum\limits_{i=1}^{\lfloor \frac{n}{d} \rfloor}[x\mid i]\sum\limits_{j=1}^{\lfloor \frac{n}{d} \rfloor}[x\mid j]}
$$

$$
\prod_{d=1}^n d\ ^{\large\sum\limits_{x=1}^{\lfloor \frac{n}{d} \rfloor} \mu(x) \lfloor \frac{n}{dx} \rfloor^2}
$$
于是整除分块套整除分块即可。

然后出题人恶意卡空间就稍微卡下空间吧，反正是线性 $O(n)$ 的，实在不行杜教筛直接上即可。

$Code:$

```cpp
#include<bits/stdc++.h>
#define ll long long
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define N 1000005
const int mod=104857601;
int pri[78505],cnt,mu[N];
bool f[N];
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
inline ll qpow(ll x,ll y,ll ret=1){
	for(;y;y>>=1,x=x*x%mod) if(y&1)
		ret=ret*x%mod; return ret;
}
inline void reduce(ll&x){
	x+=x>>31&mod;
	x-=mod,x+=x>>31&mod;
}
inline ll solve2(ll x){
	ll ret=0;
	for(ll l=1,r;l<=x;l=r+1){
		r=x/(x/l);
		ret+=1ll*(mu[r]-mu[l-1])*(x/l)*(x/l);
	}
	return ret;
}
inline ll solve1(ll x){
	ll ret=1;
	for(ll l=1,r;l<=x;l=r+1){
		r=x/(x/l);
		ll tmp=solve2(x/l);
//		std::cerr<<x/l<<' '<<tmp<<'\n';
		rep(i,l,r) ret=1ll*ret*qpow(i,tmp)%mod;
	}
	return qpow(ret*ret%mod,mod-2);
}
int main(){
	int n; std::cin>>n;
	init(n);
	ll ans=1;
	rep(i,1,n) ans=1ll*ans*i%mod;
	ll tmp=ans; ans=1;
	rep(i,1,n) ans=ans*tmp%mod*qpow(i,n)%mod;
	ans=1ll*ans*solve1(n)%mod;
	std::cout<<ans<<'\n';
	return 0;
}
```

## P6156 简单题

题目所求

$$
ans(n,k)=\sum_{i=1}^n\sum_{j=1}^n (i+j)^k \mu^2(\gcd(i,j)) \gcd(i,j)
$$

$$
=\sum_{i=1}^n\sum_{j=1}^n\sum_{d=1}^n [\gcd(i,j)=d](i+j)^k\mu^2(\gcd(i,j))\gcd(i,j)
$$

$$
=\sum_{d=1}^n d^{k}\sum_{i=1}^{\lfloor \frac{n}{d} \rfloor}\sum_{j=1}^{\lfloor \frac{n}{d} \rfloor} [\gcd(i,j)=1](i+j)^k\mu^2(d)d
$$

$$
=\sum_{d=1}^n \mu^2(d)d^{k+1}\sum_{i=1}^{\lfloor \frac{n}{d} \rfloor}\sum_{j=1}^{\lfloor \frac{n}{d} \rfloor}[\gcd(i,j)=1](i+j)^k
$$

$$
=\sum_{d=1}^n \mu^2(d)d^{k+1}\sum_{i=1}^{\lfloor \frac{n}{d} \rfloor}\sum_{j=1}^{\lfloor \frac{n}{d} \rfloor}\sum_{x\mid i,x\mid j}\mu(x)(i+j)^k
$$

$$
=\sum_{d=1}^n \mu^2(d)d^{k+1}\sum_{x=1}^{\lfloor \frac{n}{d} \rfloor}\mu(x)x^k\sum_{i=1}^{\lfloor \frac{n}{dx} \rfloor}\sum_{j=1}^{\lfloor \frac{n}{dx} \rfloor} (i+j)^k
$$

令 $S(n)=\sum\limits_{i=1}^n\sum\limits_{j=1}^n (i+j)^k$，则有
$$
ans(n,k)=\sum_{dx=t=1}^n t^kS(\lfloor \frac{n}{t} \rfloor) \sum_{d\mid t} \mu^2(d)\mu(\frac{t}{d}) d
$$
然后我们开始发掘性质。

首先我们知道 $\mu^2(x)\in\{0,1\}$，因此对于所有的有质因子平方为因数的数，他们的 $\mu$ 都是 $0$，否则都是 $1$。

然后我们考虑线性筛的时候顺带筛出幂函数和后式，记现在筛到了 $i$，枚举质数为 $p$，则：

1. $p\not\mid i$：正常线性筛转移即可，即 `f[i*pri[j]]=f[i]*(pri[j]-1)`；
2. $p\mid i$ 但 $p^2\not\mid i$：考虑到 $q=\dfrac{i}{p}$ 它不含平方质因子，则转移时应为 `f[i*pri[j]]=-f[i/p]*pri[j]`（即 $\mu=-1$）；
3. $p^2 \mid i$：考虑到它必然有平方质因子，因此本次无贡献，不需要转移。

这样我们已经求完了后面的和式了，简记为 $f(x)$。

现在我们考虑快速求出 $g(x)=f(x)x^k$，这样就能在 $O(\sqrt n)$ 内进行一次整除分块了。

我们考虑直接对于每一个 $f(x)$ 乘上线性筛时的 $x^k$ 即可。

然后现在考虑如何快速计算（比如在 $O(1)$ 的时间内求出 $S(n)$）$S$。
$$
S(n)=\sum_{i=1}^n i^k+\sum_{i=n+1}^{2n} i^k
$$
于是处理一下幂次方的二次前缀和即可 $O(1)$ 计算。

于是预处理用线性筛时间复杂度大约为 $O(n)$，单次查询时间复杂度为 $O(\sqrt{n})$，加强版卡卡常应该（/yiw）就能过去了。

$Code:$ （加强版没卡过去的代码，待卡常）

```cpp
#include<bits/stdc++.h>
#define ll long long
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define N 20000005
typedef unsigned u32;
int pri[N],cnt;
u32 f[N],pw[N];
bool fg[N];
int n,k,q;
inline u32 qpow(u32 x,int y,u32 ret=1){
	for(;y;y>>=1,x=x*x) if(y&1)
		ret=ret*x; return ret;
}
inline void init(int n){
	f[1]=pw[1]=1;
	rep(i,2,n){
		if(!fg[i]) pri[++cnt]=i,pw[i]=qpow(i,k),f[i]=i-1;
		for(int j=1;j<=cnt&&i*pri[j]<=n;++j){
			const int v=i*pri[j];
			fg[v]=1,pw[v]=pw[i]*pw[pri[j]];
			if(i%pri[j]==0){
				int x=i/pri[j];
				if(x%pri[j]) f[v]=-pri[j]*f[x];
				break;
			}
			f[v]=f[i]*(pri[j]-1);
		}
	}
	rep(i,2,n) f[i]=f[i-1]+f[i]*pw[i],pw[i]+=pw[i-1];
	rep(i,2,n) pw[i]+=pw[i-1];
}
inline u32 calc(int n){
	return pw[n<<1]-pw[n]*2;
}
template<typename T>inline void rd(T&x){
	x=0; char c=getchar();
	while(!isdigit(c)) c=getchar();
	while(isdigit(c)) x=x*10+(c&15),c=getchar();
}
inline void wr(u32 x){
	if(x>9) wr(x/10); putchar(x%10|48);
}
int main(){
	std::ios::sync_with_stdio(0);
	int T; rd(T),rd(n),rd(k);
	init(n<<1);
	while(T--){
		int newn; rd(newn);
		u32 ret=0;
		for(ll l=1,r;l<=newn;l=r+1){
			r=newn/(newn/l);
			ret+=calc(newn/l)*(f[r]-f[l-1]);
		}
		wr(ret),putchar('\n');
	}
	return 0;
}
```

然后现在卡过去了（循环变量用 u32）：

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define N 20000005
typedef unsigned u32;
int pri[N],cnt;
u32 f[N],pw[N];
bool fg[N];
int n,k,q;
inline u32 qpow(u32 x,int y,u32 ret=1){
	for(;y;y>>=1,x=x*x) if(y&1)
		ret=ret*x; return ret;
}
inline void init(int n){
	f[1]=pw[1]=1;
	rep(i,2,n>>1){
		if(!fg[i]) pri[++cnt]=i,pw[i]=qpow(i,k),f[i]=i-1;
		for(int j=1;j<=cnt&&i*pri[j]<=n;++j){
			fg[i*pri[j]]=1,pw[i*pri[j]]=pw[i]*pw[pri[j]];
			if(i%pri[j]==0){
				int x=i/pri[j];
				if(x%pri[j]) f[i*pri[j]]=-pri[j]*f[x];
				break;
			}
			f[i*pri[j]]=f[i]*(pri[j]-1);
		}
	}
	rep(i,(n>>1)+1,n){
		if(!fg[i]) pri[++cnt]=i,pw[i]=qpow(i,k);
		for(int j=1;j<=cnt&&i*pri[j]<=n;++j){
			fg[i*pri[j]]=1,pw[i*pri[j]]=pw[i]*pw[pri[j]];
			if(i%pri[j]==0) break;
		}
	}
	rep(i,2,n>>1) f[i]=f[i-1]+f[i]*pw[i];
	rep(i,2,n) pw[i]+=pw[i-1];
	rep(i,2,n) pw[i]+=pw[i-1];
}
inline u32 calc(int n){
	return pw[n<<1]-pw[n]*2;
}
template<typename T>inline void rd(T&x){
	x=0; char c=getchar();
	while(!isdigit(c)) c=getchar();
	while(isdigit(c)) x=x*10+(c&15),c=getchar();
}
inline void wr(u32 x){
	if(x>9) wr(x/10); putchar(x%10|48);
}
int main(){
	std::ios::sync_with_stdio(0);
	int T; rd(T),rd(n),rd(k);
	init(n<<1);
	while(T--){
		int newn; rd(newn);
		u32 ret=0;
		for(u32 l=1,r;l<=newn;l=r+1){
			r=newn/(newn/l);
			ret+=calc(newn/l)*(f[r]-f[l-1]);
		}
		wr(ret),putchar('\n');
	}
	return 0;
}
```

## P3704 [SDOI2017]数字表格

参考前文的 [Product](https://www.luogu.com.cn/problem/P5221)，我们同样可以用枚举 $d$ 的方式来解决此题。

考虑到原式可以转化为
$$
\sum_{d=1}^n f(d)\ ^{ \sum\limits_{i=1}^n\sum\limits_{j=1}^m[\gcd(i,j)=d]}
$$
把上面的次方提出来：
$$
\sum_{i=1}^n\sum_{j=1}^m[\gcd(i,j)=d]
$$
然后就是一个裸的 $\mathcal{Mobius}$ 反演：
$$
\sum_{x=1}^n\mu(x)\lfloor \frac{n}{dx}\rfloor \lfloor \frac{m}{dx}\rfloor
$$
考虑枚举 $t=dx$，题目需求的原式转化为：
$$
\sum_{t=1}^n\sum_{d\mid t} f(d)\ ^{\mu(\frac{t}{d})\lfloor \frac{n}{t}\rfloor \lfloor \frac{m}{t}\rfloor}
$$
令 $g(t)=\sum_{d\mid t}f(d)\ ^{\mu(\frac{t}{d})}$，则有原式：
$$
\sum_{t=1}^n g(t)^{\lfloor \frac{n}{t}\rfloor \lfloor \frac{m}{t}\rfloor}
$$
那么 $O(n\log n)$ 预处理 $g(t)$，每次询问整除分块即可。

时间复杂度：$O(n\log n+q\sqrt{n}\log \omega)$，其中 $\omega$ 为值域。

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define ll long long
#define N 1000005
const int mod=1e9+7;
int T,n,fib[N],inv[N];
int pri[N],cnt,mu[N],mx,g[N];
bool f[N];
inline int fpow(int x,int y){
	if(y==0) return 1;
	if(y==1) return fib[x];
	return inv[x];
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
	rep(i,0,n) g[i]=1;
	rep(i,1,n) for(int j=i;j<=n;j+=i)
		g[j]=1ll*g[j]*fpow(i,mu[j/i])%mod;
	rep(i,1,n) g[i]=1ll*g[i-1]*g[i]%mod;
}
inline int qpow(int x,int y,int ret=1){
	for(;y;y>>=1,x=1ll*x*x%mod) if(y&1)
		ret=1ll*ret*x%mod; return ret;
}
inline int solve(int x,int y,int ret=1){
	for(int l=1,r;l<=std::min(x,y);l=r+1){
		r=std::min(x/(x/l),y/(y/l));
		ret=1ll*ret*qpow(1ll*qpow(g[l-1],mod-2)*g[r]%mod,1ll*(x/l)*(y/l)%(mod-1))%mod;
	}
	return ret;
}
int main(){
	std::cin>>T;
	mx=1e6; fib[1]=1;
	rep(i,2,mx) fib[i]=(fib[i-1]+fib[i-2])%mod;
	rep(i,0,mx) inv[i]=qpow(fib[i],mod-2);
	init(mx);
	while(T--){
		int n,m; std::cin>>n>>m;
		std::cout<<solve(n,m)<<'\n';
	}
	return 0;
}
```

## P5518 [MtOI2019]幽灵乐团

> 幽灵乐团名副其实。/kk

精神污染。。。（迫真）

有事不闲的人千万别来做这题。。。

以下我们大致钦定 $A\leqslant B\leqslant C$。

但是对于 $type=2$ 似乎要对最大最小值特殊处理一下。

### $type=0$

$$
\prod_{i=1}^A \prod_{j=1}^B\prod_{k=1}^C \frac{ij}{\gcd(i,j)\gcd(i,k)}
$$

把分母取出来：
$$
\prod_{i=1}^A \prod_{j=1}^B\prod_{k=1}^C \gcd(i,j)\times \prod_{i=1}^A\prod_{j=1}^B\prod_{k=1}^C\gcd(i,k)
$$
考虑左右边的本质是一样的，那么单独计算左边：
$$
\prod_{i=1}^A\prod_{j=1}^B\prod_{k=1}^C \gcd(i,j)
$$

$$
=\prod_{k=1}^C \prod_{i=1}^A \prod_{j=1}^B \gcd(i,j)
$$

$$
=\left( \prod _{d=1}^{A} d\ ^{\sum\limits_{i=1}^A\sum\limits_{j=1}^B [\gcd(i,j)=d]}\right)^C
$$

$$
=\left(\prod_{d=1}^Ad\ ^{\sum\limits_{x=1}^{\lfloor \frac{A}{d} \rfloor} \mu(x)\lfloor \frac{A}{dx} \rfloor \lfloor \frac{B}{dx} \rfloor}\right)
$$
右边同理。

那么整除分块即可。

### $type=1$

同 $type=0$，先处理左式。我们有：
$$
 \prod_{k=1}^C\left(\prod_{d=1}^A d\ ^{\sum\limits_{i=1}^A\sum\limits_{j=1}^B [\gcd(i,j)=d]ij}\right)^k
$$
把后式取出来：
$$
 \left(\prod_{d=1}^A d\ ^{d^2\sum\limits_{i=1}^{\lfloor \frac{A}{d} \rfloor} \sum\limits_{j=1}^{\lfloor \frac{B}{d}\rfloor} [\gcd(i,j)=1]ij} \right)^k
$$

$$
 =\left(\prod_{d=1}^A d\ ^{d^2\sum\limits_{i=1}^{\lfloor \frac{A}{d} \rfloor} \sum\limits_{j=1}^{\lfloor \frac{B}{d}\rfloor} \sum\limits_{x\mid i,x\mid j} \mu(x)ij} \right)^k
$$


把上面的次方取出来，可以发现是一个 [Crash 的数字表格] 的形式：

令 $S(x)=\sum\limits_{i=1}^x i$，则有
$$
 d^2\sum_{x=1}^{\lfloor \frac{A}{d}\rfloor} \mu(x) x^2S(\lfloor \frac{A}{dx}\rfloor) S(\lfloor \frac{B}{dx} \rfloor)
$$
那么整除分块即可。

### $type=2$

考虑其分母的形式：
$$
 \prod_{d=1}^A d\ ^{\sum\limits_{i=1}^A\sum\limits_{j=1}^B [\gcd(i,j)=d]\sum\limits_{k=1}^C\gcd(d,k)}
$$

$$
 =\prod_{d=1}^A d\ ^{\sum\limits_{x=1}^{\lfloor \frac{A}{d} \rfloor} \mu(x) \lfloor \frac{A}{dx}\rfloor \lfloor \frac{B}{dx} \rfloor \sum\limits_{k=1}^C \gcd(d,k)}
$$
考虑幂次的右式：
$$
 \sum_{k=1}^C \gcd(d,k)
$$

$$
 =\sum_{k=1}^C\sum_{x=1}^A [\gcd(d,k)=x]x
$$

$$
 =\sum_{x\mid d}x \sum_{k=1}^{\lfloor \frac{A}{x} \rfloor} [\gcd(\frac{d}{x},k)=1]
$$

$$
 =\sum_{x\mid d}x\sum_{k=1}^{\lfloor \frac{A}{x} \rfloor}\sum_{y\mid \frac{d}{x},y\mid k} \mu(y)
$$

$$
 =\sum_{x\mid d}x\sum_{y\mid \frac{d}{x}} \mu(y)\sum_{k=1}^{\lfloor \frac{A}{x} \rfloor} [y\mid k]
$$

$$
 =\sum_{x\mid d}x\sum_{y\mid \frac{d}{x}} \mu(y) \lfloor \frac{A}{xy} \rfloor
$$

$$
 =\sum_{t=xy\mid d} \sum_{y\mid t} \mu(y) \frac{t}{y}\lfloor \frac{A}{t} \rfloor
$$

$$
 =\sum_{t\mid d} \varphi(t)\lfloor \frac{A}{t} \rfloor
$$
所以现在分母化成了这幅鬼样：
$$
 \prod_{i=1}d\ ^{\sum\limits_{x=1}^{\lfloor \frac{A}{d} \rfloor} \mu(x) \lfloor \frac{A}{dx}\rfloor \lfloor \frac{B}{dx} \rfloor \sum\limits_{t\mid d} \varphi(t)\lfloor \frac{A}{t} \rfloor}
$$
我已经不想再继续推下去了。。。

然后考虑分子的形式，我们把它提出来：
$$
 \prod_{i=1}^A \prod_{j=1}^B \prod_{k=1}^C ij\ ^{\gcd(i,j,k)}
$$
考虑到 $i$ 与 $j$ 其实可以互相独立的，因此
$$
 \prod_{i=1}^A \prod_{j=1}^B \prod_{k=1}^C i\ ^{\gcd(i,j,k)}
$$

$$
 =\prod_{i=1}^A i\ ^{\sum\limits_{j=1}^B\sum\limits_{k=1}^C\gcd(i,j,k)}
$$
然后考虑单独拆除 $i$ 的幂次：
$$
\sum\limits_{j=1}^B\sum\limits_{k=1}^C \gcd(i,j,k)
$$

$$
 =\sum_{d=1}^A d \sum_{j=1}^{\lfloor \frac{B}{d} \rfloor} \sum_{k=1}^{\lfloor \frac{C}{d} \rfloor} [\gcd(i,j,k)=1]
$$

$$
 =\sum_{d\mid i} d\sum_{x=1}^{\lfloor \frac{B}{d} \rfloor} \mu(x)  \lfloor \frac{B}{dx} \rfloor \lfloor \frac{C}{dx} \rfloor
$$

$$
 =\sum_{dx=t=1}^B \lfloor \frac{B}{t} \rfloor \lfloor \frac{C}{t} \rfloor\sum_{d\mid t,d\mid i} \mu(\frac{t}{d}) d
$$

$$
 =\sum_{dx=t=1}^B \lfloor \frac{B}{t} \rfloor \lfloor \frac{C}{t} \rfloor \varphi(\gcd(i,t))
$$
然后我们发现我们彻底地推挂了。。。

然后我们考虑换一种思维

------

我们现在先考虑分子，我们重新地推一下：
$$
 \prod_{i=1}^A \prod_{j=1}^B \prod_{k=1}^C i\ ^{\gcd(i,j,k)}
$$
我们先枚举 $\gcd(i,j,k)$，大概长这样（其中我们记 $fac(x)=\prod\limits_{i=1}^x i$）：
$$
 \prod_{d=1}^{\min\{A,B,C\}} \prod_{i=1}^A i\ ^{d\sum\limits_{j=1}^B\sum\limits_{k=1}^C [\gcd(i,j,k)=d]}
$$

$$
 =\prod_{d=1}^{\min\{A,B,C\}} \prod_{i=1}^{\lfloor \frac{A}{d} \rfloor} (di)\ ^{d\sum\limits_{j=1}^{\lfloor \frac{B}{d} \rfloor}\sum\limits_{k=1}^{\lfloor \frac{C}{d} \rfloor}[\gcd(i,j,k)=1]}
$$

$$
 =\prod_{d=1}^{\min\{A,B,C\}} \prod_{x=1}^{\min\{A,B,C\}} \prod_{i=1}^{\lfloor \frac{A}{dx} \rfloor} (dxi)\ ^{\mu(x)d \lfloor \frac{B}{dx} \rfloor \lfloor \frac{C}{dx} \rfloor}
$$

$$
 =\prod_{t=1}^{\min\{A,B,C\}} \prod_{d\mid t} (t^{\lfloor \frac{A}{t}\rfloor}\times fac(\lfloor \frac{A}{t} \rfloor)\ ^{\mu(\frac{t}{d}) d\lfloor \frac{B}{t} \rfloor \lfloor \frac{B}{t} \rfloor}
$$

$$
 =\prod_{t=1}^{\min\{A,B,C\}} (t^{\lfloor \frac{A}{t}\rfloor}\times fac(\lfloor \frac{A}{t} \rfloor)\ ^{\lfloor \frac{B}{t} \rfloor \lfloor \frac{C}{t}\rfloor\sum\limits_{d\mid t}\mu(\frac{t}{d}) d }
$$

$$
 =\prod_{t=1}^{\min\{A,B,C\}} (t^{\lfloor \frac{A}{t}\rfloor}\times fac(\lfloor \frac{A}{t} \rfloor)\ ^{\lfloor \frac{B}{t} \rfloor \lfloor \frac{C}{t}\rfloor\varphi(t)}
$$

$$
 =\prod_{t=1}^{\min\{A,B,C\}} \left(\left(t^{\lfloor \frac{A}{t}\rfloor}\times fac(\lfloor \frac{A}{t} \rfloor)\right)\ ^{\varphi(t)} \right)^{\lfloor \frac{B}{t} \rfloor \lfloor \frac{C}{t}\rfloor}
$$
然后整除分块即可。

我们现在来重新考虑分母：
$$
 \prod_{i=1}^A \prod_{j=1}^B \prod_{k=1}^C \gcd(i,j)\ ^{\gcd(i,j,k)}
$$
考虑枚举 $d=\gcd(i,j,k)$ 得到：
$$
 \prod_{d=1}^{\min\{A,B,C\}} \prod_{i=1}^A\prod_{j=1}^B \gcd(i,j)\ ^{\sum\limits_{k=1}^C [\gcd(i,j,k)=d]}
$$

$$
 =\prod_{d=1}^{\min\{A,B,C\}} \prod_{i=1}^{\lfloor \frac{A}{d} \rfloor} \prod_{j=1}^{\lfloor \frac{B}{d} \rfloor} (d\gcd(i,j))^{d\sum\limits_{k=1}^C [\gcd(i,j,k)=1]}
$$

$$
 =\prod_{d=1}^{\min\{A,B,C\}} \prod_{d_1=1}^{\min\{\lfloor \frac{A}{d} \rfloor,\lfloor \frac{B}{d} \rfloor\}}\prod_{i=1}^{\lfloor \frac{A}{dd_1}\rfloor} \prod_{j=1}^{\lfloor \frac{B}{dd_1}\rfloor} \left(dd_1\gcd(i,j)\right)^{d\mu(d_1)\lfloor \frac{C}{dd_1} \rfloor}
$$
将 $dd_1$ 与 $\gcd(i,j)$ 的乘积分离开来，我们先计算左边的乘积得到：
$$
 \prod_{d=1}^{\min\{A,B,C\}} \prod_{d_1=1}^{\min\{\lfloor \frac{A}{d} \rfloor,\lfloor \frac{B}{d} \rfloor\}} \left(\left(dd_1\right)\ ^{d\mu(d_1)\lfloor \frac{C}{dd_1} \rfloor}\right)^{\lfloor \frac{A}{dd_1}\rfloor \lfloor \frac{B}{dd_1}\rfloor}
$$

$$
 =\prod_{dd_1=t=1}^{\min\{A,B,C\}} \prod_{x\mid t} t\ ^{x\mu(\frac{t}{x})\lfloor \frac{C}{dd_1} \rfloor \lfloor \frac{A}{dd_1}\rfloor \lfloor \frac{B}{dd_1}\rfloor}
$$

$$
 =\prod_{dd_1=t=1}^{\min\{A,B,C\}} t\ ^{\varphi(t)\lfloor \frac{C}{t} \rfloor \lfloor \frac{A}{t}\rfloor \lfloor \frac{B}{t}\rfloor}
$$
然后这个应该是已经可以整除分块了。

然后我们考虑另一个右边的乘积：
$$
 \prod_{d=1}^{\min\{A,B,C\}} \prod_{d_1=1}^{\min\{\lfloor \frac{A}{d} \rfloor,\lfloor \frac{B}{d} \rfloor\}}\prod_{i=1}^{\lfloor \frac{A}{dd_1}\rfloor} \prod_{j=1}^{\lfloor \frac{B}{dd_1}\rfloor} \gcd(i,j)\ ^{d\mu(d_1)\lfloor \frac{C}{dd_1} \rfloor}
$$

$$
 =\prod_{t=1}^C \left(\prod_{t_1=1}^{\min\{\frac{A}{t},\frac{B}{t}\}} \left( \prod_{d\mid t_1} d^{\mu(\frac{t_1}{d})}\right)^{\lfloor \frac{A}{tt_1}\rfloor \lfloor \frac{B}{tt_1}\rfloor}\right)^{\lfloor \frac{C}{t} \rfloor \varphi(t)}
$$
然后整除分块套整除分块即可。

时间复杂度大概依赖于 $type=2$，所以总时间大概为 $O(n\log n+Tn^{\frac{3}{4}}\log n)$。

总代码：

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define ll long long
#define int long long
#define N 100005
int T,mod,A,B,C,mu[N],mu2[N],pri[N],cnt,pw[N],pw1[N],pw2[N];
int phi[N],pwphi[N],dmu[N],inv[N],dmuinv[N];
bool f[N];
inline int qpow(int x,int y,int ret=1){
	for(;y;y>>=1,x=1ll*x*x%mod) if(y&1)
		ret=1ll*ret*x%mod; return ret;
}
inline void reduce(int&x,const int mod){
	x+=x>>31&mod,x-=mod,x+=x>>31&mod;
}
inline int upd(int x,int p){
	return x+(x>>31&p);
}
inline void init(int n){
	mu[1]=phi[1]=1;
	inv[0]=1;
	rep(i,1,n) inv[i]=qpow(i,mod-2);
	rep(i,2,n){
		if(!f[i]) pri[++cnt]=i,mu[i]=-1,phi[i]=i-1;
		for(int j=1;j<=cnt&&i*pri[j]<=n;++j){
			f[i*pri[j]]=1;
			if(i%pri[j]==0){mu[i*pri[j]]=0; phi[i*pri[j]]=phi[i]*pri[j]; break;}
			mu[i*pri[j]]=-mu[i],phi[i*pri[j]]=phi[i]*(pri[j]-1);
		}
	}
	const int p=mod-1;
	pwphi[0]=dmu[0]=1;
	rep(i,1,n) pwphi[i]=qpow(i,phi[i]),dmu[i]=1;
	rep(i,1,n) reduce(phi[i]+=phi[i-1],p),pwphi[i]=1ll*phi[i]*pwphi[i-1]%mod;
	rep(i,1,n) if(mu[i]){
		for(int j=i<<1;j<=n;j+=i)
			if(mu[i]==1) dmu[j]=1ll*dmu[j]*(j/i)%mod;
			else dmu[j]=1ll*dmu[j]*inv[j/i]%mod;
	}
	rep(i,1,n) dmu[i]=1ll*dmu[i]*dmu[i-1]%mod;
	dmuinv[0]=1;
	rep(i,1,n) dmuinv[i]=qpow(dmu[i],mod-2);
	rep(i,1,n) mu2[i]=1ll*mu[i]*i%p*i%p;
	rep(i,1,n) reduce(mu[i]+=mu[i-1],p),reduce(mu2[i]+=mu2[i-1],p);
	rep(i,1,n) reduce(mu2[i]+=p,p);
	pw[0]=pw1[0]=pw2[0]=1;
	rep(i,1,n) pw[i]=i,pw1[i]=qpow(i,i),pw2[i]=qpow(i,1ll*i*i%p);
	rep(i,1,n){
		pw[i]=1ll*pw[i-1]*pw[i]%mod;
		pw1[i]=1ll*pw1[i-1]*pw1[i]%mod;
		pw2[i]=1ll*pw2[i-1]*pw2[i]%mod;
	}
}

namespace type0{
inline int solve2(int x,int y){
	ll ret=0;
	const int p=mod-1;
	for(int l=1,r;l<=std::min(x,y);l=r+1){
		r=std::min(x/(x/l),y/(y/l));
		(ret+=1ll*upd(mu[r]-mu[l-1],p)*(x/l)%p*(y/l))%=p;
	}
	return ret;
}
inline int solve1(int x,int y){
	int ret=1;
	for(int l=1,r;l<=std::min(x,y);l=r+1){
		r=std::min(x/(x/l),y/(y/l));
		int tmp=solve2(x/l,y/l);
		ret=1ll*ret*qpow(1ll*pw[r]*qpow(pw[l-1],mod-2)%mod,tmp)%mod;
	}
	return ret;
}
void brute(){
	int ret=1ll*qpow(pw[A],B)*qpow(pw[B],A)%mod;
	ret=qpow(ret,C);
	int inv=1ll*qpow(solve1(A,B),C)*qpow(solve1(A,C),B)%mod;
	inv=qpow(inv,mod-2);
	std::cout<<1ll*ret*inv%mod<<' ';
}
}
namespace type1{
inline int g(int x){
	return (x*(x+1ll)>>1)%(mod-1);
}
inline int solve2(int x,int y){
	int ret=0;
	const int p=mod-1;
	for(int l=1,r;l<=std::min(x,y);l=r+1){
		r=std::min(x/(x/l),y/(y/l));
		(ret+=1ll*(mu2[r]-mu2[l-1])*g(x/l)%p*g(y/l)%p)%=p;
	}
	return (ret%p+p)%p;
}
inline int solve1(int x,int y){
	int ret=1;
	for(int l=1,r;l<=std::min(x,y);l=r+1){
		r=std::min(x/(x/l),y/(y/l));
		int tmp=solve2(x/l,y/l);
		ret=1ll*ret*qpow(1ll*pw2[r]*qpow(pw2[l-1],mod-2)%mod,tmp)%mod;
	}
	return ret;
}
void brute(){
	int ret=1ll*qpow(pw1[A],g(B))*qpow(pw1[B],g(A))%mod;
	ret=qpow(ret,g(C));
	int inv=1ll*qpow(solve1(A,B),g(C))*qpow(solve1(A,C),g(B))%mod;
	inv=qpow(inv,mod-2);
	ret=1ll*ret*inv%mod;
	std::cout<<ret<<' ';
}
}
namespace type2{
inline int min(int x,int y,int z){
	return std::min(std::min(x,y),z);
}
inline int up(int x,int y,int z){
	int ret=1;
	const int p=mod-1;
	for(int l=1,r;l<=min(x,y,z);l=r+1){
		r=min(x/(x/l),y/(y/l),z/(z/l));
		ret=1ll*ret*qpow(pw[x/l],1ll*(y/l)*(z/l)%p*upd(phi[r]-phi[l-1],p)%p)%mod;
	}
	return ret;
}
inline int down2(int x,int y){
	int ret=1;
	const int p=mod-1;
	for(int l=1,r;l<=std::min(x,y);l=r+1){
		r=std::min(x/(x/l),y/(y/l));
		ret=1ll*ret*qpow(1ll*dmuinv[l-1]*dmu[r]%mod,1ll*(x/l)*(y/l)%p)%mod;
	}
	return ret;
}
inline int down1(int x,int y,int z){
	int ret=1;
	const int p=mod-1;
	for(int l=1,r;l<=min(x,y,z);l=r+1){
		r=min(x/(x/l),y/(y/l),z/(z/l));
		ret=1ll*ret*qpow(down2(x/l,y/l),1ll*(z/l)*upd(phi[r]-phi[l-1],p)%p)%mod;
	}
	return ret;
} 
void brute(){
	int ret=1ll*up(A,B,C)*up(B,A,C)%mod;
	int tmp=1ll*down1(A,B,C)*down1(A,C,B)%mod;
	tmp=qpow(tmp,mod-2);
	std::cout<<1ll*ret*tmp%mod<<'\n'; 
}
}
signed main(){
	std::cin>>T>>mod;
	init(1e5);
	rep(i,1,T){
		std::cin>>A>>B>>C;
		type0::brute();
		type1::brute();
		type2::brute();
	}
	return 0;
}
```

