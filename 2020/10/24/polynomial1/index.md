[$\color{purple}\text{back to top}$](https://cyn2006.github.io)

# 多项式学习笔记

## $\mathscr{START}$

### $\mathfrak{FFT}$

#### $\mathfrak{C}\mathcal{onvolution}$

卷积（$\mathfrak{convolution}$），其广义定义为：
$$
h(x)=\int_{-\infty}^{+\infty} g(i)\cdot f(x-i) di
$$
我们称 $h(x)$ 为 $f(x)$ 和 $g(x)$ 的卷积，在这里记作 $h=f*g$ 。

在整次数多项式域内我们有化简版本：
$$
h(x)=(f*g)(x)=\sum_{i=0}^n \sum_{j=0}^i a_jb_{i-j} x^{i}
$$
显然在这里的多项式均为 $n$ 次多项式，则卷之后为 $2n-1$ 次多项式。

#### 点值表示法

在平面上取 $n+1$ 个点能够确定一个 $n$ 次的函数。

那么我们取 $\{x_0,\ldots x_n\}$ 及 $\{y_0,\ldots y_n\}$ 和 $\{y'_0,\ldots y'_n\}$。

由于多项式乘法后的次数为 $2n-1$，因此我们把剩下的位数补齐可以得到：
$$
\{(x_0,y_0),\ldots (x_{2n-1},y_{2n-1})\}\\
\{(x_0,y'_0),\ldots (x_{2n-1},y'_{2n-1})\}
$$
乘起来得到
$$
\{(x_0,y_0y'_0),\ldots (x_{2n-1},y_{2n-1}y'_{2n-1})\}
$$
然后我们将它转化成系数表达式即可。

#### $\mathfrak{DFT\&IDFT}$

鸽了。

然后这只老鸽子原理不懂就开始瞎做题了。。。

#### 板子

NTT：

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define ll long long
#define N 5000005
const int mod=998244353,G=3,_inv=mod+1>>1;
inline void reduce(int&x){x-=mod,x+=x>>31&mod;}
//inline int upd(int x){x-=mod; return x+=x>>31&mod;}
inline int upd(int x){return x+=x>>31&mod;}
int wn[N],rev[N];
inline int qpow(int x,int y,int ret=1){
	for(;y;y>>=1,x=1ll*x*x%mod) if(y&1)
		ret=1ll*ret*x%mod; return ret;
}
namespace ntt{
	inline void install(int n){
		for(int mid=1;mid<n;mid<<=1){
			int w=qpow(G,mod/mid>>1);
			wn[mid]=1;
			rep(i,1,mid-1) wn[mid+i]=1ll*wn[mid+i-1]*w%mod;
		}
	}
	int lim,invlim;
	inline void init(int n){
		for(lim=invlim=1;lim<n;lim<<=1)
			invlim=1ll*invlim*_inv%mod;
		rep(i,1,lim-1) rev[i]=rev[i>>1]>>1|(lim>>1&-(i&1));
	}
	inline void ntt(int *a,int ty){
		static int b[N];
		rep(i,0,lim-1) b[i]=a[rev[i]];
		for(int mid=1;mid<lim;mid<<=1){
			for(int i=0;i<lim;i+=mid<<1){
				rep(j,0,mid-1){
					const int x=1ll*b[i+j+mid]*wn[mid+j]%mod;
					b[i+j+mid]=upd(b[i+j]-x);
					reduce(b[i+j]+=x);
				}
			}
		}
		if(!ty){
			rep(i,0,lim-1) a[i]=1ll*b[i]*invlim%mod;
			std::reverse(a+1,a+lim);
		} else {
			std::memcpy(a,b,lim<<2);
		}
	}
	inline void _ntt(int *a,int *b){
		ntt(a,1),ntt(b,1);
		rep(i,0,lim-1) a[i]=1ll*a[i]*b[i]%mod;
		ntt(a,0);
	}
}
int a[N],b[N],n,m;
int main(){
	std::ios::sync_with_stdio(0),std::cin.tie(0),std::cout.tie(0);
	std::cin>>n>>m;
	rep(i,0,n) std::cin>>a[i];
	rep(i,0,m) std::cin>>b[i];
	ntt::install(n+m+2<<1),ntt::init(n+m+1);
	ntt::_ntt(a,b);
	rep(i,0,n+m) std::cout<<a[i]<<" \n"[i==n+m];
	return 0;
}
```

FFT：

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define db double
#define ll long long
#define N 5000005
const db pi=acos(-1);
struct poly{
	db x,y;
	inline poly operator+(const poly&p)const{return poly{x+p.x,y+p.y};}
	inline poly operator-(const poly&p)const{return poly{x-p.x,y-p.y};}
	inline poly operator*(const poly&p)const{return poly{x*p.x-y*p.y,x*p.y+y*p.x};}
} a[N],b[N];

namespace fft{
	int lim,rev[N];
	inline void install(int n){
		int len=0; lim=1;
		while(lim<=n) lim<<=1,++len;
		rep(i,0,lim-1) rev[i]=(rev[i>>1]>>1)|((i&1)<<(len-1));
	}
	inline void fft(poly *a,int ty){
		rep(i,0,lim-1) if(i<rev[i]) std::swap(a[i],a[rev[i]]);
		for(int mid=1;mid<lim;mid<<=1){
			poly wn=poly{cos(pi/mid),ty*sin(pi/mid)};
			for(int i=mid<<1,j=0;j<lim;j+=i){
				poly w=poly{1,0};
				for(int k=0;k<mid;++k,w=w*wn){
					poly x=a[j+k],y=w*a[j+mid+k];
					a[j+k]=x+y,a[j+mid+k]=x-y;
				}
			}
		}
		if(!~ty){
			rep(i,0,lim) a[i].x/=lim;
		}
	}
	inline void _fft(poly *a,poly *b){
		fft(a,1),fft(b,1);
		rep(i,0,lim) a[i]=a[i]*b[i];
		fft(a,-1);
	}
}
int n,m;
int main(){
	std::ios::sync_with_stdio(0),std::cin.tie(0),std::cout.tie(0);
	std::cin>>n>>m;
	rep(i,0,n) std::cin>>a[i].x;
	rep(i,0,m) std::cin>>b[i].x;
	fft::install(n+m);
	fft::_fft(a,b);
	rep(i,0,n+m) std::cout<<(int)(a[i].x+0.5)<<" \n"[i==n+m];
	return 0;
}
```

## $\mathfrak{Problemset}$

## $\mathfrak{Part 1}$

#### P3338 [ZJOI2014]力

$$
F_i=q_i\left(\sum_{j=1}^i \dfrac{q_j}{(i-j)^2} -\sum_{j=i}^n \dfrac{q_j}{(i-j)^2} \right)\\
E_i=\dfrac{F_i}{q_i}
$$

令 $g(x)=\dfrac{1}{x^2}$，则有：
$$
E_i=\sum_{j=1}^i q_jg_{i-j} - \sum_{j=i}^n q_j g_{j-i}
$$
前面的式子已经可以卷积了，考虑分析后面的式子：
$$
\sum_{j=i}^n q_j g_{j-i}\\
=\sum_{j=0}^{n-i} q_{j+i}g_{j}
$$
令 $f_i=q_{n-i}$，则有：
$$
\sum_{j=0}^{n-i} f_{n-(i+j)} g_{j}\\
=\sum_{j=0}^{n-i} f_{(n-i)-j} g_{j}\\
=\sum_{j=0}^{k=n-i} f_{k-j} g_j
$$
然后分别计算一下卷积即可，即将 $g$ 与 $f$、$q$ 分别做一边多项式乘法。

$Code:$

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define db double
#define ll long long
#define N 500005
const db pi=std::acos(-1);
struct poly{
	db x,y;
	inline poly operator+(const poly&p)const{return poly{x+p.x,y+p.y};}
	inline poly operator-(const poly&p)const{return poly{x-p.x,y-p.y};}
	inline poly operator*(const poly&p)const{return poly{x*p.x-y*p.y,x*p.y+y*p.x};}
} q[N],g[N],f[N],cpy[N];

namespace fft{
	int lim,rev[N];
	inline void install(int n){
		int len=0; lim=1;
		while(lim<=n) lim<<=1,++len;
		rep(i,0,lim-1) rev[i]=(rev[i>>1]>>1)|((i&1)<<(len-1));
	}
	inline void fft(poly *a,int ty){
		rep(i,0,lim-1) if(i<rev[i]) std::swap(a[i],a[rev[i]]);
		for(int mid=1;mid<lim;mid<<=1){
			poly wn=poly{std::cos(pi/mid),ty*std::sin(pi/mid)};
			for(int i=mid<<1,j=0;j<lim;j+=i){
				poly w=poly{1,0};
				for(int k=0;k<mid;++k,w=w*wn){
					poly x=a[j+k],y=w*a[j+mid+k];
					a[j+k]=x+y,a[j+mid+k]=x-y;
				}
			}
		}
		if(!~ty){
			rep(i,0,lim) a[i].x/=lim;
		}
	}
	inline void _fft(poly *a,poly *b){
		fft(a,1),fft(b,1);
		rep(i,0,lim) a[i]=a[i]*b[i];
		fft(a,-1);
	}
}
int n;
int main(){
	scanf("%d",&n);
	rep(i,1,n) scanf("%lf",&q[i].x);
	fft::install(n+n);
	rep(i,1,n) g[i].x=1.0/i/i;
	rep(i,1,n) f[n-i]=q[i];
	using fft::fft;
	fft(g,1),fft(q,1),fft(f,1);
	rep(i,0,fft::lim-1) q[i]=q[i]*g[i],f[i]=f[i]*g[i];
	fft(q,-1),fft(f,-1);
	rep(i,1,n) printf("%.6lf\n",q[i].x-f[n-i].x);
	return 0;
}
```

#### P3723 [AH2017/HNOI2017]礼物

考虑一下：增加不同手环的亮度和手环循环等价于最小化：
$$
\sum_{i=1}^n (x_i-y_i+c)^2
$$
其中 $c$ 为增加的亮度，$y_i$ 为原先 $y$ 的一个循环数组。

考虑展开：
$$
\sum_{i=1}^n (x_i-y_i+c)^2\\
=\sum_{i=1}^n (x_i^2+y_i^2+c^2-2x_iy_i-2y_ic+2x_ic)\\
=\sum_{i-1}^n (x_i^2+y_i^2) + \left(nc^2+2c\sum_{i=1}^n (x_i-y_i)\right) -2\sum_{i=1}^n x_i y_i
$$
第一项值为定值，第二项直接二次函数求顶点，考虑转化第三项。

我们注意到 $x_i y_i$ 比较像一个卷积的形式，于是把 $y_i$ 翻转，然后做一遍卷积。

我们注意到 $y_i$ 是一个循环数组的形式，于是考虑可以倍长 $y_i$。

于是作一遍卷积即可。

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define ll long long
#define db double
#define N 5000005
const int mod=998244353,G=3,_inv=mod+1>>1;
inline void reduce(int&x){x-=mod,x+=x>>31&mod;}
//inline int upd(int x){x-=mod; return x+=x>>31&mod;}
inline int upd(int x){return x+=x>>31&mod;}
int wn[N],rev[N];
inline int qpow(int x,int y,int ret=1){
	for(;y;y>>=1,x=1ll*x*x%mod) if(y&1)
		ret=1ll*ret*x%mod; return ret;
}
namespace ntt{
	inline void install(int n){
		for(int mid=1;mid<n;mid<<=1){
			int w=qpow(G,mod/mid>>1);
			wn[mid]=1;
			rep(i,1,mid-1) wn[mid+i]=1ll*wn[mid+i-1]*w%mod;
		}
	}
	int lim,invlim;
	inline void init(int n){
		for(lim=invlim=1;lim<n;lim<<=1)
			invlim=1ll*invlim*_inv%mod;
		rep(i,1,lim-1) rev[i]=rev[i>>1]>>1|(lim>>1&-(i&1));
	}
	inline void ntt(int *a,int ty){
		static int b[N];
		rep(i,0,lim-1) b[i]=a[rev[i]];
		for(int mid=1;mid<lim;mid<<=1){
			for(int i=0;i<lim;i+=mid<<1){
				rep(j,0,mid-1){
					const int x=1ll*b[i+j+mid]*wn[mid+j]%mod;
					b[i+j+mid]=upd(b[i+j]-x);
					reduce(b[i+j]+=x);
				}
			}
		}
		if(!ty){
			rep(i,0,lim-1) a[i]=1ll*b[i]*invlim%mod;
			std::reverse(a+1,a+lim);
		} else {
			std::memcpy(a,b,lim<<2);
		}
	}
	inline void _ntt(int *a,int *b){
		ntt(a,1),ntt(b,1);
		rep(i,0,lim-1) a[i]=1ll*a[i]*b[i]%mod;
		ntt(a,0);
	}
}
int x[N],y[N],n,m;
int main(){
	std::ios::sync_with_stdio(0),std::cin.tie(0),std::cout.tie(0);
	std::cin>>n>>m;
	rep(i,1,n) std::cin>>x[i];
	rep(i,1,n) std::cin>>y[i];
	ll ans=0; int b=0;
	rep(i,1,n) ans+=x[i]*x[i]+y[i]*y[i],b+=x[i]-y[i];
	int l=floor(-b/(db)n),r=ceil(-b/(db)n);
	ans+=std::min(1ll*n*l*l+2ll*l*b,1ll*n*r*r+2ll*r*b);
	std::cerr<<ans<<'\n';
	std::reverse(y+1,y+n+1);
	rep(i,n+1,2*n) y[i]=y[i-n];
	ntt::install((n+1)*6),ntt::init((n+1)*3);
	ntt::_ntt(x,y);
	int tmp=0;
	rep(i,1,n) tmp=std::max(tmp,x[i+n]);
	std::cout<<(ans-tmp*2)<<'\n';
	return 0;
}
```

#### P5641 【CSGRound2】开拓者的卓识

我们考虑计算贡献。

考虑最终的对于一个数 $a_i$，设它在 $k$ 个选定多区间嵌套内的出现次数为 $c_i$，则它对答案的贡献为 $a_ic_i$。

那么现在我们考虑一个区间最终区间 $[1,r]$，它的其中答案来源 $a_i$，来源于若干 $k$ 个区间 $[l_i,r_i]$ 满足 $l_k=1\leqslant l_{k-1}\leqslant l_{k-2}\leqslant \ldots\leqslant l_1\leqslant i\leqslant r_1\leqslant \ldots\leqslant r_{k-1}\leqslant r_k$。

那么现在考虑我们从每一个 $i=1\ldots r$ 开始走，不断向左边或者右边走一步或者不走，这样容易发现左端点 $l$ 和右端点 $r$ 的贡献是互相独立的，因此我们分开考虑。

我们先考虑左端点。

其实就是插板法，左边一共可以走 $k-1$ 次，答案为 $\binom{(i+(k-1))-1}{k-1}=\binom{i+k-2}{k-1}$。

然后右边同理，答案为 $\binom{((r-i+1)+(k-1))-1}{k-1}=\binom{r-i+k-1}{k-1}$。

答案即为
$$
sum_{k,1,r}=\sum_{i=1}^r a_r \binom{i+k-2}{k-1} \binom{r-i+k-1}{k-1}
$$
我们考虑优化。

显然我们需要 NTT 来优化。

令 $p(x)=\binom{x}{k-1}$。

则
$$
sum_{k,1,r}=\sum_{i=1}^r a_i p(i+k-2) p(r-i+k-1)
$$
令 $f(i)=a_ip(i+k-2),g(i)=p(i+k-1)$，则
$$
sum_{k,1,r}=\sum_{i=1}^r f(i) g(r-i)
$$
卷一下即可。

对于组合数部分的递推：
$$
\binom{i+k-2}{k-1}=\dfrac{(i+k-2)!}{(k-1)!(i+k-2-k+1)!}\\
=\dfrac{(i+k-2)!}{(k-1)!(i-1)!}\\
=\dfrac{k\times (k+1)\times \ldots\times (i+k-2)}{(i-1)!}
$$
于是这部分递推即可。
$$
\binom{i+k-1}{k-1}=\dfrac{k\times (k+1)\times \ldots \times (i+k-1)}{i!}
$$
于是预处理到 $n$ 的阶乘和从 $k$ 开始的乘积即可（或者直接递推也可以）。

时间复杂度：$O(n\log n)$。

$Code:$

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define ll long long
#define N 5000005
const int mod=998244353,G=3,_inv=mod+1>>1;
inline void reduce(int&x){x-=mod,x+=x>>31&mod;}
//inline int upd(int x){x-=mod; return x+=x>>31&mod;}
inline int upd(int x){return x+=x>>31&mod;}
int wn[N],rev[N];
inline int qpow(int x,int y,int ret=1){
	for(;y;y>>=1,x=1ll*x*x%mod) if(y&1)
		ret=1ll*ret*x%mod; return ret;
}
namespace ntt{
	inline void install(int n){
		for(int mid=1;mid<n;mid<<=1){
			int w=qpow(G,mod/mid>>1);
			wn[mid]=1;
			rep(i,1,mid-1) wn[mid+i]=1ll*wn[mid+i-1]*w%mod;
		}
	}
	int lim,invlim;
	inline void init(int n){
		for(lim=invlim=1;lim<n;lim<<=1)
			invlim=1ll*invlim*_inv%mod;
		rep(i,1,lim-1) rev[i]=rev[i>>1]>>1|(lim>>1&-(i&1));
	}
	inline void ntt(int *a,int ty){
		static int b[N];
		rep(i,0,lim-1) b[i]=a[rev[i]];
		for(int mid=1;mid<lim;mid<<=1){
			for(int i=0;i<lim;i+=mid<<1){
				rep(j,0,mid-1){
					const int x=1ll*b[i+j+mid]*wn[mid+j]%mod;
					b[i+j+mid]=upd(b[i+j]-x);
					reduce(b[i+j]+=x);
				}
			}
		}
		if(!ty){
			rep(i,0,lim-1) a[i]=1ll*b[i]*invlim%mod;
			std::reverse(a+1,a+lim);
		} else {
			std::memcpy(a,b,lim<<2);
		}
	}
	inline void _ntt(int *a,int *b){
		ntt(a,1),ntt(b,1);
		rep(i,0,lim-1) a[i]=1ll*a[i]*b[i]%mod;
		ntt(a,0);
	}
}
int a[N],n,k,inv[N];
int f[N],g[N];
int main(){
	std::ios::sync_with_stdio(0),std::cin.tie(0),std::cout.tie(0);
	std::cin>>n>>k;
	rep(i,1,n) std::cin>>a[i];
	inv[1]=1;
	rep(i,2,n) inv[i]=1ll*inv[mod%i]*(mod-mod/i)%mod;
	g[0]=1;
	rep(i,1,n) g[i]=1ll*g[i-1]*(i+k-1)%mod*inv[i]%mod;
	rep(i,1,n) f[i]=1ll*g[i-1]*a[i]%mod;
	ntt::install(n+1<<2); ntt::init(n+1<<1);
	ntt::_ntt(f,g);
	rep(i,1,n) std::cout<<f[i]<<" \n"[i==n];
	return 0;
}
```

#### P4245 【模板】任意模数多项式乘法

crt 即可。最后合并答案。

不过似乎常数巨大？~~不管了不管了~~

$Code:$

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define ll long long
#define N 500005
inline int qpow(int x,int y,int mod,int ret=1){
	for(;y;y>>=1,x=1ll*x*x%mod) if(y&1)
		ret=1ll*ret*x%mod; return ret;
}
int mod;
const int mod1=998244353,mod2=1004535809,mod3=469762049,G=3;
const ll mod_1_2=1ll*mod1*mod2;
const int inv_1=qpow(mod1,mod2-2,mod2),inv_2=qpow(mod_1_2%mod3,mod3-2,mod3);

struct crt{
	int a,b,c;
	crt(){}
	crt(int _a){a=b=c=_a;}
	crt(int _a,int _b,int _c){a=_a,b=_b,c=_c;}
	inline friend crt upd(const crt&p){
		return crt(p.a+(p.a>>31&mod1),p.b+(p.b>>31&mod2),p.c+(p.c>>31&mod3));
	}
	inline crt operator+(const crt&p)const{
		return upd(crt(a+p.a-mod1,b+p.b-mod2,c+p.c-mod3));
	}
	inline crt operator-(const crt&p)const{
		return upd(crt(a-p.a,b-p.b,c-p.c));
	}
	inline crt operator*(const crt&p)const{
		return crt(1ll*a*p.a%mod1,1ll*b*p.b%mod2,1ll*c*p.c%mod3);
	}
	inline int norm(){
		const ll x=1ll*(b-a+mod2)%mod2*inv_1%mod2*mod1+a;
		return (1ll*(c-x%mod3+mod3)%mod3*inv_2%mod3*(mod_1_2%mod)%mod+x)%mod;
	}
} wn[N];
int rev[N];
namespace ntt{
	inline void install(int n){
		for(int mid=1;mid<n;mid<<=1){
			crt w=crt(qpow(G,mod1/mid>>1,mod1),qpow(G,mod2/mid>>1,mod2),qpow(G,mod3/mid>>1,mod3));
			wn[mid]=crt(1);
			rep(i,1,mid-1) wn[mid+i]=wn[mid+i-1]*w;
		}
	}
	int lim;
	inline void init(int n){
		for(lim=1;lim<n;lim<<=1);
		rep(i,1,lim-1) rev[i]=rev[i>>1]>>1|(lim>>1&-(i&1));
	}
	inline void ntt(crt *a,int ty){
		static crt b[N];
		rep(i,0,lim-1) b[i]=a[rev[i]];
		for(int mid=1;mid<lim;mid<<=1){
			for(int i=0;i<lim;i+=mid<<1){
				rep(j,0,mid-1){
					const crt x=b[i+j+mid]*wn[mid+j];
					b[i+j+mid]=b[i+j]-x;
					b[i+j]=b[i+j]+x;
				}
			}
		}
		if(!ty){
			crt invlim=crt(qpow(lim,mod1-2,mod1),qpow(lim,mod2-2,mod2),qpow(lim,mod3-2,mod3));
			rep(i,0,lim-1) a[i]=b[i]*invlim;
			std::reverse(a+1,a+lim);
		} else {
			//rep(i,0,lim-1) a[i]=b[i];
			std::memcpy(a,b,lim*12);
		}
	}
	inline void _ntt(crt *a,crt *b){
		ntt(a,1),ntt(b,1);
		rep(i,0,lim-1) a[i]=a[i]*b[i];
		ntt(a,0);
	}
}

crt a[N],b[N];
int n,m;
int main(){
	std::ios::sync_with_stdio(0),std::cin.tie(0),std::cout.tie(0);
	std::cin>>n>>m>>mod;
	rep(i,0,n){
		int x; std::cin>>x;
		a[i]=crt(x%mod);
	}
	rep(i,0,m){
		int x; std::cin>>x;
		b[i]=crt(x%mod);
	}
	ntt::install(n+m+1<<1),ntt::init(n+m+1);
	ntt::_ntt(a,b);
	rep(i,0,n+m) std::cout<<a[i].norm()<<" \n"[i==n+m];
	return 0;
}
```

#### P5488 差分与前缀和

前缀和：

考虑它的具体意义。

同开拓者的卓识，我们令 $c_i$ 为 $a_i$ 对 $k$ 维前缀和时某一位的贡献。假设那一位为 $n$。

那么转化成开拓者的卓识不难得到 $c_n=\binom{n+k-1}{k}$。

差分：

考虑差分的定义式：$b_i=a_i-a_{i-1}$。

假设某一位为 $n$，则对其贡献为 $c_n=\binom{k}{n}\times (-1)^n$。

于是递推求出 $\{c\}$，然后卷一下即可。

$Code:$

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define ll long long
#define N 5000005
const int mod=1004535809,G=3,_inv=mod+1>>1;
inline void reduce(int&x){x-=mod,x+=x>>31&mod;}
//inline int upd(int x){x-=mod; return x+=x>>31&mod;}
inline int upd(int x){return x+=x>>31&mod;}
int wn[N],rev[N];
inline int qpow(int x,int y,int ret=1){
	for(;y;y>>=1,x=1ll*x*x%mod) if(y&1)
		ret=1ll*ret*x%mod; return ret;
}
namespace ntt{
	inline void install(int n){
		for(int mid=1;mid<n;mid<<=1){
			int w=qpow(G,mod/mid>>1);
			wn[mid]=1;
			rep(i,1,mid-1) wn[mid+i]=1ll*wn[mid+i-1]*w%mod;
		}
	}
	int lim,invlim;
	inline void init(int n){
		for(lim=invlim=1;lim<n;lim<<=1)
			invlim=1ll*invlim*_inv%mod;
		rep(i,1,lim-1) rev[i]=rev[i>>1]>>1|(lim>>1&-(i&1));
	}
	inline void ntt(int *a,int ty){
		static int b[N];
		rep(i,0,lim-1) b[i]=a[rev[i]];
		for(int mid=1;mid<lim;mid<<=1){
			for(int i=0;i<lim;i+=mid<<1){
				rep(j,0,mid-1){
					const int x=1ll*b[i+j+mid]*wn[mid+j]%mod;
					b[i+j+mid]=upd(b[i+j]-x);
					reduce(b[i+j]+=x);
				}
			}
		}
		if(!ty){
			rep(i,0,lim-1) a[i]=1ll*b[i]*invlim%mod;
			std::reverse(a+1,a+lim);
		} else {
			std::memcpy(a,b,lim<<2);
		}
	}
	inline void _ntt(int *a,int *b){
		ntt(a,1),ntt(b,1);
		rep(i,0,lim-1) a[i]=1ll*a[i]*b[i]%mod;
		ntt(a,0);
	}
}
int a[N],b[N],n,k,type;

inline int rd(){
	int x=0; char c;
	for(;!isdigit(c=getchar()););
	for(x=c&15;isdigit(c=getchar());x=(x*10ll+(c&15))%mod);
	return x;
}

int main(){
//	std::ios::sync_with_stdio(0),std::cin.tie(0),std::cout.tie(0);
	n=rd();
	k=rd();
	std::cerr<<k<<'\n';
	type=rd();
	b[0]=1;
	rep(i,0,n-1) a[i]=rd();
	if(type==0){
		rep(i,1,n-1) b[i]=1ll*b[i-1]*qpow(i,mod-2)%mod*(k+i-1ll)%mod;
	} else {
		rep(i,1,n-1) b[i]=1ll*b[i-1]*qpow(i,mod-2)%mod*(k-i+1ll)%mod;
		rep(i,1,n-1) if(i&1) b[i]=upd(-b[i]);
	}
	ntt::install(n<<2),ntt::init(n<<1);
	ntt::_ntt(a,b);
	rep(i,0,n-1) std::cout<<a[i]<<" \n"[i==n-1];
	return 0;
}
```

### $\mathfrak{Part 2}$

#### P4238 【模板】多项式乘法逆

设多项式 $f$ 的逆元为 $g$，那么有 $f *g\equiv 1\pmod{x^n} $。

显然如果 $f$ 只有一个常数项 $f_0$（这里令 $f_i=[x^i]f$），那么$f_0 * g_0 \equiv 1\pmod{x^0}$。即 $g_0$ 为 $f_0$ 的逆元。

那么加假设我们已经知道了 $f * h\equiv 1\pmod{x^{\lceil \frac{n}{2} \rceil}}$，（显然 $f * g\equiv 1\pmod{x^{\lceil \frac{n}{2} \rceil}}$）那么有：
$$
\begin{cases} f * g\equiv 1\pmod{x^{\lceil \frac{n}{2} \rceil}} \\ f * h\equiv 1\pmod{x^{\lceil \frac{n}{2} \rceil}} \end{cases} \Longleftrightarrow f * (g-h)\equiv 0\pmod{x^{\lceil \frac{n}{2} \rceil}}\\
\Longrightarrow g-h\equiv 0\pmod{x^{\lceil \frac{n}{2} \rceil}}\\
\Longrightarrow (g-h)^2\equiv 0\pmod{x^n}\\
\Longrightarrow g^2 -2g * h+h^2 \equiv 0\pmod{x^n}\\
\Longrightarrow f * (g^2 -2g * h+h^2) \equiv 0\pmod{x^n}\\
\Longrightarrow g-2h+f * h^2\equiv 0\pmod{x^n}\\
\Longrightarrow g\equiv 2h-f * h^2 \pmod{x^2}
$$
于是递归求解即可。

根据主定理，时间复杂度为 $O(n\log n)$。

#### P4725 【模板】多项式对数函数（多项式 ln）

> 前置知识：复合函数 $f(g(x))$ 的求导公式为 $(f(g(x)))'=f'(g(x))g'(x)$。$\ln(x)$ 的导数为 $\dfrac{1}{x}$。

令  $f(x)=A(x),g(x)=B(x)$。

则
$$
g(x)\equiv  \ln(f(x))\pmod{x^n}\\
g'(x)=(\ln(f(x)))'\pmod{x^n}\\
g'(x)=\ln'(f(x))f'(x)\pmod{x^n}\\
g'(x)=\dfrac{f'(x)}{f(x)}\pmod{x^n}
$$
于是先求出 $f$ 的导数，然后求出 $f$ 的逆，再作一边卷积后求出 $g'$ 的原函数即可。

```cpp
#include<bits/stdc++.h>
#define rep(i,x,y) for(int i=x,i##end=y;i<=i##end;++i)
#define _rep(i,x,y) for(int i=x,i##end=y;i>=i##end;--i)
#define ll long long
#define N 5000005
const int mod=998244353,G=3,_inv=mod+1>>1;
inline void reduce(int&x){x-=mod,x+=x>>31&mod;}
//inline int upd(int x){x-=mod; return x+=x>>31&mod;}
inline int upd(int x){return x+=x>>31&mod;}
int wn[N],rev[N];
inline int qpow(int x,int y,int ret=1){
	for(;y;y>>=1,x=1ll*x*x%mod) if(y&1)
		ret=1ll*ret*x%mod; return ret;
}
namespace ntt{
	inline void install(int n){
		for(int mid=1;mid<n;mid<<=1){
			int w=qpow(G,mod/mid>>1);
			wn[mid]=1;
			rep(i,1,mid-1) wn[mid+i]=1ll*wn[mid+i-1]*w%mod;
		}
	}
	int lim,invlim;
	inline void init(int n){
		for(lim=invlim=1;lim<n;lim<<=1)
			invlim=1ll*invlim*_inv%mod;
		rep(i,1,lim-1) rev[i]=rev[i>>1]>>1|(lim>>1&-(i&1));
	}
	inline void ntt(int *a,int ty){
		static int b[N];
		rep(i,0,lim-1) b[i]=a[rev[i]];
		for(int mid=1;mid<lim;mid<<=1){
			for(int i=0;i<lim;i+=mid<<1){
				rep(j,0,mid-1){
					const int x=1ll*b[i+j+mid]*wn[mid+j]%mod;
					b[i+j+mid]=upd(b[i+j]-x);
					reduce(b[i+j]+=x);
				}
			}
		}
		if(!ty){
			rep(i,0,lim-1) a[i]=1ll*b[i]*invlim%mod;
			std::reverse(a+1,a+lim);
		} else {
			std::memcpy(a,b,lim<<2);
		}
	}
	inline void _ntt(int *a,int *b){
		ntt(a,1),ntt(b,1);
		rep(i,0,lim-1) a[i]=1ll*a[i]*b[i]%mod;
		ntt(a,0);
	}
}
int f[N],g[N],n;

inline void change(int *a,int *b,int n,int ty){
	if(ty==0){
		rep(i,1,n-1) b[i-1]=1ll*a[i]*i%mod;
	} else {
		rep(i,1,n-1) b[i]=1ll*a[i-1]*qpow(i,mod-2)%mod;
	}
}

inline void solve(int len,int *f,int *g){
	if(len==1) return g[0]=qpow(f[0],mod-2),void();
	solve(len+1>>1,f,g);
	ntt::install(len<<2),ntt::init(len<<1);
	static int c[N];
	rep(i,0,len-1) c[i]=f[i];
	rep(i,len,ntt::lim-1) c[i]=0;
	ntt::ntt(g,1),ntt::ntt(c,1);
	rep(i,0,ntt::lim-1) g[i]=upd(2ll*g[i]%mod-1ll*c[i]*g[i]%mod*g[i]%mod);
	ntt::ntt(g,0);
	rep(i,len,ntt::lim-1) g[i]=0;
}
int a[N],b[N];
inline void getln(int *f,int *g){
	change(f,a,n,0);
	solve(n,f,b);
	ntt::install(n<<2),ntt::init(n<<1);
	ntt::_ntt(a,b);
	change(a,g,n,1);
}

int main(){
	std::ios::sync_with_stdio(0),std::cin.tie(0),std::cout.tie(0);
	std::cin>>n;
	rep(i,0,n-1) std::cin>>f[i];
	getln(f,g);
	rep(i,0,n-1) std::cout<<g[i]<<" \n"[i==n-1];
	return 0;
}
```

