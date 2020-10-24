我太菜了。

## $\mathfrak{Problemset}$

### $\mathfrak{Part2}$

#### P5205 【模板】多项式开根

假设现在求出了 $f$ 模 $x^{\lceil \frac{n}{2} \rceil}$ 意义下的多项式 $g$ 满足 $g^2(x)\equiv f(x)\pmod{x^{\lceil \frac{n}{2} \rceil}}$。

令 $h$ 为 $f$ 模 $x^n$ 意义下的多项式，那么显然有
$$
\begin{cases} g^2(x)\equiv f(x)\pmod{x^{\lceil \frac{n}{2} \rceil}} \\h^2(x)\equiv f(x)\pmod{x^{\lceil \frac{n}{2} \rceil}} \end{cases}
$$
两式相减得到
$$
g(x)^2 \equiv h^2(x) \pmod{x^{\lceil \frac{n}{2} \rceil}}\\
g(x)^2 -h^2(x) \equiv 0\pmod{x^{\lceil \frac{n}{2} \rceil}}\\
(g(x)+h(x))(g(x)-h(x))\equiv 0\pmod{x^{\lceil \frac{n}{2} \rceil}}
$$
这里要求的是系数最小，那么我们取较小的正根，那么
$$
g(x)-h(x)\equiv 0\pmod{x^{\lceil \frac{n}{2} \rceil}}\\
g^2(x)-2h(x)g(x) +h^2(x) \equiv 0\pmod{x^n}\\
2h(x)g(x)=g^2(x)+f(x)\pmod{x^n}\\
h(x)=\dfrac{g^2(x)+f(x)}{2g(x)} \pmod{x^n}\\
h(x)=\dfrac{1}{2} \left(g(x) +\dfrac{f(x)}{g(x)}\right)
$$
然后每次将 $\dfrac{1}{g(x)}$ 和 $f(x)$ 分别转化成点值的形式，加起来即可。

#### P5277 【模板】多项式开根（加强版）

考虑常数项无法直接开根了。

然后考虑取模意义下的常数开根，即二次剩余。

于是 $\mathcal{BSGS}$ 或者 $\mathcal{Cipolla}$ 求出最小的二次剩余即可。

$Code:$（如果是弱化版，常数项为 $1$，把二次剩余部分去掉即可）

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
}
int f[N],g[N],n;

inline void solve(int len,int *f,int *g){
	if(len==1) return g[0]=qpow(f[0],mod-2),void();
	solve(len+1>>1,f,g);
//	ntt::install(len<<2);
	ntt::init(len<<1);
	static int c[N];
	rep(i,0,len-1) c[i]=f[i];
	rep(i,len,ntt::lim-1) c[i]=0;
	ntt::ntt(g,1),ntt::ntt(c,1);
	rep(i,0,ntt::lim-1) g[i]=upd(2ll*g[i]%mod-1ll*c[i]*g[i]%mod*g[i]%mod);
	ntt::ntt(g,0);
	rep(i,len,ntt::lim-1) g[i]=0;
}
int w;
struct complex{
	int x,y;
	complex(int _x=0,int _y=0){x=_x,y=_y;}
	inline complex operator*(const complex&p)const{
		return complex((1ll*x*p.x%mod+1ll*w*y%mod*p.y%mod)%mod,(1ll*x*p.y%mod+1ll*y*p.x%mod)%mod);
	}
};

inline int qpow(complex x,int y){
	complex ret=complex(1,0);
	for(;y;y>>=1,x=x*x) if(y&1) ret=ret*x;
	return ret.x;
}

inline int cipolla(int x){
	if(qpow(x,mod-1>>1)==mod-1) return -1;
	while(1){
		int a=(1ll*rand()<<15|rand())%mod;
		w=(1ll*a*a%mod-x+mod)%mod;
		if(qpow(w,mod-1>>1)==mod-1){
			int ret=qpow(complex(a,1),mod+1>>1);
			return std::min(ret,mod-ret);
		}
	}
}

inline void _sqrt(int len,int *f,int *g){
	if(len==1) return g[0]=cipolla(f[0]),std::cerr<<g[0]<<'\n',void();
	_sqrt(len+1>>1,f,g);
	static int c[N],d[N];
	memset(c,0,len<<3);
	solve(len,g,c);
//	std::cerr<<"c:\n";
//	rep(i,0,7) std::cerr<<g[i]<<" \n"[i==7];
	std::memcpy(d,f,len<<2);
//	ntt::install(len<<1);
	ntt::init(len<<1);
	rep(i,len,ntt::lim-1) d[i]=0;
	ntt::ntt(c,1),ntt::ntt(d,1);
	rep(i,0,ntt::lim-1) c[i]=1ll*c[i]*d[i]%mod;
	ntt::ntt(c,0);
	rep(i,0,len-1) g[i]=1ll*(g[i]+c[i])*_inv%mod;
//	rep(i,0,7) std::cerr<<g[i]<<" \n"[i==7];
}

int main(){
//	std::ios::sync_with_stdio(0),std::cin.tie(0),std::cout.tie(0);
	srand(time(NULL));
	std::cin>>n; ntt::install(n<<2);
	rep(i,0,n-1) std::cin>>f[i];
	_sqrt(n,f,g);
	rep(i,0,n-1) std::cout<<g[i]<<" \n"[i==n-1];
	return 0;
}
```

