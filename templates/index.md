[$\color{purple}{\text{back to top}}$](https://cyn2006.github.io/main.html)

<div>
    <font size="5",style="font-family:SontTi">
    	istream
	</font>
</div>


```cpp
struct istream{
	char *nxt,buf[1<<25];
	istream&init(FILE *f=stdin){
		fread(buf,1,sizeof(buf),f);
		return nxt=buf,*this;
	}
	template<typename T>inline istream&operator>>(T&x){
		x=0; int f=0; char c;
		for(;!isdigit(c=*nxt++);f=c) if(!~c) return *this;
		for(x=c&15;isdigit(c=*nxt++);x=x*10+(c&15));
		return (f==45&&(x=-x)),*this;
	}
	inline istream&operator>>(char&c){
		for(;!isalpha(c=*nxt++););
		return *this;
	}
} cin;
```

<div>
	<font size="5",style="font-family:SontTi">
        ostream
    </font>
</div>

```cpp
struct ostream{
	char *nxt,buf[1<<23],_buf[34];
	ostream(){nxt=buf;}
	inline void flush(FILE *f=stdout){
		fwrite(buf,1,nxt-buf,f),nxt=buf;
	}
	template<typename T>inline ostream&operator<<(T x){
		if(!x) return pc(48),*this;
		if(x<0) x=-x,pc(48);
		int cnt=0;
		for(;x;x/=10) _buf[++cnt]=x%10|48;
		for(;cnt;--cnt) pc(_buf[cnt]);
		return *this;
	}
	inline ostream&operator<<(char c){return pc(c),*this;}
	inline void pc(char c){*nxt++=c;}
} cout;
```

<div>
	<font size="5",style="font-family:SontTi">
        并查集
    </font>
</div>

```cpp
inline int ancestor(int x){
    return x==fa[x]?x:fa[x]=ancestor(fa[x]);
}
inline void merge(int x,int y){
    x=ancestor(x),y=ancestor(y);
    if(x==y) return;
    if(sz[x]>sz[y]) fa[y]=x,sz[x]+=sz[y];
    else fa[x]=y,sz[y]+=sz[x];
}
```

<div>
	<font size="5",style="font-family:SontTi">
        存图建边
    </font>
</div>

```cpp
std::vector<int> e[N];
// std::vector<pii> e[N];
inline void add(int u,int v){
    e[u].push_back(v);
	// e[u].emplace_back(v,w);
}
```

<div>
	<font size="5",style="font-family:SontTi">
        ODT
    </font>
</div>

```cpp
namespace ODT{
struct node{
	int l,r; mutable int v;
	node(int L,int R=-1,int V=0):l(L),r(R),v(V){}
	inline bool operator<(const node&x)const{return l<x.l;}
};
#define It set<node>::iterator
set<node> S;
inline It split(int p){
	It it=S.lower_bound(node(p));
	if(it!=S.end()&&it->l==p) return it; --it;
	int L=it->l,R=it->r,val=it->v;
	S.erase(it),S.insert(node(L,p-1,val));
	return S.insert(node(p,R,val)).fi;
}
inline void assign(int l,int r,int v){
	It R=split(r+1),L=split(l); S.erase(L,R),S.insert(node(l,r,v));
}
inline void query(int l,int r){
// Queries
}
}
```

<div>
	<font size="5",style="font-family:SontTi">
        树状数组
    </font>
</div>

```cpp
struct Array_tree{
int c[N],n;
inline void add(int x,int v){for(;x<=n;x+=x&-x) c[x]+=v;}
inline int qry(int x,int ret=0){for(;x;x-=x&-x) ret+=c[x]; return ret;}
} T;
```

<div>
	<font size="5",style="font-family:SontTi">
        Trie
    </font>
</div>

```cpp
struct Trie{
	int t[N*32][2],tot,cnt[N*32];
	void clear(){
		memset(t,0,tot+1<<3);
		memset(cnt,0,tot+1<<3);
	}
	void insert(const int x){
		register int now=0;
		_rep(i,20,0){
			int tmp=x>>i&1;
			if(!t[now][tmp]) t[now][tmp]=++tot;
			now=t[now][tmp],++cnt[now];
		}
	}
	void erase(const int x){
		register int now=0;
		_rep(i,20,0){
			int tmp=x>>i&1;
			now=t[now][tmp],--cnt[now];
		}
	}
	int qry(const int x){
		register int ans=0,now=0;
		_rep(i,20,0){
			int tmp=x>>i&1;
			if(t[now][tmp^1]&&cnt[t[now][tmp^1]]>0)
				ans|=1<<i,now=t[now][tmp^1];
			else now=t[now][tmp];
		}
		return ans;
	}
} T;
```

<div>
	<font size="5",style="font-family:SontTi">
        NTT
    </font>
</div>

```cpp
namespace NTT{
const int N=5e6+5,mod=998244353,G=3,_inv=mod+1>>1;
inline int qmod(const int x){return x>=mod?x-mod:x;}
inline void reduce(int&x){x>=mod?x-=mod:0;}
using rp=int[N];
rp wn,rev;
int n,m,a[N],b[N];
#define ll long long
inline int qpow(int x,int y,int ret=1){
	for(;y;y>>=1,x=1ll*x*x%mod) if(y&1)
		ret=1ll*ret*x%mod; return ret; 
}
void initall(int n){
	for(int mid=1;mid<n;mid<<=1){
		int w=qpow(G,mod/mid>>1);
		wn[mid]=1;
		for(int i=1;i<mid;++i)
			wn[mid+i]=1ll*wn[mid+i-1]*w%mod;
	}
}
int lim,invlim;
void init(int n){
	for(lim=invlim=1;lim<n;lim<<=1)
		invlim=1ll*invlim*_inv%mod;
	for(int i=1;i<lim;++i)
		rev[i]=rev[i>>1]>>1|(lim>>1&-(i&1));
}
inline void fft(int*a,int ty){
	static int b[N];
	for(int i=0;i<lim;++i) b[i]=a[rev[i]];
	for(int mid=1;mid<lim;mid<<=1){
		for(int i=0;i<lim;i+=mid<<1){
			for(int j=0;j<mid;++j){
				register ll x=1ll*b[i+j+mid]*wn[mid+j]%mod;
				b[i+j+mid]=qmod(b[i+j]+mod-x),reduce(b[i+j]+=x);
			}
		}
	}
	if(!ty){
		for(int i=0;i<lim;++i) a[i]=1ll*b[i]*invlim%mod;
		reverse(a+1,a+lim);
	} else {
		for(int i=0;i<lim;++i) a[i]=b[i];
	}
}

inline void main_solve(){
	ios::sync_with_stdio(0),cin.tie(0),cout.tie(0);
	cin>>n>>m;
	for(int i=0;i<=n;++i) cin>>a[i];
	for(int i=0;i<=m;++i) cin>>b[i];
	initall((n+m+2)<<1),init(n+m+1);
	fft(a,1),fft(b,1);
	for(int i=0;i<lim;++i) a[i]=1ll*a[i]*b[i]%mod;
	fft(a,0);
	for(int i=0;i<=n+m;++i) cout<<a[i]<<' '; cout<<'\n';
}

}
```

