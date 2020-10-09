[$\color{purple}{\text{back to top}}$](https://cyn2006.github.io/main/html)

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

