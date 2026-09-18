#include "../common.inc"
#include <climits>
// access-set / uncovered-region sizes for top-K in-degree hubs
static vector<uint32_t> fh,bh; static vector<int32_t> fto,bto; static vector<uint32_t> fwt,bwt;
static void csr(const vector<int>&a,const vector<int>&b,vector<uint32_t>&h,vector<int>&to,vector<uint32_t>&w){
  h.assign(V+1,0); for(int i=0;i<E;i++)h[a[i]+1]++; for(int v=0;v<V;v++)h[v+1]+=h[v];
  vector<std::pair<uint32_t,int>> t(E); vector<uint32_t> c(h.begin(),h.end()-1);
  for(int i=0;i<E;i++)t[c[a[i]]++]={ew[i],b[i]}; for(int v=0;v<V;v++)std::sort(t.begin()+h[v],t.begin()+h[v+1]);
  to.resize(E);w.resize(E); for(int i=0;i<E;i++){w[i]=t[i].first;to[i]=t[i].second;}
}
const int64_t INF=LLONG_MAX/4;
vector<int> rk;
// returns (uncovered settled, uncovered hubs, settled until stop)
static void run(int s,int K,const vector<uint32_t>&h,const vector<int>&to,const vector<uint32_t>&w,double*acc){
  static vector<int64_t> d; static vector<char> cov; static vector<int> tch; if(d.empty()){d.assign(V,INF);cov.assign(V,0);}
  Heap4<int64_t> H; d[s]=0; cov[s]=0; tch.push_back(s); H.push(0,s); long open_unc=1; long unc=0,unch=0,sett=0;
  vector<char> done(0);
  static vector<char> st; if(st.empty())st.assign(V,0);
  while(!H.empty()&&open_unc>0){auto it=H.pop(); int u=it.v; if(it.d>d[u]||st[u])continue; st[u]=1; sett++;
    if(!cov[u]){open_unc--; unc++; if(rk[u]<K&&u!=s)unch++;}
    bool c=cov[u]||(rk[u]<K&&u!=s);
    for(uint32_t k=h[u];k<h[u+1];k++){int x=to[k]; int64_t nd=it.d+w[k]; if(st[x])continue;
      if(nd<d[x]){ if(d[x]==INF){tch.push_back(x); if(!c)open_unc++;} else { if(cov[x]&&!c)open_unc++; if(!cov[x]&&c)open_unc--; }
        d[x]=nd; cov[x]=c; H.push(nd,x);}
      else if(nd==d[x]&&c&&!cov[x]){cov[x]=1;open_unc--;}
    }
  }
  for(int x:tch){d[x]=INF;cov[x]=0;st[x]=0;} tch.clear();
  acc[0]+=unc; acc[1]+=unch; acc[2]+=sett;
}
int main(int argc,char**argv){
  read_graph(argv[1]); read_queries(argv[2]);
  csr(eu,ev,fh,fto,fwt); csr(ev,eu,bh,bto,bwt);
  vector<int> id(V),od(V); for(int i=0;i<E;i++){id[ev[i]]++;od[eu[i]]++;}
  vector<int> ord; for(int v=0;v<V;v++) if(od[v]>0) ord.push_back(v); std::stable_sort(ord.begin(),ord.end(),[&](int a,int b){return id[a]>id[b];});
  rk.assign(V,INT_MAX); for(size_t k=0;k<ord.size();k++)rk[ord[k]]=k;
  int Ks[]={300,1000,3000,10000};
  for(int K:Ks){ double fr[3]={0},br[3]={0},fh_[3]={0},bh_[3]={0}; int n=200;
    for(int i=0;i<n;i++){int y=qs[i*131%Q]; run(y,K,fh,fto,fwt,fr); run(y,K,bh,bto,bwt,br);
      int hb=ord[(i*7919)%K]; run(hb,K,fh,fto,fwt,fh_); run(hb,K,bh,bto,bwt,bh_);}
    printf("K=%d random v: fwd unc %.0f acc %.1f settled %.0f | bwd unc %.0f acc %.1f settled %.0f\n",K,fr[0]/n,fr[1]/n,fr[2]/n,br[0]/n,br[1]/n,br[2]/n);
    printf("K=%d hub   v: fwd unc %.0f acc %.1f settled %.0f | bwd unc %.0f acc %.1f settled %.0f\n",K,fh_[0]/n,fh_[1]/n,fh_[2]/n,bh_[0]/n,bh_[1]/n,bh_[2]/n);
    fflush(stdout);
  }
}
