#include "../common.inc"
#include <climits>
// For K hubs (top in-degree): UB via hubs (brute force), then hub-avoiding bidir with mu=UB. Count work.
static vector<uint32_t> fh,bh; static vector<int32_t> fto,bto; static vector<uint32_t> fwt,bwt;
static void csr(const vector<int>&a,const vector<int>&b,vector<uint32_t>&h,vector<int>&to,vector<uint32_t>&w){
  h.assign(V+1,0); for(int i=0;i<E;i++)h[a[i]+1]++; for(int v=0;v<V;v++)h[v+1]+=h[v];
  vector<std::pair<uint32_t,int>> t(E); vector<uint32_t> c(h.begin(),h.end()-1);
  for(int i=0;i<E;i++)t[c[a[i]]++]={ew[i],b[i]}; for(int v=0;v<V;v++)std::sort(t.begin()+h[v],t.begin()+h[v+1]);
  to.resize(E);w.resize(E); for(int i=0;i<E;i++){w[i]=t[i].first;to[i]=t[i].second;}
}
const int64_t INF=LLONG_MAX/4;
static void full(int s,const vector<uint32_t>&h,const vector<int>&to,const vector<uint32_t>&w,vector<int64_t>&d){
  std::fill(d.begin(),d.end(),INF); Heap4<int64_t> H; d[s]=0; H.push(0,s);
  while(!H.empty()){auto it=H.pop(); if(it.d>d[it.v])continue; int u=it.v; for(uint32_t k=h[u];k<h[u+1];k++){int64_t nd=it.d+w[k]; if(nd<d[to[k]]){d[to[k]]=nd;H.push(nd,to[k]);}}}
}
int main(int argc,char**argv){
  read_graph(argv[1]); read_queries(argv[2]);
  csr(eu,ev,fh,fto,fwt); csr(ev,eu,bh,bto,bwt);
  vector<int> id(V); for(int i=0;i<E;i++)id[ev[i]]++;
  vector<int> ord(V); for(int v=0;v<V;v++)ord[v]=v; std::sort(ord.begin(),ord.end(),[&](int a,int b){return id[a]>id[b];});
  vector<int> rk(V); for(int k=0;k<V;k++)rk[ord[k]]=k;
  int Ks[]={0,30,100,300,1000,3000,10000}; const int NK=7; int NS=400;
  double sF[NK]={0},sB[NK]={0},rF[NK]={0},rB[NK]={0}; int exact[NK]={0}; int n=0;
  vector<int64_t> ds(V),dt(V),dF(V,INF),dB(V,INF); vector<int> tF,tB; Heap4<int64_t> hF,hB;
  for(int qi=0;qi<NS;qi++){int q=qi*97%Q; int s=qs[q],t=qt[q]; if(s==t)continue; n++;
    full(s,fh,fto,fwt,ds); full(t,bh,bto,bwt,dt);
    for(int ki=0;ki<NK;ki++){int K=Ks[ki];
      int64_t ub=INF; for(int k=0;k<K;k++){int h=ord[k]; if(ds[h]<INF&&dt[h]<INF) ub=std::min(ub,ds[h]+dt[h]);}
      int64_t truth=ds[t]; if(ub==truth)exact[ki]++;
      if(rk[s]<K||rk[t]<K) continue;
      int64_t mu=ub; hF.clear();hB.clear(); dF[s]=0;tF.push_back(s);hF.push(0,s); dB[t]=0;tB.push_back(t);hB.push(0,t);
      uint64_t wF=0,wB=0;
      while(!hF.empty()&&!hB.empty()){
        int64_t kf=hF.top_key(),kb=hB.top_key(); if(kf+kb>=mu)break;
        int tf=hF.a[0].v,tb=hB.a[0].v; bool fwd=wF+(fh[tf+1]-fh[tf])<=wB+(bh[tb+1]-bh[tb]);
        auto&hp=fwd?hF:hB; auto&d1=fwd?dF:dB; auto&d2=fwd?dB:dF; auto&tc=fwd?tF:tB;
        const uint32_t*H=fwd?fh.data():bh.data(); const int*TO=fwd?fto.data():bto.data(); const uint32_t*WT=fwd?fwt.data():bwt.data();
        int64_t other=fwd?kb:kf; auto it=hp.pop(); int u=it.v; int64_t d=it.d; if(d>d1[u])continue;
        (fwd?sF:sB)[ki]++;
        uint32_t k=H[u],e=H[u+1];
        for(;k<e;k++){int64_t nd=d+WT[k]; if(nd>=mu)break; int x=TO[k]; if(rk[x]<K)continue; if(nd<d1[x]){if(d1[x]==INF)tc.push_back(x); d1[x]=nd; if(d2[x]<INF&&nd+d2[x]<mu)mu=nd+d2[x]; if(nd+other<mu)hp.push(nd,x);}}
        (fwd?rF:rB)[ki]+=k-H[u]; (fwd?wF:wB)+=(k-H[u])+1;
      }
      for(int v:tF)dF[v]=INF; for(int v:tB)dB[v]=INF; tF.clear();tB.clear();
      if(mu!=truth&&!(mu>=INF/2&&truth>=INF/2)) printf("MISMATCH K=%d\n",K);
    }
  }
  for(int ki=0;ki<NK;ki++) printf("K=%5d exactUB %.3f settleF %.0f settleB %.0f relaxF %.0f relaxB %.0f\n",Ks[ki],(double)exact[ki]/n,sF[ki]/n,sB[ki]/n,rF[ki]/n,rB[ki]/n);
}
