#include "../common.inc"
#include <climits>
// instrumented bidir (copy of solver's) to count work
static vector<uint32_t> fh,bh; static vector<int32_t> fto,bto; static vector<uint32_t> fwt,bwt;
static void csr(const vector<int>&a,const vector<int>&b,vector<uint32_t>&h,vector<int>&to,vector<uint32_t>&w){
  h.assign(V+1,0); for(int i=0;i<E;i++)h[a[i]+1]++; for(int v=0;v<V;v++)h[v+1]+=h[v];
  vector<std::pair<uint32_t,int>> t(E); vector<uint32_t> c(h.begin(),h.end()-1);
  for(int i=0;i<E;i++)t[c[a[i]]++]={ew[i],b[i]}; for(int v=0;v<V;v++)std::sort(t.begin()+h[v],t.begin()+h[v+1]);
  to.resize(E);w.resize(E); for(int i=0;i<E;i++){w[i]=t[i].first;to[i]=t[i].second;}
}
int main(int argc,char**argv){
  read_graph(argv[1]); read_queries(argv[2]);
  csr(eu,ev,fh,fto,fwt); csr(ev,eu,bh,bto,bwt);
  const int64_t INF=LLONG_MAX/4; vector<int64_t> dF(V,INF),dB(V,INF); vector<int> tF,tB; Heap4<int64_t> hF,hB;
  double sF=0,sB=0,rF=0,rB=0; int n=0;
  for(int q=0;q<Q;q+=5){int s=qs[q],t=qt[q]; if(s==t)continue; n++;
    int64_t mu=INF; hF.clear();hB.clear(); dF[s]=0;tF.push_back(s);hF.push(0,s); dB[t]=0;tB.push_back(t);hB.push(0,t);
    uint64_t wF=0,wB=0;
    while(!hF.empty()&&!hB.empty()){
      int64_t kf=hF.top_key(),kb=hB.top_key(); if(kf+kb>=mu)break;
      int tf=hF.a[0].v,tb=hB.a[0].v; bool fwd=wF+(fh[tf+1]-fh[tf])<=wB+(bh[tb+1]-bh[tb]);
      auto&hp=fwd?hF:hB; auto&d1=fwd?dF:dB; auto&d2=fwd?dB:dF; auto&tc=fwd?tF:tB;
      const uint32_t*H=fwd?fh.data():bh.data(); const int*TO=fwd?fto.data():bto.data(); const uint32_t*WT=fwd?fwt.data():bwt.data();
      int64_t other=fwd?kb:kf; auto it=hp.pop(); int u=it.v; int64_t d=it.d; if(d>d1[u])continue;
      (fwd?sF:sB)++;
      uint32_t k=H[u],e=H[u+1];
      for(;k<e;k++){int64_t nd=d+WT[k]; if(nd>=mu)break; int x=TO[k]; if(nd<d1[x]){if(d1[x]==INF)tc.push_back(x); d1[x]=nd; if(d2[x]<INF&&nd+d2[x]<mu)mu=nd+d2[x]; if(nd+other<mu)hp.push(nd,x);}}
      (fwd?rF:rB)+=k-H[u]; (fwd?wF:wB)+=(k-H[u])+1;
    }
    for(int v:tF)dF[v]=INF; for(int v:tB)dB[v]=INF; tF.clear();tB.clear();
  }
  printf("n=%d settleF %.0f settleB %.0f relaxF %.0f relaxB %.0f\n",n,sF/n,sB/n,rF/n,rB/n);
}
