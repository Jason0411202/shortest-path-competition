#include "../common.inc"
#include <map>
int main(int argc,char**argv){
  read_graph(argv[1]); read_queries(argv[2]);
  vector<int> od(V),id(V);
  for(int i=0;i<E;i++){od[eu[i]]++;id[ev[i]]++;}
  std::map<int,int> odh,idh; for(int v=0;v<V;v++){odh[od[v]]++; idh[std::min(id[v],20)]++;}
  printf("outdeg hist:"); for(auto&p:odh)printf(" %d:%d",p.first,p.second); printf("\n");
  printf("indeg hist(cap20):"); for(auto&p:idh)printf(" %d:%d",p.first,p.second); printf("\n");
  vector<int> ord(V); for(int v=0;v<V;v++)ord[v]=v;
  std::sort(ord.begin(),ord.end(),[&](int a,int b){return id[a]>id[b];});
  printf("top indeg:"); for(int k=0;k<20;k++)printf(" %d",id[ord[k]]); printf(" ... k300=%d k1000=%d k3000=%d\n",id[ord[299]],id[ord[999]],id[ord[2999]]);
  long long cum=0; int ks[]={1,10,100,300,1000,3000,10000,30000}; int ki=0;
  for(int k=0;k<V&&ki<8;k++){cum+=id[ord[k]]; if(k+1==ks[ki]){printf("top%d share of arcs %.3f\n",ks[ki],(double)cum/E);ki++;}}
  vector<int> rk(V); for(int k=0;k<V;k++)rk[ord[k]]=k;
  // CSR
  vector<int> h(V+1),to(E); vector<unsigned> w(E);
  for(int i=0;i<E;i++)h[eu[i]+1]++; for(int v=0;v<V;v++)h[v+1]+=h[v];
  {vector<int> c(h.begin(),h.end()-1); for(int i=0;i<E;i++){int k=c[eu[i]]++;to[k]=ev[i];w[k]=ew[i];}}
  // sample queries: full dijkstra with parent
  vector<long long> d(V); vector<int> par(V);
  int NS=300; long long maxd=0; vector<int> hopsh(64); int cover[8]={0}; int reach=0;
  std::map<int,int> minrk;
  for(int qi=0;qi<NS;qi++){
    int s=qs[qi*37%Q], t=qt[qi*37%Q];
    std::fill(d.begin(),d.end(),LLONG_MAX); Heap4<long long> H; d[s]=0; par[s]=-1; H.push(0,s);
    while(!H.empty()){auto it=H.pop(); if(it.d>d[it.v])continue; int u=it.v; for(int k=h[u];k<h[u+1];k++){long long nd=it.d+w[k]; if(nd<d[to[k]]){d[to[k]]=nd;par[to[k]]=u;H.push(nd,to[k]);}}}
    for(int v=0;v<V;v++) if(d[v]!=LLONG_MAX) maxd=std::max(maxd,d[v]);
    if(d[t]==LLONG_MAX||s==t) continue; reach++;
    int hop=0, mr=V; for(int x=t;x!=-1;x=par[x]){hop++; if(x!=s&&x!=t) mr=std::min(mr,rk[x]);}
    hopsh[std::min(hop,63)]++;
    for(int j=0;j<8;j++) if(mr<ks[j]) cover[j]++;
  }
  printf("maxdist seen %lld reach %d\n",maxd,reach);
  printf("hops:"); for(int i=0;i<64;i++) if(hopsh[i])printf(" %d:%d",i,hopsh[i]); printf("\n");
  for(int j=0;j<8;j++)printf("path interior touches top%d: %.3f\n",ks[j],(double)cover[j]/reach);
}
