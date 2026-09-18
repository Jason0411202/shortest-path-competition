#include <cstdio>
#include <vector>
#include <chrono>
#include <random>
#include <numeric>
#include <algorithm>
int main(){
  for(size_t MB : {1,4,8,16,64,256}){
    size_t n=MB*1024*1024/64; std::vector<size_t> perm(n); std::iota(perm.begin(),perm.end(),0);
    std::shuffle(perm.begin(),perm.end(),std::mt19937_64(1));
    struct alignas(64) L{size_t nx; char pad[56];}; std::vector<L> a(n);
    for(size_t i=0;i<n;i++) a[perm[i]].nx=perm[(i+1)%n];
    size_t p=0; auto t0=std::chrono::steady_clock::now(); size_t it=20000000;
    for(size_t i=0;i<it;i++) p=a[p].nx;
    double s=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    printf("%zu MB: %.1f ns/access (%zu)\n",MB,s/it*1e9,p);
  }
}
