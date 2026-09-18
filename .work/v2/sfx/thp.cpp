#include <cstdio>
#include <vector>
#include <cstdint>
#include <cstring>
int main(){
  std::vector<uint16_t> a(300000*64, 1);
  std::vector<uint32_t> b(300000*8, 1);
  char line[256]; FILE* f=fopen("/proc/self/smaps_rollup","r"); while(fgets(line,256,f)) if(strstr(line,"Anon")||strstr(line,"Rss")) fputs(line,stdout);
  printf("%p %p\n",(void*)a.data(),(void*)b.data());
}
