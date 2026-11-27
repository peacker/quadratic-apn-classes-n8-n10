// Read LUTs of (8,8)-functions (256 integers per line, any separators).
// Output per line: idx \t status \t repr(ortho_invariant)   status: Q=quadratic APN, N=not quadratic, X=quadratic non-APN, B=bad line
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#define N 256
static int par(unsigned x){return __builtin_parity(x);}
static void pr_tuple(char *buf, long *h, int maxv){
  int first=1, cnt=0; char *p=buf; p+=sprintf(p,"(");
  for(int v=0;v<=maxv;v++) if(h[v]){ if(!first)p+=sprintf(p,", "); p+=sprintf(p,"(%d, %ld)",v,h[v]); first=0; cnt++; }
  if(cnt==1) p+=sprintf(p,",");
  sprintf(p,")");
}
int main(int argc,char**argv){
  long stride=1, shard=0; if(argc>2){shard=atol(argv[1]); stride=atol(argv[2]);}
  static char line[1<<16]; int F[N]; long idx=-1;
  while(fgets(line,sizeof line,stdin)){
    idx++; if(idx%stride!=shard) continue;
    int k=0; char*s=line;
    while(*s && k<=N){ while(*s && !isdigit((unsigned char)*s)) s++; if(!*s)break; F[k++]=strtol(s,&s,10);} 
    if(k!=N){ printf("%ld\tB\t-\n",idx); continue; }
    int ok=1;
    for(int i=0;i<8&&ok;i++)for(int j=i+1;j<8&&ok;j++){int a=1<<i,b=1<<j,c=F[0]^F[a]^F[b]^F[a^b];
      for(int x=0;x<N;x++) if((F[x]^F[x^a]^F[x^b]^F[x^a^b])!=c){ok=0;break;}}
    if(!ok){ printf("%ld\tN\t-\n",idx); continue; }
    int pi[N]; pi[0]=0;
    for(int a=1;a<N&&ok;a++){
      int img[8]; for(int i=0;i<8;i++){int x=1<<i; img[i]=F[x]^F[x^a]^F[a]^F[0];}
      int found=-1,nf=0;
      for(int c=1;c<N;c++){int z=1; for(int i=0;i<8;i++) if(par(c&img[i])){z=0;break;} if(z){found=c;nf++;}}
      if(nf!=1){ok=0;break;} pi[a]=found;
    }
    if(!ok){ printf("%ld\tX\t-\n",idx); continue; }
    long dh[N+1]; memset(dh,0,sizeof dh);
    int cnt[N];
    for(int a=1;a<N;a++){ memset(cnt,0,sizeof cnt); for(int x=0;x<N;x++) cnt[pi[x^a]^pi[x]]++; for(int b=0;b<N;b++) dh[cnt[b]]++; }
    long wh[N+1]; memset(wh,0,sizeof wh);
    int W[N];
    for(int c=1;c<N;c++){
      for(int x=0;x<N;x++) W[x]=par(c&pi[x])?-1:1;
      for(int h=1;h<N;h<<=1) for(int i=0;i<N;i+=2*h) for(int j=i;j<i+h;j++){int u=W[j],v=W[j+h];W[j]=u+v;W[j+h]=u-v;}
      for(int u=0;u<N;u++) wh[abs(W[u])]++;
    }
    char b1[8192],b2[8192]; pr_tuple(b1,dh,N); pr_tuple(b2,wh,N);
    printf("%ld\tQ\t(%s, %s)\n",idx,b1,b2);
  }
  return 0;
}
