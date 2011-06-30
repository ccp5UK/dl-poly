#include <stdio.h>

/* dl_poly decryption program */
/* copyright daresbury laboratory */
/* author w.smith */

int ir;
int jr;
float ww;
float cd;
float cm;
float u[97];
char key[12];

void setrnd()
{
  int i,j,ii,jj,kk,ll,mm;
  float s,t;
  ir=96;
  jr=32;
  ii=1+(((unsigned int)key[0])*((unsigned int)key[1])*((unsigned int)key[2]))%177;
  jj=1+(((unsigned int)key[3])*((unsigned int)key[4])*((unsigned int)key[5]))%177;
  kk=1+(((unsigned int)key[6])*((unsigned int)key[7])*((unsigned int)key[8]))%177;
  ll=(((unsigned int)key[9])*((unsigned int)key[10])*((unsigned int)key[11]))%168;
  for(i=0;i<97;i++)
    {
      s=0.0;
      t=0.5;
      for(j=0;j<24;j++)
	{
	  mm=(((ii*jj)%179)*kk)%179;
	  ii=jj;
	  jj=kk;
	  kk=mm;
	  ll=(53*ll+1)%169;
	  if((ll*mm)%64>=32)s=s+t;
	  t=0.5*t;
	}
      u[i]=s;
    }
  ww=  362436.0/16777216.0;
  cd= 7654321.0/16777216.0;
  cm=16777213.0/16777216.0;
}

float randum()
{
  float uni;
  uni=u[ir]-u[jr];
  if(uni<0.0)uni+=1.0;
  u[ir]=uni;
  ir--;
  if(ir<0)ir=96;
  jr--;
  if(jr<0)jr=96;
  ww-=cd;
  if(ww<0.0)ww+=cm;
  uni-=ww;
  if(uni<0.0)uni+=1.0;
  return uni;
}

main(int  argc, char *argv[])
{
  char fname[100],ename[100];
  FILE *fpi,*fpo,*fopen();
  int kkk[1000];
  int i,j,k,c,n;
  if(argc == 1)
    {
      printf("\nEnter the password: ");
      gets(key);
      printf("\nEnter the file name for decryption: ");
      gets(fname);
    }
  else
    {
      strcpy(key,argv[1]);
      strcpy(fname,argv[2]);
    }
  n=1000;
  setrnd();
  for(i=0;i<n;i++)
    kkk[i]=(int)(256.0*randum());
  strcpy(ename,fname);
  ename[strlen(fname)-1]='z';
  printf("\nThe output file will be named: %s\n",ename);
  if((fpi=fopen(fname,"r"))==NULL)
    {
      printf("\nError - file %s not found",fname);
      exit(1);
    }
  if((fpo=fopen(ename,"w"))==NULL)
    {
      printf("\nError - file %s not opened",ename);
      exit(1);
    }
  j=0;
  while((c = getc(fpi)) != EOF)
    {
      if(j  == n) j=0;
      k=(unsigned int)c-kkk[j];
      if(k < 0)k+=256;
      c=(char)k;
      putc(c,fpo);
      j++;
    }
}

