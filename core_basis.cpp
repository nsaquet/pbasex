/*
    Filename: core_basis.c
    To be used with pBaseCore.py, as an imported library. Use gcc to compile.
    gcc -shared -Wall -o libCoreBasis.so -fPIC core_basis.c /Developer/extralibs/lib/libgsl.a /Developer/extralibs/lib/libgslcblas.a -I/Developer/extralibs/include/ -arch i386 -arch x86_64
*/
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

using namespace std;

#define NR 256 // Radial binning 256
#define NTH 256 // Angular binning (180 degrees) 256
#define NK 128 // Number of basis functions 128

/* Function prototypes */
extern "C" float lpoly(float,float *,int);
extern "C" double theta_f(double,double);
extern "C" double model_f(double,void*);
extern "C" double model_c(int,int,double ,double,gsl_integration_workspace *);
extern "C" void write_forward(int,int);

/* Evaluates the Forward Abel Integral for a given basis function (ik,l) */
/* Where ik is the radial function and l is the Pl Legendre polynomial   */

double model_f(double r,void *par)
{
  float lpoly(float,float [],int);
  float *pl;
  double *xz=(double*)par;
  double theta_f(double,double);
  double thta=theta_f(r,xz[1]);
  double R=r/sin(thta),f;
  int k=(int)xz[2]*2; //change for different pixel widths
  int l=(int)xz[3];

  pl=(float*)calloc(l+1,sizeof(float));

  f=exp(-(R-k)*(R-k)/2.0);
  f*=lpoly(thta,pl,l);
  free(pl);

  return f*r/sqrt(r*r-xz[0]*xz[0]);
}

double model_c(int k,int l,double x,double z,gsl_integration_workspace *w)
{
  double par[4]={x,z,k,l};
  gsl_function f;
  //size_t neval;
  double error,result;
  f.function=&model_f;
  f.params=&par;
  gsl_integration_qag(&f,fabs(x),300.0,0.00001,0.00005,100000,6,w,&result,&error);
  return result;
}

/* Return Legendre polynomial at X */
float lpoly(float X, float pl[], int npl)
{
  int j;
  float twox,f2,f1,d,x;

  x=(float)cos(X);
  pl[0]=1.0;
  pl[1]=x;
  if (npl >= 2) {
    twox=2.0*x;
    f2=x;
    d=1.0;
    for (j=2;j <= npl ;j++) {
      f1=d++;
      f2 += twox;
      pl[j]=(f2*pl[j-1]-f1*pl[j-2])/d;
    }
  }
  return pl[npl];
}

/* Return azimut angle from coordinates X,Y */
double theta_f(double x,double y)
{
  double thta;
  thta=atan(fabs(x)/fabs(y));
  if (x<0&&y>0) thta=2.0*M_PI-thta;
  if (x>0&&y<0) thta=M_PI-thta;
  if (x<0&&y<0) thta=M_PI+thta;
  if(x==0&&y>0) thta=0.0;
  if(x==0&&y<0) thta=M_PI;
  if(y==0&&x>0) thta=M_PI/2.0;
  if(y==0&&x<0) thta=3.0*M_PI/2.0;
  if(x==0&&y==0) thta=2.0*M_PI;
  return thta;
}

/***********************************************************************/
/* This function creates the matrix N,M of forward basis functions     */
/* where N=number of radial functions x number of legendre polynomials */
/* Then it decomposes the matrix using the SVD method and writes the   */
/* matrices U, V, and S to files                                       */
/***********************************************************************/


void write_forward(int l,int odd)
{
  int ir,it,ik,il,nl,basis_index,point_index,NL,N,M,i;
  double p,dx,dz,thta,rad;
  float cpu1,cpu0;
  gsl_matrix *Basis,*BasisT,*V,*X;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(100000); // Allocate space for integration
  gsl_vector *W,*S;
  FILE *fs,*fu,*fv;
  char fileV[50],fileU[50],fileS[50],file[50],str[5];

  if(odd) NL=l+1;
  else NL=l/2+1;
  N=NL*NK;
  M=NR*NTH;
  Basis=gsl_matrix_alloc(N,M); //Allocate Basis matrix
  if (Basis==NULL){
    printf("Couldn't allocate memory for basis functions\n");
    exit(0);
  }
  for(ik=0;ik<NK;ik++){
    for(il=0;il<NL;il++){
      basis_index=ik*NL+il;     
      cpu0=(float)clock()/(float)CLOCKS_PER_SEC;
      if(odd) nl=il;
      else nl=il*2;
      for(ir=0;ir<NR;ir++){
	    for(it=0;it<NTH;it++){
	        point_index=ir*NTH+it;
	        rad=(double)ir;
	        thta=it/(double)(NTH-1)*M_PI;
	        dx=rad*sin(thta);
	        dz=rad*cos(thta);
	        p=model_c(ik,nl,dx,dz,w);
	        gsl_matrix_set(Basis,basis_index,point_index,p);
	        }
        }
      cpu1=(float)clock()/(float)CLOCKS_PER_SEC;
      printf("Basis set %d,%d (index %d of %d) done in %f seconds\n",ik,il,basis_index,NK*NL,cpu1-cpu0);
    }
  }
  gsl_integration_workspace_free(w);

  //printf("Decomposing matrix, wait please ...\n");

  /***********************/
  /* Decomposition stuff */
  /***********************/

  /* Transpose it */

  BasisT=gsl_matrix_alloc(M,N);
  gsl_matrix_transpose_memcpy(BasisT,Basis);
  gsl_matrix_free(Basis);
  
  /*Allocate working space for decomposition */

  S=gsl_vector_alloc(N);
  V=gsl_matrix_alloc(N,N);
  W=gsl_vector_alloc(N);
  X=gsl_matrix_alloc(N,N);

  /* Decompose */
  gsl_linalg_SV_decomp_mod(BasisT,X,V,S,W);
  gsl_matrix_free(X);
  gsl_vector_free(W);

  /* Describe matrices' names */

  memset(fileV,NULL,50);
  memset(fileU,NULL,50);
  memset(fileS,NULL,50);
  memset(file,NULL,50);

  for(i=0;i<=l;i++)
  {
    if((i%2))
    {
      sprintf(str,"P%d",i);
      if(odd) strcat(file,str);
    }
    else
    {
      sprintf(str,"P%d",i);
      strcat(file,str);
    }
  }

  strcat(file,".dat"); // Add extension

  sprintf(fileU,"U");
  sprintf(fileV,"V");
  sprintf(fileS,"S");
  
  strcat(fileU,file);
  strcat(fileV,file);
  strcat(fileS,file);

  /* Write everything to file */

  fu=fopen(fileU,"w");
  gsl_matrix_fwrite(fu,BasisT);
  gsl_matrix_free(BasisT);
  fclose(fu);

  fv=fopen(fileV,"w");
  gsl_matrix_fwrite(fv,V);
  gsl_matrix_free(V);
  fclose(fv);

  fs=fopen(fileS,"w");
  gsl_vector_fwrite(fs,S);
  gsl_vector_free(S);
  fclose(fs);
  
}