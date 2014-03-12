/*
    Filename: core_basis.c
    So far this is a standalone application for generating basis.
    g++ -O4 -Wall -o makebasis.exe core_basis.cpp /Developer/usr/lib/libgsl.a /Developer/usr/lib/libgslcblas.a -I/Developer/usr/include/

*/
#include <time.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

using namespace std;

//Gustavo is using a 440 by 440 matrix
#define NR 440 // Radial binning 256
#define NTH 440 // Angular binning (180 degrees) 256
#define NK 220 // Number of basis functions 128
#define BSPACE	2.0
#define BWIDTH	1.0

/* Function prototypes */
float lpoly(float,float *,int);
double theta_f(double,double);
double model_f(double,void*);
double model_c(int,int,double ,double,gsl_integration_workspace *);
void write_forward(int,int);

/* Evaluates the Forward Abel Integral for a given basis function (ik,l) */
/* Where ik is the radial function and l is the Pl Legendre polynomial   */

double model_f(double r,void *par)
{
  float lpoly(float,float [],int);
  float *pl;
  double *xz=(double*)par;
  double R,f,thta;
  thta=theta_f(r,xz[1]);
  R=sqrt(r*r + xz[1]*xz[1]);
  double k0,width; 
  int l;
  k0=xz[2]; //change for different pixel widths
  l=(int)xz[3];
  width=xz[4];
  
  pl=(float*)calloc(l+1,sizeof(float));

  f=exp(-(R-k0)*(R-k0)/width);
  f*=lpoly(thta,pl,l);
  free(pl);
  return f*r/sqrt(r*r-xz[0]*xz[0]);
}

double model_c(double k,int l,double x,double z,double width, gsl_integration_workspace *w)
{
  double par[5]={x,z,k,l,width};
  gsl_function f;
  //size_t neval;
  int status;
  double epsrel;
  epsrel=0.00005;
  double error,result;
  f.function=&model_f;
  f.params=&par;
  while(1)
  {
    status=gsl_integration_qag(&f,fabs(x),NR+50.,0.0,epsrel,100000,4,w,&result,&error);
    switch(status)
    {
        case 0:
            return result;
            break;
        case GSL_EFAILED:
        case GSL_EROUND:
            epsrel*=1.2;
            break;
        default:
            return result;
            break;
    }
  }
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
  int ir,it,ik,il,nl,basis_index,point_index,NL,N,M,i,nth;
  double p,dx,dz,thta,rad,space,width,sigma,k0;
  space=BSPACE;
  sigma=BWIDTH;
  float cpu1,cpu0;
  gsl_matrix *Basis,*V,*X;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(100000); // Allocate space for integration
  gsl_vector *W,*S;
  FILE *fs,*fu,*fv;
  char fileV[50],fileU[50],fileS[50],file[50],str[5];

  if(odd) NL=l+1;
  else NL=l/2+1;
  N=NL*NK;
  M=NR*NTH;
  Basis=gsl_matrix_alloc(M,N); //Allocate Basis matrix
  if (Basis==NULL){
    printf("Couldn't allocate memory for basis functions\n");
    exit(0);
  }
  
  gsl_set_error_handler_off(); // I will take care of integration errors myself
  
  width=2*sigma*sigma;
  for(ik=0;ik<NK;ik++){
    for(il=0;il<NL;il++){
      basis_index=ik*NL+il;     
      cpu0=(float)clock()/(float)CLOCKS_PER_SEC;
      if(odd) nl=il;
      else nl=il*2;
      point_index=0;
      for(ir=0;ir<NR;ir++){
	    nth=2*ir+1;
	    for(it=0;it<nth;it++){
	        rad=(double)ir;
	        thta=(double)it/(double)(nth)*M_PI;
	        dx=rad*sin(thta);
	        dz=rad*cos(thta);
	        k0=ik*space;
	        p=model_c(k0,nl,dx,dz,width,w);
	        gsl_matrix_set(Basis,point_index,basis_index,p);
	        point_index++;
	        }
        }
      cpu1=(float)clock()/(float)CLOCKS_PER_SEC;
      printf("Basis set %d,%d (index %d of %d) done in %f seconds\n",ik,il,basis_index,NK*NL,cpu1-cpu0);
    }
  }
  gsl_integration_workspace_free(w);

  //printf("Decomposing matrix, wait please ...\n");
  
  //SAVING the basis BEFORE decomposition
  memset(fileU,0,50);
  memset(file,0,50);
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
  strcat(fileU,file);
  fu=fopen(fileU,"w");
  gsl_matrix_fwrite(fu,Basis);
  gsl_matrix_free(Basis);
  fclose(fu);

  /***********************/
  /* Decomposition stuff */
  /***********************/
  
  /*Allocate working space for decomposition */
  /*
  S=gsl_vector_alloc(N);
  V=gsl_matrix_alloc(N,N);
  W=gsl_vector_alloc(N);
  X=gsl_matrix_alloc(N,N);
  */
  /* Decompose */
  /*
  gsl_linalg_SV_decomp_mod(Basis,X,V,S,W);
  gsl_matrix_free(X);
  gsl_vector_free(W);
  */

  /* Describe matrices' names */
  /*
  memset(fileV,0,50);
  memset(fileU,0,50);
  memset(fileS,0,50);
  memset(file,0,50);

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
  */
  /* Write everything to file */
  /*
  fu=fopen(fileU,"w");
  gsl_matrix_fwrite(fu,Basis);
  gsl_matrix_free(Basis);
  fclose(fu);

  fv=fopen(fileV,"w");
  gsl_matrix_fwrite(fv,V);
  gsl_matrix_free(V);
  fclose(fv);

  fs=fopen(fileS,"w");
  gsl_vector_fwrite(fs,S);
  gsl_vector_free(S);
  fclose(fs);
  */
}

int main ()
{
	write_forward(6,0);
	write_forward(4,1);
	return 0;
}