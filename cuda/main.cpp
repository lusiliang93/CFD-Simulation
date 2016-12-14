#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "cuda_util.h"

int imax, jmax;

#define get_index(i,j) ((jmax+2)*i+j)

int read_parameter(char *inputfile){
    return 0;
}

double *RMATRIX(int nrl,int nrh, int ncl,int nch){
/* reserve memory for matrix of size [nrl,nrh]x[ncl,nch]*/
    double* m;
    m = (double*)malloc((unsigned)(nrh-nrl+1)*(nch-ncl+1)*sizeof(double));
    if(m==NULL){
        printf("no more memory \n");
        exit(0);
    }
    return m;
}

void FREE_RMATRIX(double *m){
    free(m);
    return;
}

extern double *cudaDevice_u;
extern double *cudaDevice_v;
extern double *cudaDevice_p;
extern double *cudaDevice_f;
extern double *cudaDevice_g;
extern double *cudaDevice_rhs;
extern double *cudaDevice_u2;
extern double *cudaDevice_v2;
extern double *cudaDevice_p2;
extern double *cudaDevice_f2;
extern double *cudaDevice_g2;
extern double *cudaDevice_rhs2;

int main(int argc,char* argv[]){
    int n=0;
    double t=0;
    double xlength,ylength;
    double tend,tau,itermax,eps,omg,gamma;
    double Re,GX,GY,UI,VI,PI;
    int wW,wE,wS,wN;
    double delx,dely,delt;
    FILE *input;
    clock_t t1,t2;
    double  total_t;
    input = fopen("input.txt","r");
    if(!input){
        printf("Unable to open the file.\n");
        return -1;
    }

    fscanf(input,"%lf %lf\n",&xlength,&ylength);
    fscanf(input,"%d %d\n",&imax,&jmax);
    fscanf(input,"%lf %lf %lf %lf %lf %lf\n",&Re,&UI,&VI,&PI,&GX,&GY);
    /* tau must be postive*/
    fscanf(input,"%lf %lf %lf %lf %lf %lf \n",&tend,&tau,&itermax,&eps,&omg,&gamma);
    fscanf(input,"%d %d %d %d\n",&wW,&wE,&wN,&wS);
    delx=xlength/imax;
    dely=ylength/jmax;
    printf("xlengh:%f ylength:%f jmax:%d imax:%d Re:%f UI:%f VI:%f PI:%f GX:%f GY:%f tend:%f tau:%f itermax:%f eps:%f omg:%f gamma:%f\n",xlength,ylength,jmax,imax,Re,UI,VI,PI,GX,GY,tend,tau,itermax,eps,omg,gamma);

    init_uvp(imax,jmax,UI,VI,PI);
    printf("success init_uvp\n");
    t1=clock();
    cuda_init(imax, jmax);
    printf("success cuda_init\n");
    while(t<tend){
    	copy_matrix(imax, jmax);
    	printf("success copy_matrix\n");
    	printf("pointers: %p %p %p %p\n", cudaDevice_u, cudaDevice_u2, cudaDevice_v, cudaDevice_v2);
    	printf("check max_vector: %lf\n", max_vector(cudaDevice_u, (imax+2)*(jmax+2)));
    	printf("check max_vector: %lf\n", max_vector(cudaDevice_u, (imax+2)*(jmax+2)));
        if(n==0){
            delt=0.02;
            // setbound(imax,jmax,wW,wE,wN,wS);
            // comp_fg(imax,jmax,delt,delx,dely,GX,GY,gamma,Re);
            // comp_rhs(imax,jmax,delt,delx,dely);
            // poisson(imax,jmax,delx,dely,eps,itermax,omg);
            // adap_uv(imax,jmax,delt,delx,dely);
            t=t+delt;
            n++;
            printf("The current t:%f\n",t);
        }else{
        	printf("defore comp_delt\n");
            delt = comp_delt(imax,jmax,delx,dely,Re,tau);
            printf("after comp_delt\n");
            // setbound(imax,jmax,wW,wE,wN,wS);
            // comp_fg(imax,jmax,delt,delx,dely,GX,GY,gamma,Re);
            // comp_rhs(imax,jmax,delt,delx,dely);
            // poisson(imax,jmax,delx,dely,eps,itermax,omg);
            // adap_uv(imax,jmax,delt,delx,dely);
            t=t+delt;
            n++;
            printf("The current t:%f\n",t);
        }

    }
    t2=clock();
    total_t=(double)(t2-t1)/CLOCKS_PER_SEC;
    printf("Time elapsed:%f\n",total_t);
    return 0;
}
