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

void adap_uv(double *u,double *v,double *f,double *g,double *p,double delt,double delx,double dely){
    int i,j;
    for(j=1;j<jmax+1;j++){
        for(i=1;i<imax;i++){
            u[get_index(j,i)] = f[get_index(j,i)] - delt/delx*(p[get_index(j,i+1)]-p[get_index(j,i)]);
        }
    }
    for(j=1;j<jmax;j++){
        for(i=1;i<imax+1;i++){
            v[get_index(j,i)] = g[get_index(j,i)] - delt/dely*(p[get_index(j+1,i)]-p[get_index(j,i)]);
        }
    }
    return;
}


double max(double* u){
    int i,j;
    double max=0;
    for(j=0;j<jmax+1;j++){
        for(i=0;i<imax+1;i++){
            if(u[get_index(j,i)]>max){
                max = u[get_index(j,i)];
            }
        }
    }
    return max;
}

int main(int argc,char* argv[]){
    int n=0;
    double t=0;
    double xlength,ylength;
    double tend,tau,itermax,eps,omg,gamma;
    double Re,GX,GY,UI,VI,PI;
    int wW,wE,wS,wN;
    double x,y,x1,y1,x2,y2,u1,u2,u3,u4,v1,v2,v3,v4;
    double delx,dely,delt;
    int i,j,ii,jj;
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
            delt = comp_delt(delx,dely,Re,tau);
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
