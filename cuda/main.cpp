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

void gen_output(){
	double x,y,x1,y1,x2,y2,u1,u2,u3,u4,v1,v2,v3,v4;
	int i,j;
	double* u=RMATRIX(0,imax+1,0,jmax+1);
	double* v=RMATRIX(0,imax+1,0,jmax+1);
	cudaMemcpy(u,cudaDevice_u,(imax+2)*(jmax+2)*sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(v,cudaDevice_v,(imax+2)*(jmax+2)*sizeof(double),cudaMemcpyDeviceToHost);
	double* uu=RMATRIX(0,imax,0,jmax);
	double* vv=RMATRIX(0,imax,0,jmax);
	double* xx=(double *)malloc((imax+1)*sizeof(double));
	double* yy=(double *)malloc((jmax+1)*sizeof(double));
	outputu = fopen("outputu.txt","w+");
	outputv = fopen("outputv.txt","w+");
	for(j=0;j<jmax+2;j++){
	    for(i=0;i<imax+2;i++){
	        fprintf(outputu,"%f ", u[get_index(j,i)]);
	    }
	    fprintf(outputu,"\n");
	}
	for(j=0;j<jmax+2;j++){
	    for(i=0;i<imax+2;i++){
	        fprintf(outputv,"%f ", v[get_index(j,i)]);
	    }
	    fprintf(outputv,"\n");
	}
	fclose(outputu);
	fclose(outputv);
	printf("u into file\n");
	printf("v into file\n");

	/* post for visualization*/
	for(i=0;i<imax+1;i++){
	    xx[i]=delx*i;
	}
	for(j=0;j<jmax+1;j++){
	    yy[j]=dely*j;
	}
	for(j=0;j<jmax+1;j++){/*Is that right?*/
	    for(i=0;i<imax+1;i++){/*Is that right?*/
	        x=xx[i];/*Is that right?*/
	        y=yy[j];/*Is that right?*/
	        ii=floor(x/delx)+1;
	        jj=floor((y+dely/2)/dely)+1;
	        x1=(ii-1)*delx;
	        y1=((jj-1)-0.5)*dely;
	        x2=ii*delx;
	        y2=(jj-1/2)*dely;
	        u1 = u[get_index(jj-1,ii-1)];
	        u2 = u[get_index(jj-1,ii)];
	        u3 = u[get_index(jj,ii-1)];
	        u4 = u[get_index(jj,ii)];
	        uu[j*(jmax+1)+i] = 1/(delx*dely)*((x2-x)*(y2-y)*u1+(x-x1)*(y2-y)*u2+(x2-x)*(y-y1)*u3+(x-x1)*(y-y1)*u4);
	    }
	}
	/* vv*/
	for(j=0;j<jmax+1;j++){/*Is that right?*/
	    for(i=0;i<imax+1;i++){
	        x=xx[i];
	        y=yy[j];
	        jj=floor(y/dely)+1;
	        ii=floor((x+delx/2)/delx)+1;
	        y1=(jj-1)*dely;
	        x1=((ii-1)-0.5)*delx;
	        y2=jj*dely;
	        x2=(ii-0.5)*delx;
	        v1 = v[get_index(jj-1,ii-1)];
	        v2 = v[get_index(jj-1,ii)];
	        v3 = v[get_index(jj,ii-1)];
	        v4 = v[get_index(jj,ii)];
	        vv[j*(jmax+1)+i] = 1/(delx*dely)*((x2-x)*(y2-y)*v1+(x-x1)*(y2-y)*v2+(x2-x)*(y-y1)*v3+(x-x1)*(y-y1)*v4);
	    }
	}

	FILE* outputu1 = fopen("post_outputu.txt","w+");
	FILE* outputv1 = fopen("post_outputv.txt","w+");
	for(j=0;j<jmax+1;j++){
	    for(i=0;i<imax+1;i++){
	        fprintf(outputu1,"%f ", uu[j*(jmax+1)+i]);
	    }
	    fprintf(outputu1,"\n");
	}
	for(j=0;j<jmax+1;j++){
	    for(i=0;i<imax+1;i++){
	        fprintf(outputv1,"%f ", vv[j*(jmax+1)+i]);
	    }
	    fprintf(outputv1,"\n");
	}
	fclose(outputu1);
	fclose(outputv1);
	printf("uu into file\n");
	printf("vv into file\n");

	FREE_RMATRIX(u);
	FREE_RMATRIX(v);
	FREE_RMATRIX(uu);
	FREE_RMATRIX(vv);
	free(xx);
	free(yy);
}

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
    t1=clock();
    cuda_init(imax, jmax);
    while(t<tend){
    	copy_matrix(imax, jmax);
        if(n==0){
            delt=0.02;
            setbound(imax,jmax,wW,wE,wN,wS);
            comp_fg(imax,jmax,delt,delx,dely,GX,GY,gamma,Re);
            comp_rhs(imax,jmax,delt,delx,dely);
            poisson(imax,jmax,delx,dely,eps,itermax,omg);
            adap_uv(imax,jmax,delt,delx,dely);
            t=t+delt;
            n++;
            printf("The current t:%f\n",t);
        }else{
            delt = comp_delt(imax,jmax,delx,dely,Re,tau);
            setbound(imax,jmax,wW,wE,wN,wS);
            comp_fg(imax,jmax,delt,delx,dely,GX,GY,gamma,Re);
            comp_rhs(imax,jmax,delt,delx,dely);
            poisson(imax,jmax,delx,dely,eps,itermax,omg);
            adap_uv(imax,jmax,delt,delx,dely);
            t=t+delt;
            n++;
            printf("The current t:%f\n",t);
        }
    }
    t2=clock();
    total_t=(double)(t2-t1)/CLOCKS_PER_SEC;
    printf("Time elapsed:%f\n",total_t);
    gen_output();
    return 0;
}
