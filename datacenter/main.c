#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "init.h"
#include "boundary.h"
#include "uvp.h"

double max(double** u,int imax,int jmax){
	int i,j;
	double max=0;
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
			if (u[j][i]>max)
				max=u[j][i];
		}
	}
    /*printf("max u:%f\n",max);*/
	return max;
}

void comp_delt(double* delt,int imax,int jmax,double delx,double dely,double **u,double **v,double Re,double tau){
    double first,second,third,min;
	double delta = 1/(delx*delx)+1/(dely*dely);
	first = Re/2/delta;
	min=first;
	second = delx/abs(max(u,imax,jmax));
	third= dely/abs(max(v,imax,jmax));
	if(min>second){
		min=second;
		if(min>third)
			min=third;
	}
	else{
		if(min>third)
			min=third;
	}
	*delt=tau*min;
	return;
}

int main(int argc,char* argv[]){
	int opt=0,n=0;
	double t=0;
	double xlength,ylength,width,height,uin,vin,inlet,outlet;
	double tend,tau,itermax,eps,omg,gamma;
	double Re,GX,GY,UI,VI,PI;
	int imax,jmax;
	int wW,wE,wS,wN;
	double x,y,x1,y1,x2,y2,u1,u2,u3,u4,v1,v2,v3,v4;
	double xB,xC,xD,xG,yI;
	double **u;
	double **v;
	double **p;
	double **f;
	double **g;
	double **rhs;
	double **uu;
	double **vv;
	double *xx;
	double *yy;
	double delx,dely,delt;
	int i,j,ii,jj;
	int iB,iC,iD,iG,jI;
    FILE *input;
    FILE *outputu;
    FILE *outputv;
    FILE *outputu1;
    FILE *outputv1;
    char *outputfilenameu;
    char *outputfilenamev;
    char *outputfilenameu1;
    char *outputfilenamev1;
    char *inputname=NULL;
    clock_t t1,t2;
    double  total_t;
	/*Read command line arguments*/
	do {
		opt = getopt(argc,argv,"f:"); /*why type -f?*/
		switch(opt){
			case 'f':
                inputname=optarg;
			break;

			default:
            break;
		}
	}while (opt!=-1);

	/* read inputs*/
	if (inputname==NULL){
		printf("Inputname is incorrect.\n");
		return -1;
	}
    input = fopen(inputname,"r");
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
	fscanf(input,"%lf %lf %lf %lf\n",&width,&height,&inlet,&outlet);
	fscanf(input,"%lf %lf\n",&uin,&vin);

	delx=xlength/imax;
	dely=ylength/jmax;
    printf("xlength:%f ylength:%f server width:%f server height:%f\n",xlength,ylength,width,height);
    printf("inlet diameter:%f outlet diameter:%f\n",inlet,outlet);
    printf("jmax:%d imax:%d Re:%f UI:%f VI:%f PI:%f GX:%f GY:%f tend:%f tau:%f itermax:%f eps:%f omg:%f gamma:%f\n",jmax,imax,Re,UI,VI,PI,GX,GY,tend,tau,itermax,eps,omg,gamma);
    printf("inlet u:%f inlet v:%f\n",uin,vin);

	/* assign initial values to u,v,p,f,g,rhs,uu,vv*/
	u=RMATRIX(0,jmax+1,0,imax+1);
	v=RMATRIX(0,jmax+1,0,imax+1);
	p=RMATRIX(0,jmax+1,0,imax+1);
	f=RMATRIX(0,jmax+1,0,imax+1);
	g=RMATRIX(0,jmax+1,0,imax+1);
	rhs=RMATRIX(0,jmax+1,0,imax+1);
	uu=RMATRIX(0,jmax,0,imax);
	vv=RMATRIX(0,jmax,0,imax);
	/* allocate memory to xx, yy*/
	xx=(double *)malloc((imax+1)*sizeof(double));
	yy=(double *)malloc((jmax+1)*sizeof(double));
	/* create geometric model of data center*/
	xB=width/2; xC=xB+0.3; xD=xC+inlet;xG=outlet/2;yI=height;
	iB=round(xB/delx);iC=round(xC/delx);iD=round(xD/delx);iG=round(xG/delx);jI=round(yI/dely);
	init_uvp(u,v,p,imax,jmax,UI,VI,PI);
    t1=clock();
	while(t<tend){
        if(n==0){
		/*comp_delt(&delt,imax,jmax,delx,dely,u,v,Re,tau);*/
        delt=0.02;
		setbound(u,v,imax,jmax,wW, wE,wN,wS,uin,vin,iB,iC,iD,iG,jI);
		comp_fg(u,v,f,g,imax,jmax,delt,delx,dely,GX,GY,gamma,Re,iB,iC,iD,iG,jI);
		comp_rhs(f, g,rhs,imax,jmax,delt,delx,dely);
		poisson(p,rhs,imax,jmax,delx,dely,eps,itermax,omg,iB,iC,iD,iG,jI);
		adap_uv(u,v,f,g,p,imax,jmax,delt,delx,dely,iB,jI);
		t=t+delt;
		n++;
		printf("The current t:%f\n",t);
        }else{
            comp_delt(&delt,imax,jmax,delx,dely,u,v,Re,tau);
            setbound(u,v,imax,jmax,wW,wE,wN,wS,uin,vin,iB,iC,iD,iG,jI);
            comp_fg(u,v,f,g,imax,jmax,delt,delx,dely,GX,GY,gamma,Re,iB,iC,iD,iG,jI);
            comp_rhs(f,g,rhs,imax,jmax,delt,delx,dely);
            poisson(p,rhs,imax,jmax,delx,dely,eps,itermax,omg,iB,iC,iD,iG,jI);
            adap_uv(u,v,f,g,p,imax,jmax,delt,delx,dely,iB,jI);
            t=t+delt;
            n++;
            /*printf("The current delt:%f\n",delt);*/
            printf("The current t:%f\n",t);
        }
    }
    t2=clock();
    total_t=(double)(t2-t1)/CLOCKS_PER_SEC;
    printf("Time elapsed:%f\n",total_t);

	outputfilenameu="outputu.txt";
    outputfilenamev="outputv.txt";
	outputu = fopen(outputfilenameu,"w+");
	outputv = fopen(outputfilenamev,"w+");

	for(j=0;j<jmax+2;j++){
		for(i=0;i<imax+2;i++){
		  fprintf(outputu,"%f ",u[j][i]);
		}
		fprintf(outputu,"\n");
	}
	for(j=0;j<jmax+2;j++){
		for(i=0;i<imax+2;i++){
            fprintf(outputv,"%f ",v[j][i]);
		}
		fprintf(outputv,"\n");
	}
	fclose(outputu);
	fclose(outputv);
	printf("u into file:%s\n",outputfilenameu);
	printf("v into file:%s\n",outputfilenamev);

	/** post for visualization*/
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
			u1=u[jj-1][ii-1];
			u2=u[jj-1][ii];
			u3=u[jj][ii-1];
			u4=u[jj][ii];
			uu[j][i]=1/(delx*dely)*((x2-x)*(y2-y)*u1+(x-x1)*(y2-y)*u2+(x2-x)*(y-y1)*u3+(x-x1)*(y-y1)*u4);
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
			v1=v[jj-1][ii-1];
			v2=v[jj-1][ii];
			v3=v[jj][ii-1];
			v4=v[jj][ii];
			vv[j][i]=1/(delx*dely)*((x2-x)*(y2-y)*v1+(x-x1)*(y2-y)*v2+(x2-x)*(y-y1)*v3+(x-x1)*(y-y1)*v4);
		}
	}

	outputfilenameu1="post_outputu.txt";
    outputfilenamev1="post_outputv.txt";
	outputu1 = fopen(outputfilenameu1,"w+");
	outputv1 = fopen(outputfilenamev1,"w+");
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
			fprintf(outputu1,"%f ",uu[j][i]);
		}
		fprintf(outputu1,"\n");
	}
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
            fprintf(outputv1,"%f ",vv[j][i]);
		}
		fprintf(outputv1,"\n");
	}
	fclose(outputu1);
	fclose(outputv1);
	printf("uu into file:%s\n",outputfilenameu1);
	printf("vv into file:%s\n",outputfilenamev1);

	FREE_RMATRIX(u,0,jmax+1,0,imax+1);
	FREE_RMATRIX(v,0,jmax+1,0,imax+1);
	FREE_RMATRIX(p,0,jmax+1,0,imax+1);
	FREE_RMATRIX(f,0,jmax+1,0,imax+1);
	FREE_RMATRIX(g,0,jmax+1,0,imax+1);
	FREE_RMATRIX(rhs,0,jmax+1,0,imax+1);
	FREE_RMATRIX(uu,0,jmax,0,imax);
	FREE_RMATRIX(vv,0,jmax,0,imax);
	free(xx);
	free(yy);
    return 0;
}
