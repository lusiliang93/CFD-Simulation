#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include "init.h"
#include "boundary.h"
#include "uvp.h"
float max(float** u,int imax,int jmax){
	int i,j;
	float max=0;
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
			if (u[j][i]>max)
				max=u[j][i];
		}
	}
	return max;
}

void comp_delt(float* delt,int imax,int jmax,float delx,float dely,float **u,float **v,float Re,float tau){
    float first,second,third,min;
	float delta = 1/(delx*delx)+1/(dely*dely);
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
	float t=0;
	float xlength,ylength;
	float tend,tau,itermax,eps,omg,gamma;
	float Re,GX,GY,UI,VI,PI;
	int imax,jmax;
	int wW,wE,wS,wN;
	float **u;
	float **v;
	float **p;
	float **f;
	float **g;
	float **rhs;
	float delx,dely,delt;
	int i,j;
    FILE *input;
    FILE *outputu;
    FILE *outputv;
    char *outputfilenameu;
    char *outputfilenamev;
    char *inputname=NULL;
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

	fscanf(input,"%f %f\n",&xlength,&ylength);
	fscanf(input,"%d %d\n",&imax,&jmax);
	fscanf(input,"%f %f %f %f %f %f\n",&Re,&UI,&VI,&PI,&GX,&GY);
	/* tau must be postive*/
	fscanf(input,"%f %f %f %f %f %f \n",&tend,&tau,&itermax,&eps,&omg,&gamma);
	fscanf(input,"%d %d %d %d\n",&wW,&wE,&wN,&wS);
	delx=xlength/imax;
	dely=ylength/jmax;

	/* assign initial values to u,v,p,f,g,rhs*/
	u=RMATRIX(0,imax+1,0,jmax+1);
	v=RMATRIX(0,imax+1,0,jmax+1);
	p=RMATRIX(0,imax+1,0,jmax+1);
	f=RMATRIX(0,imax+1,0,jmax+1);
	g=RMATRIX(0,imax+1,0,jmax+1);
	rhs=RMATRIX(0,imax+1,0,jmax+1);
	init_uvp(u,v,p,imax,jmax,UI,VI,PI);
	while(t<tend){
		comp_delt(&delt,imax,jmax,delx,dely,u,v,Re,tau); /*Correct?*/
		setbound(u,v,imax,jmax,wW, wE,wN,wS);
		comp_fg(u,v,f,g, imax,jmax,delt,delx,dely,GX,GY,gamma,Re);
		comp_rhs(f, g,rhs,imax,jmax,delt,delx,dely);
		poisson(p,rhs,imax,jmax,delx,dely,eps,itermax,omg);
		adap_uv(u,v,f,g,p,imax,jmax,delt,delx,dely);
		t=t+delt;
		n++;
		printf("The current t:%f\n",t);
	}

	outputfilenameu="outputu.txt";
    outputfilenamev="outputv.txt";
	outputu = fopen(outputfilenameu,"w+");
	outputv = fopen(outputfilenamev,"w+");
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
			fprintf(outputu,"%f ",u[j][i]);
		}
		fprintf(outputu,"\n");
	}
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
			fprintf(outputv,"%f ",v[j][i]);
		}
		fprintf(outputv,"\n");
	}
	fclose(outputu);
	fclose(outputv);
	printf("u into file:%s\n",outputfilenameu);
	printf("v into file:%s\n",outputfilenamev);

	FREE_RMATRIX(u,0,imax+1,0,jmax+1);
	FREE_RMATRIX(v,0,imax+1,0,jmax+1);
	FREE_RMATRIX(p,0,imax+1,0,jmax+1);
	FREE_RMATRIX(f,0,imax+1,0,jmax+1);
	FREE_RMATRIX(g,0,imax+1,0,jmax+1);
	FREE_RMATRIX(rhs,0,imax+1,0,jmax+1);
	/*post for visualization*/
    return 0;
}
