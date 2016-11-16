#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
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
	double xlength,ylength;
	double tend,tau,itermax,eps,omg,gamma;
	double Re,GX,GY,UI,VI,PI;
	int imax,jmax;
	int wW,wE,wS,wN;
	double **u;
	double **v;
	double **p;
	double **f;
	double **g;
	double **rhs;
	double delx,dely,delt;
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


	fscanf(input,"%lf %lf\n",&xlength,&ylength);
	fscanf(input,"%d %d\n",&imax,&jmax);
	fscanf(input,"%lf %lf %lf %lf %lf %lf\n",&Re,&UI,&VI,&PI,&GX,&GY);
	/* tau must be postive*/
	fscanf(input,"%lf %lf %lf %lf %lf %lf \n",&tend,&tau,&itermax,&eps,&omg,&gamma);
	fscanf(input,"%d %d %d %d\n",&wW,&wE,&wN,&wS);
	delx=xlength/imax;
	dely=ylength/jmax;
    printf("xlengh:%f ylength:%f jmax:%d imax:%d Re:%f UI:%f VI:%f PI:%f GX:%f GY:%f tend:%f tau:%f itermax:%f eps:%f omg:%f gamma:%f\n",xlength,ylength,jmax,imax,Re,UI,VI,PI,GX,GY,tend,tau,itermax,eps,omg,gamma);

	/* assign initial values to u,v,p,f,g,rhs*/
	u=RMATRIX(0,imax+2,0,jmax+2);
	v=RMATRIX(0,imax+2,0,jmax+2);
	p=RMATRIX(0,imax+2,0,jmax+2);
	f=RMATRIX(0,imax+2,0,jmax+2);
	g=RMATRIX(0,imax+2,0,jmax+2);
	rhs=RMATRIX(0,imax+2,0,jmax+2);

	init_uvp(u,v,p,imax,jmax,UI,VI,PI);

	while(t<tend){
        if(n==0){
		/*comp_delt(&delt,imax,jmax,delx,dely,u,v,Re,tau); */
        delt=0.02;
		setbound(u,v,imax,jmax,wW, wE,wN,wS);
		comp_fg(u,v,f,g, imax,jmax,delt,delx,dely,GX,GY,gamma,Re);
		comp_rhs(f, g,rhs,imax,jmax,delt,delx,dely);
		poisson(p,rhs,imax,jmax,delx,dely,eps,itermax,omg);
		adap_uv(u,v,f,g,p,imax,jmax,delt,delx,dely);
		t=t+delt;
		n++;
		printf("The current t:%f\n",t);
        }else{
            comp_delt(&delt,imax,jmax,delx,dely,u,v,Re,tau);
            setbound(u,v,imax,jmax,wW,wE,wN,wS);
            comp_fg(u,v,f,g,imax,jmax,delt,delx,dely,GX,GY,gamma,Re);
            comp_rhs(f,g,rhs,imax,jmax,delt,delx,dely);
            poisson(p,rhs,imax,jmax,delx,dely,eps,itermax,omg);
            adap_uv(u,v,f,g,p,imax,jmax,delt,delx,dely);
            t=t+delt;
            n++;
            printf("test u:%f\n",u[64][64]);
            printf("The current delt:%f\n",delt);
            printf("The current t:%f\n",t);
        }

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

	FREE_RMATRIX(u,0,imax+2,0,jmax+2);
	FREE_RMATRIX(v,0,imax+2,0,jmax+2);
	FREE_RMATRIX(p,0,imax+2,0,jmax+2);
	FREE_RMATRIX(f,0,imax+2,0,jmax+2);
	FREE_RMATRIX(g,0,imax+2,0,jmax+2);
	FREE_RMATRIX(rhs,0,imax+2,0,jmax+2);
	/*post for visualization*/
    return 0;
}
