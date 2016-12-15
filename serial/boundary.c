#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "boundary.h"
void setbound(double **u,double **v,int imax,int jmax,int wW, int wE,int wN,int wS){
	int j,i;
	int us=1;
	for(j=1;j<jmax+1;j++){
		u[j][0]=0;
		u[j][imax]=0;
		v[j][0]=-v[j][1];
		v[j][imax+1]=-v[j][imax];
	}
	for(i=1;i<imax+1;i++){
		v[0][i]=0;
		v[jmax][i]=0;
		u[0][i]=-u[1][i];
		u[jmax+1][i]=2*us-u[jmax][i];
	}
    /*printf("setbnd test u:%f\n",u[jmax+1][64]);*/
	return;
}

void setspecbcond(double **u,double **v,int imax,int jmax,char **problem){
	return;
}
