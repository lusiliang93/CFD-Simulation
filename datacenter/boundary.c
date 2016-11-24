#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "boundary.h"
void setbound(double **u,double **v,int imax,int jmax,int wW, int wE,int wN,int wS,double uin,double vin,int iB,int iC,int iD,int iG,int jI){
	int j,i;
	/*int us=1;
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
		u[jmax+1][i]=2*us-u[jmax][i];//change
	}
	*/
	/* BC,no-slip*/
	for(i=iB+1;i<iC+1;i++){
		u[0][i]=-u[1][i];
		v[0][i]=0;
	}
	/* CD, inflow*/
	for(i=iC+1;i<iD+1;i++){
		u[0][i]=2*uin-u[1][i];
		v[0][i]=vin;
	}
	/* DE,no-slip*/
	for(i=iD+1;i<imax+1;i++){
		u[0][i]=-u[1][i];
		v[0][i]=0;
	}
	/*EF,no-slip*/
	for(j=1;j<jmax+1;j++){
		u[j][imax]=0;
		v[j][imax+1]=-v[j][imax];
	}
	/*GF,no-slip*/
	for(i=iG+1;i<imax+1;i++){
		u[jmax+1][i]=-u[jmax][i];
		v[jmax][i]=0;
	}
	/*HG,outflow*/
	for(i=1;i<iG+1;i++){
		u[jmax+1][i]=u[jmax][i];
		v[jmax][i]=v[jmax-1][i];
	}
	/*HI,symmetry*/
	for(j=jI+1;j<jmax+1;j++){
		u[j][0]=0;
		v[j][0]=v[j][1];
	}
	/*IJ,no-slip,server top*/
	for(i=1;i<iB+1;i++){
		u[jI][i]=-u[jI+1][i];
		v[jI][i]=0;
	}
	/* JB,no-slip,server right*/
	for(j=1;j<jI+1;j++){
		u[j][iB]=0;
		v[j][iB]=-v[j][iB+1];
	}
	return;
}

void setspecbcond(double **u,double **v,int imax,int jmax,char **problem){
	return;
}
