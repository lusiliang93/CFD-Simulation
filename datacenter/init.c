#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "init.h"
int read_parameter(char *inputfile){
	return 0;
}
double **RMATRIX(int nrl,int nrh, int ncl,int nch){
/** reserve memory for matrix of size [nrl,nrh]x[ncl,nch]**/
	int i;
	double **m;
	/** allocate row pointers **/
	if((m=(double**)malloc((unsigned)(nrh-nrl+1)*sizeof(double*)))==NULL){
		printf("no more memory \n");
		exit(0);
	}
	m -= nrl;
	/** allocate rows and set previously allocated row pointers to point to these**/
	for (i=nrl;i<=nrh;i++){
		if((m[i]=(double*)malloc((unsigned)(nch-ncl+1)*sizeof(double)))==NULL){
			printf("no more memory \n");
			exit(0);
		}
		m[i] -= ncl;
	}
	return m;
}
void init_uvp(double **u,double **v,double **p,int imax,int jmax,double UI,double VI,double PI){
	int i,j;
	for(j=0;j<jmax+2;j++){
		for(i=0;i<imax+2;i++){
			u[j][i]=UI;
			v[j][i]=VI;
			p[j][i]=PI;
		}
	}
	return;
}

void FREE_RMATRIX(double **m,int nrl,int nrh,int ncl,int nch){
/** frees memory of matrix allocated by RMATRIX **/
	int i;
	for(i=nrh;i>nrl;i--)
		free((char*)(m[i]+ncl));
	free((char*)(m+nrl));
    return;
}
