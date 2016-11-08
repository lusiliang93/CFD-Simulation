#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "uvp.h"
void comp_fg(float **u,float **v,float **f,float **g, int imax,int jmax,float delt,float delx,float dely,float gx,float gy,float gamma,float Re){
	return;
}
void comp_rhs(float **f, float **g,float **rhs,int imax,int jmax,float delt,float delx,float dely){
	return;
}
int poisson(float **p,float **rhs,int imax,int jmax,float delx,float dely,float eps,int itermax,float omg,float res){
	return 0;
}
void adap_uv(float **u,float **v,float **f,float **g,float **p,int imax,int jmax,float delt,float delx,float dely){
	return;
}