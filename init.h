#ifndef _INIT_H
#define _INIT_H
int read_parameter(char *inputfile);
void init_uvp(float **u,float **v,float **p,int imax,int jmax,float UI,float VI,float PI);
float **RMATRIX(int nrl,int nrh, int ncl,int nch);
void FREE_RMATRIX(float **m,int nrl,int nrh,int ncl,int nch);
#endif
