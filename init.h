#ifndef _INIT_H
#define _INIT_H
int read_parameter(char *inputfile);
void init_uvp(double **u,double **v,double **p,int imax,int jmax,double UI,double VI,double PI);
double **RMATRIX(int nrl,int nrh, int ncl,int nch);
void FREE_RMATRIX(double **m,int nrl,int nrh,int ncl,int nch);
#endif
