#ifndef _BOUNDARY_H
#define _BOUNDARY_H
void setbound(double **u,double **v,int imax,int jmax,int wW, int wE,int wN,int wS);
void setspecbcond(double **u,double **v,int imax,int jmax,char **problem);
#endif
