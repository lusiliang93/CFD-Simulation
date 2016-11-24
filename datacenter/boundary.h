#ifndef _BOUNDARY_H
#define _BOUNDARY_H
void setbound(double **u,double **v,int imax,int jmax,int wW, int wE,int wN,int wS,double uin,double vin,int iB,int iC,int iD,int iG,int jI);
void setspecbcond(double **u,double **v,int imax,int jmax,char **problem);
#endif
