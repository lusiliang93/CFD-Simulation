#ifndef _BOUNDARY_H
#define _BOUNDARY_H
void setbound(double **u,double **v,int imax,int jmax,int wW, int wE,int wN,int wS);
void setspecbcond(double **u,double **v,int imax,int jmax,char **problem);
void comp_fg(double **u,double **v,double **f,double **g, int imax,int jmax,double delt,double delx,double dely,double gx,double gy,double gamma,double Re);
#endif
