#ifndef _BOUNDARY_H
#define _BOUNDARY_H
void setbound(float **u,float **v,int imax,int jmax,int wW, int wE,int wN,int wS);
void setspecbcond(float **u,float **v,int imax,int jmax,char **problem);
void comp_fg(float **u,float **v,float **f,float **g, int imax,int jmax,float delt,float delx,float dely,float gx,float gy,float gamma,float Re);
#endif