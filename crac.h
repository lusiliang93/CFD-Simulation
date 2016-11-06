/**
 *  serial version of 2D Computaional Fluid Dynamics 
 *  Jin Hu, Siliang Lu
 */
#ifndef _CRAC_H
#define _CRAC_H

int read_parameter(char *inputfile);
void init_uvp(float **u,float **v,float **p,int imax,int jmax,float UI,float VI,float PI);
void comp_delt(float delt,int imax, int jmax,float delx,float dely, float **u, float **v, float Re, float tau);
void setbound(float **u,float **v,int imax,int jmax,int wW, int wE,int wN,int wS);
void comp_fg(float **u,float **v,float **f,float **g, int imax,int jmax,float delt,float delx,float dely,float gx,float gy,float gamma,float Re);
void comp_rhs(float **f, float **g,float **rhs,int imax,int jmax,float delt,float delx,float dely);
int poisson(float **p,float **rhs,int imax,int jmax,float delx,float dely,float eps,int itermax,float omg,float res);
void adap_uv(float **u,float **v,float **f,float **g,float **p,int imax,int jmax,float delt,float delx,float dely);
//float **rmatrix(int nrl, int nrh, int ncl, int nch);
#endif 