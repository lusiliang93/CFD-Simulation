#ifndef _UVP_H
#define _UVP_H
void comp_fg(float **u,float **v,float **f,float **g, int imax,int jmax,float delt,float delx,float dely,float gx,float gy,float gamma,float Re);
void comp_rhs(float **f, float **g,float **rhs,int imax,int jmax,float delt,float delx,float dely);
int poisson(float **p,float **rhs,int imax,int jmax,float delx,float dely,float eps,int itermax,float omg,float res);
void adap_uv(float **u,float **v,float **f,float **g,float **p,int imax,int jmax,float delt,float delx,float dely);
#endif