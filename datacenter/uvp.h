#ifndef _UVP_H
#define _UVP_H
void comp_fg(double **u,double **v,double **f,double **g, int imax,int jmax,double delt,double delx,double dely,double gx,double gy,double gamma,double Re,int iB,int iC,int iD,int iG,int jI);
void comp_rhs(double **f, double **g,double **rhs,int imax,int jmax,double delt,double delx,double dely);
int poisson(double **p,double **rhs,int imax,int jmax,double delx,double dely,double eps,int itermax,double omg,int iB,int iC,int iD,int iG,int jI);
void adap_uv(double **u,double **v,double **f,double **g,double **p,int imax,int jmax,double delt,double delx,double dely);
#endif
