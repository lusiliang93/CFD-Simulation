#ifndef _UVP_H
#define _UVP_H
void comp_fg(double **u,double **v,double **f,double **g, int imax,int jmax,double delt,double delx,double dely,double gx,double gy,double gamma,double Re,int rid,int cid);
void comp_rhs(double **f, double **g,double **rhs,int imax,int jmax,double delt,double delx,double dely);
int poisson(double **p,double **rhs,int imax,int jmax,double delx,double dely,double eps,int itermax,double omg,int iw,int ie,int is,int in,
	int tw,int te,int ts,int tn,int rid,int cid,int procID,int nproc);
void adap_uv(double **p,double **rhs,int imax,int jmax,double delx,double dely,double eps,int itermax,double omg,int iw,int ie,int is,int in,
	int tw,int te,int ts,int tn,int rid,int cid,int iproc,int jproc);
#endif
