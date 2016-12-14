#ifndef _UVP_H
#define _UVP_H
void comp_fg(double **u,double **v,double **f,double **g, int imax,int jmax,double delt,double delx,double dely,double gx,double gy,double gamma,double Re,int rid,int cid,int iproc,int jproc);
void comp_rhs(double **f, double **g,double **rhs,int imax,int jmax,double delt,double delx,double dely,int procID);
int poisson(double **p,double **rhs,int imax,int jmax,double delx,double dely,double eps,int itermax,double omg,int iw,int ie,int js,int jn,
	int tw,int te,int ts,int tn,int rid,int cid,int iproc,int jproc,int procID,int nproc);
void adap_uv(double **u,double **v,double **f,double **g,double **p,int imax,int jmax,double delt,double delx,double dely,
	int tw,int te,int ts,int tn,int nproc,int rid,int cid,int iproc,int jproc,int procID);
#endif
