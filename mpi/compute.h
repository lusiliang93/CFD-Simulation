#ifndef _COMPUTE_H
#define _COMPUTe_H
double max(double** u,int imax,int jmax);
void comp_delt(double* delt,int imax,int jmax,double delx,double dely,double **u,double **v,double Re,double tau,int procID,int nproc,int iproc,int jproc);
int compute(int procID,int nproc,char *inputname);
#endif
