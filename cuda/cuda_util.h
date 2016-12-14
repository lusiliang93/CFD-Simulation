void setbound(double *u,double *v,int imax,int jmax,int wW, int wE,int wN,int wS);
void cuda_init(int imax, int jmax);

void comp_fg(int imax, int jmax,double delt,double delx,double dely,double gx,double gy,double gamma,double Re);
void comp_rhs(int imax, int jmax,double delt,double delx,double dely);
int poisson(int imax, int jmax,double delx,double dely,double eps,int itermax,double omg);
void init_uvp(int imax, int jmax,int UI, int VI, int PI);
void adap_uv(int imax, int jmax, double delt, double delx, double dely);