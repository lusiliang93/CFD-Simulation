void cuda_init(int imax, int jmax);
void copy_matrix(int imax, int jmax);

void init_uvp(int imax, int jmax,int UI, int VI, int PI);

void setbound(int imax,int jmax,int wW, int wE,int wN,int wS);

double comp_delt(int imax, int jmax,double delx,double dely,double Re,double tau);

void comp_fg(int imax, int jmax,double delt,double delx,double dely,double gx,double gy,double gamma,double Re);
void comp_rhs(int imax, int jmax,double delt,double delx,double dely);
int poisson(int imax, int jmax,double delx,double dely,double eps,int itermax,double omg);
void adap_uv(int imax, int jmax, double delt, double delx, double dely);
