#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "init.h"
#include "uvp.h"
void comp_fg(double **u,double **v,double **f,double **g, int imax,int jmax,double delt,double delx,double dely,double gx,double gy,double gamma,double Re){
	int j,i;
	double a,b,c,d,e,ff,gg,h,va,vb,u2x,uvy,u2x2,u2y2;
	double ua,ub,uvx,v2y,v2x2,v2y2;
    for(j=0;j<jmax+2;j++){
        for(i=0;i<imax+2;i++){
            f[j][i]=0;
            g[j][i]=0;
        }
    }
	for(j=1;j<jmax+1;j++){
		f[j][0]=u[j][0];
		f[j][imax]=u[j][imax];
	}
	for(i=1;i<imax+1;i++){
		g[0][i]=v[0][i];
		g[jmax][i]=v[jmax][i];
	}
	for(j=1;j<jmax+1;j++){
		for(i=1;i<imax;i++){
			a=u[j][i]+u[j][i+1];
			b=u[j][i-1]+u[j][i];
			c=u[j][i]-u[j][i+1];
			d=u[j][i-1]-u[j][i];
			e=u[j][i]+u[j+1][i];
			ff=u[j-1][i]+u[j][i];
			gg=u[j][i]-u[j+1][i];
			h=u[j-1][i]-u[j][i];
			va=v[j][i]+v[j][i+1];
			vb=v[j-1][i]+v[j-1][i+1];
			u2x=1/delx*((a/2)*(a/2)-(b/2)*(b/2))
			+gamma*1/delx*(abs(a)/2*c/2-abs(b)/2*d/2);
			uvy=1/dely*(va/2*e/2-vb/2*ff/2)
			+gamma*1/dely*(abs(va)/2*gg/2-abs(vb)/2*h/2);
			u2x2=(u[j][i+1]-2*u[j][i]+u[j][i-1])/(delx*delx);
			u2y2=(u[j+1][i]-2*u[j][i]+u[j-1][i])/(dely*dely);
			f[j][i]=u[j][i]+delt*(1/Re*(u2x2+u2y2)-u2x-uvy+gx);
		}
	}
	for(j=1;j<jmax;j++){
		for(i=1;i<imax+1;i++){
			a=v[j][i]+v[j][i+1];
			b=v[j][i-1]+v[j][i];
			c=v[j][i]-v[j][i+1];
			d=v[j][i-1]-v[j][i];
			e=v[j][i]+v[j+1][i];
			ff=v[j-1][i]+v[j][i];
			gg=v[j][i]-v[j+1][i];
			h=v[j-1][i]-v[j][i];
			ua=u[j][i]+u[j+1][i];
			ub=u[j][i-1]+u[j+1][i-1];
			uvx=1/delx*(ua/2*a/2-ub/2*b/2)
			+gamma*1/delx*(abs(ua)/2*c/2-abs(ub)/2*d/2);
			v2y=1/dely*((e/2)*(e/2)-(ff/2)*(ff/2))
			+gamma*1/dely*(abs(e)/2*gg/2-abs(ff)/2*h/2);
			v2x2=(v[j][i+1]-2*v[j][i]+v[j][i-1])/(delx*delx);
			v2y2=(v[j+1][i]-2*v[j][i]+v[j-1][i])/(dely*dely);
			g[j][i]=v[j][i]+delt*(1/Re*(v2x2+v2y2)-uvx-v2y+gy);

		}
	}
    /*printf("test f:%f\n",f[64][64]);
    printf("test g:%f\n",g[64][64]);*/
	return;
}
void comp_rhs(double **f, double **g,double **rhs,int imax,int jmax,double delt,double delx,double dely){
	int j,i;
    for(j=0;j<jmax+2;j++){
        for(i=0;i<imax+2;i++){
            rhs[j][i]=0;
        }
    }
	for(j=1;j<jmax+1;j++){
		for(i=1;i<imax+1;i++){
			rhs[j][i]=1/delt*((f[j][i]-f[j][i-1])/delx+(g[j][i]-g[j-1][i])/dely);
        }
    }
    /*printf("test rhs:%f\n",rhs[32][1]);*/
	return;
}
int poisson(double **p,double **rhs,int imax,int jmax,double delx,double dely,double eps,int itermax,double omg){
	int it,j,i,eiw,eie,ejs,ejn;
    double sum;
	double **r;
    double res;
	for(it=0;it<itermax;it++){
		for(j=1;j<jmax+1;j++){
			p[j][0]=p[j][1];
			p[j][imax+1]=p[j][imax];
		}
		for(i=1;i<imax+1;i++){
			p[0][i]=p[1][i];
			p[jmax+1][i]=p[jmax][i];
		}
		r=RMATRIX(0,jmax+1,0,imax+1); /* Is that right??*/
        for(j=0;j<jmax+2;j++){
            for(i=0;i<imax+2;i++){
                r[j][i]=0;
            }
        }
        sum=0;
		for(j=1;j<jmax+1;j++){
			for(i=1;i<imax+1;i++){
                eiw=1;eie=1;ejs=1;ejn=1;
                double a1 = (1-omg)*p[j][i];
                double a2 = omg/
					(
						(eie+eiw)/(delx*delx)
						+(ejn+ejs)/(dely*dely)
					);
				double aa1 = eie*p[j][i+1];
				double aa2 = eiw*p[j][i-1];
				double a4 = (aa1+aa2)/(delx*delx);
				double aa3 = ejn*p[j+1][i];
				double aa4 = ejs*p[j-1][i];
				double a5 = (aa3+aa4)/(dely*dely);
				double a6 = rhs[j][i];
				double a3 = (a4+a5-a6);
				p[j][i]=a1+a2*a3;

				r[j][i]=
				(eie*(p[j][i+1]-p[j][i])-eiw*(p[j][i]-p[j][i-1]))
				/(delx*delx)
				+
				(ejn*(p[j+1][i]-p[j][i])-ejs*(p[j][i]-p[j-1][i]))
				/(dely*dely)-rhs[j][i];
                sum=sum+r[j][i]*r[j][i];
			}
		}
        FREE_RMATRIX(r,0,jmax+1,0,imax+1);
        res=sqrt(sum/(imax*jmax));
        /*printf("res is %f\n",res);*/
        if(res<eps){
            printf("Converged...%f\n",res);
        	break;
        }
	}
    /*printf("pressure test:%f\n",p[jmax+2][64]);*/
    /*printf("number of iteration:%d\n",it);*/
	return it;
}

void adap_uv(double **u,double **v,double **f,double **g,double **p,int imax,int jmax,double delt,double delx,double dely){
    int i,j;
    for(j=1;j<jmax+1;j++){
    	for(i=1;i<imax;i++){
    		u[j][i]=f[j][i]-delt/delx*(p[j][i+1]-p[j][i]);
    	}
    }
    for(j=1;j<jmax;j++){
    	for(i=1;i<imax+1;i++){
    		v[j][i]=g[j][i]-delt/dely*(p[j+1][i]-p[j][i]);
    	}
    }
    /*printf("adap test u:%f\n",u[64][64]);
    printf("adap test v:%f\n",v[64][64]);*/
	return;
}
