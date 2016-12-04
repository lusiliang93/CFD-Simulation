#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "init.h"
#include "uvp.h"
#include "mpi.h"
void comp_fg(double **u,double **v,double **f,double **g, int imax,int jmax,double delt,double delx,double dely,double gx,double gy,double gamma,double Re,int rid,int cid){
	int j,i;
	double a,b,c,d,e,ff,gg,h,va,vb,u2x,uvy,u2x2,u2y2;
	double ua,ub,uvx,v2y,v2x2,v2y2;
    // for(j=0;j<jmax+2;j++){
    //      for(i=0;i<imax+2;i++){
    //          f[j][i]=0;
    //          g[j][i]=0;
    //      }
    // }
    /* west */
    if(cid==0){
    	for(j=1;j<jmax+1;j++){
    		f[j][1]=u[j][1]; /* is that right?*/	 
    	}
    }
    /* east */
    if(cid==iproc-1){
    	for(j=1;j<jmax+1;j++){
    		f[j][imax+1]=u[j][imax+1];
    	}
    }
    /* south */
    if(rid==0){
    	for(i=1;i<imax+1;i++){
    		g[1][i]=v[1][i];		
    	}
    }
    /* north */
    if(rid==jproc-1){
    	for(i=1;i<imax+1;i++){
    		g[jmax+1][i]=v[jmax+1][i];
    	}
    }

	for(j=1;j<jmax+1;j++){ /* changed to fit the dimension */
		for(i=1;i<imax;i++){
			a=u[j][i+1]+u[j][i+1+1];
			b=u[j][i-1+1]+u[j][i+1];
			c=u[j][i+1]-u[j][i+1+1];
			d=u[j][i-1+1]-u[j][i+1];
			e=u[j][i+1]+u[j+1][i+1];
			ff=u[j-1][i+1]+u[j][i+1];
			gg=u[j][i+1]-u[j+1][i+1];
			h=u[j-1][i+1]-u[j][i+1];
			va=v[j+1][i]+v[j+1][i+1];
			vb=v[j-1+1][i]+v[j-1+1][i+1];
			u2x=1/delx*((a/2)*(a/2)-(b/2)*(b/2))
			+gamma*1/delx*(abs(a)/2*c/2-abs(b)/2*d/2);
			uvy=1/dely*(va/2*e/2-vb/2*ff/2)
			+gamma*1/dely*(abs(va)/2*gg/2-abs(vb)/2*h/2);
			u2x2=(u[j][i+1+1]-2*u[j][i+1]+u[j][i-1+1])/(delx*delx);
			u2y2=(u[j+1][i+1]-2*u[j][i+1]+u[j-1][i+1])/(dely*dely);
			f[j][i+1]=u[j][i+1]+delt*(1/Re*(u2x2+u2y2)-u2x-uvy+gx);
		}
	}
	for(j=1;j<jmax;j++){
		for(i=1;i<imax+1;i++){
			a=v[j+1][i]+v[j+1][i+1];
			b=v[j+1][i-1]+v[j+1][i];
			c=v[j+1][i]-v[j+1][i+1];
			d=v[j+1][i-1]-v[j+1][i];
			e=v[j+1][i]+v[j+1+1][i];
			ff=v[j-1+1][i]+v[j+1][i];
			gg=v[j+1][i]-v[j+1+1][i];
			h=v[j-1+1][i]-v[j+1][i];
			ua=u[j][i+1]+u[j+1][i+1];
			ub=u[j][i-1+1]+u[j+1][i-1+1];
			uvx=1/delx*(ua/2*a/2-ub/2*b/2)
			+gamma*1/delx*(abs(ua)/2*c/2-abs(ub)/2*d/2);
			v2y=1/dely*((e/2)*(e/2)-(ff/2)*(ff/2))
			+gamma*1/dely*(abs(e)/2*gg/2-abs(ff)/2*h/2);
			v2x2=(v[j+1][i+1]-2*v[j+1][i]+v[j+1][i-1])/(delx*delx);
			v2y2=(v[j+1+1][i]-2*v[j+1][i]+v[j-1+1][i])/(dely*dely);
			g[j+1][i]=v[j+1][i]+delt*(1/Re*(v2x2+v2y2)-uvx-v2y+gy);

		}
	}
	return;
}
void comp_rhs(double **f, double **g,double **rhs,int imax,int jmax,double delt,double delx,double dely){
	int j,i;
    // for(j=0;j<jmax+2;j++){
    //     for(i=0;i<imax+2;i++){
    //         rhs[j][i]=0;
    //     }
    // }
    /* changed to fit the dimension*/
	for(j=1;j<jmax+1;j++){
		for(i=1;i<imax+1;i++){
			rhs[j][i]=1/delt*((f[j][i+1]-f[j][i-1+1])/delx+(g[j+1][i]-g[j-1+1][i])/dely);
        }
    }
	return;
}
int poisson(double **p,double **rhs,int imax,int jmax,double delx,double dely,double eps,int itermax,double omg,int iw,int ie,int is,int in,
	int tw,int te,int ts,int tn,int rid,int cid,int iproc,int jproc){
	int it,j,i,eiw,eie,ejs,ejn,partial_sum,sender;
    double sum;
	double r;
    double res;
    double *w,*e,*s,*n;
    MPI_Status stauts;
    w_s=(double*)malloc(jmax*sizeof(double));
    e_s=(double*)malloc(jmax*sizeof(double));
    s_s=(double*)malloc(imax*sizeof(double));
    n_s=(double*)malloc(imax*sizeof(double));
	for(it=0;it<itermax;it++){ 
		/* west */
		if(cid==0){
			for(j=1;j<jmax+1;j++){
				p[j][0]=p[j][1];
			}
		}
		/* east */
		if(cid==iproc-1){
			for(j=1;j<jmax+1;j++){
				p[j][imax+1]=p[j][imax];
			}
	    }
	    /* south */
	    if(rid==0){
	    	for(i=1;i<imax+1;i++){
	    		p[0][i]=p[1][i];		
	    	}
	    }
	    /* north */
	    if(rid==jproc-1){
	    	for(i=1;i<imax+1;i++){
	    		p[jmax+1][i]=p[jmax][i];
	    	}
	    }
	    /* subdomain exchange */	    
	    if(cid!=0){
	    	/* send to the west(left) */
	    	for(j=1;j<jmax+1;j++){
	    		w_s[j-1]=p[j][1];
	    	}
	    	MPI_Send(&w_s,jmax,MPI_DOUBLE,tw,0,MPI_COMM_WORLD);
	    	/* receive from the west(right) */
	    	MPI_Recv(&w_r,jmax,MPI_DOUBLE,tw,0,MPI_COMM_WORLD,&status);
	    	for(j=1;j<jmax+1;j++){
	    		p[j][0]=w_r[j-1];
	    	}

	    }
	    if(cid!=iproc-1){
	    	/*send to the east(right)*/
	    	for(j=1;j<jmax+1;j++){
	    		e_s[j-1]=p[j][imax];
	    	}
	    	MPI_Send(&e_s,jmax,MPI_DOUBLE,te,0,MPI_COMM_WORLD);
	    	/* receive form the east(left) */
	    	MPI_Recv(&e_r,jmax,MPI_DOUBLE,te,0,MPI_COMM_WORLD,&status);
	    	for(j=1;j<jmax+1;j++){
	    		p[j][imax+1]=w_r[j-1];
	    	}
	    }
	    if(rid!=0){
	    	/*send to the south(bottom)*/
	    	for(i=1;i<imax+1;i++){
	    		s_s[i-1]=p[1][i];
	    	}
	    	MPI_Send(&s_s,imax,MPI_DOUBLE,ts,0,MPI_COMM_WORLD);
	    	/*receive from the south(top) */
	    	MPI_Recv(&s_r,imax,MPI_DOUBLE,ts,0,MPI_COMM_WORLD,&status);
	    	for(i=1;i<imax+1;i++){
	    		p[jmax+1][i]=s_r[i-1];
	    	}
	    }
	    if(rid!=jproc-1){
	    	/*send to the north(top)*/
	    	for(i=1;i<imax+1;i++){
	    		n_s[i-1]=p[jmax][i];
	    	}
	    	MPI_Send(&n_s,imax,MPI_DOUBLE,tn,0,MPI_COMM_WORLD);
	    	/*receive from the north(bottom) */
	    	MPI_Recv(&n_r,imax,MPI_DOUBLE,tn,0,MPI_COMM_WORLD,&status);
	    	for(i=1;i<imax+1;i++){
	    		p[0][i]=n_r[i-1];
	    	}
	    }

        sum=0;
		for(j=1;j<jmax+1;j++){
			for(i=1;i<imax+1;i++){
                eiw=1;eie=1;ejs=1;ejn=1;
				p[j][i]=(1-omg)*p[j][i]
				+omg/((eie+eiw)/(delx*delx)+(ejn+ejs)/(dely*dely))
				*((eie*p[j][i+1]+eiw*p[j][i-1])/(delx*delx)
				+(ejn*p[j+1][i]+ejs*p[j-1][i])/(dely*dely)-rhs[j][i]);
				r=(eie*(p[j][i+1]-p[j][i])-eiw*(p[j][i]-p[j][i-1]))/(delx*delx)
				+(ejn*(p[j+1][i]-p[j][i])-ejs*(p[j][i]-p[j-1][i]))/(dely*dely)-rhs[j][i];
                sum=sum+r*r;
			}
		}
		/* is that right???*/
		if(procID==0){
			for(i=1;i<nproc;i++){
				MPI_Recv(&partial_sum,1,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WOLRD,&status);
				sender=status.MPI_SOURCE;
				printf("partial sum %d returned from prcessor %d\n",partial_sum,sender);
				sum+=partial_sum;
			}
			res=sqrt(sum/(imax*jmax*nproc));
			MPI_Bcast(&res,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&iter,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		}
		else 
			MPI_Send(&sum,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);

	    if(res<eps){
	    	printf("Converged...%f\n",res);
			break;
		} 
	}
    free(e_s);
    free(w_s);
    free(s_s);
    free(n_s);
    free(e_r);
    free(w_r);
    free(n_r);
    free(s_r);
	return it;
}

void adap_uv(double **u,double **v,double **f,double **g,double **p,int imax,int jmax,double delt,double delx,double dely
	int tw,int te,int ts,int tn,int nproc,int rid,int cid,int iproc,int jproc,int tw,int te,int ts,int tn){
    int i,j;
    double *v_w_s,*v_e_s,*v_s_s,*v_n_s;
    double *v_w_r,*v_e_r,*v_s_r,*v_n_r;
    double *u_w_s,*u_e_s,*u_s_s,*u_n_s;
    double *u_w_r,*u_e_r,*u_s_r,*u_n_r;
    MPI_Status status;
    /* changed to fit the dimension*/
    for(j=1;j<jmax+1;j++){
    	for(i=1;i<imax;i++){
    		u[j][i+1]=f[j][i+1]-delt/delx*(p[j][i+1]-p[j][i]);
    	}
    }
    for(j=1;j<jmax;j++){
    	for(i=1;i<imax+1;i++){
    		v[j+1][i]=g[j+1][i]-delt/dely*(p[j+1][i]-p[j][i]);
    	}
    }
    /* subdomain exchange */
    /* west */
    if(cid!=0){
    	/*send to the west(left) */
    	for(j=1;j<jmax+2;j++){
    		v_w_s[j-1]=v[j][1];
    		if(j<=jmax){
    			u_w_s[j-1]=u[j][2];
    		}
    	}
    	MPI_Send(&v_w_s,jmax+1,MPI_DOUBLE,tw,0,MPI_COMM_WORLD);
    	MPI_Send(&u_w_s,jmax,MPI_DOUBLE,tw,0,MPI_COMM_WORLD);
    	/* receive from the west(right) */
    	/* how to make sure which one is which?*/
    	MPI_Recv(&v_w_r,jmax+1,MPI_DOUBLE,tw,0,MPI_COMM_WORLD,&status);
    	MPI_Recv(&u_w_r,jmax,MPI_DOUBLE,tw,0,MPI_COMM_WORLD,&status);
    	for(j=1;j<jmax+2;j++){
    		v[j][0]=v_w_r[j-1];
    		if(j<=jmax){
    			u[j][0]=u_w_r[j-1];
    		}
    	}
    }
    /* east */
    if(cid!=iproc-1){
    	/* send to east(right) */
    	for(j=1;j<jmax+2;j++){
    		v_e_s[j-1]=v[j][imax];
    		if(j<=jmax+1){
    			u_e_s[j-1]=u[j][imax]; /* is that right?*/
    		}
    	}
    	MPI_Send(&v_e_s,jmax+1,MPI_DOUBLE,te,0,MPI_COMM_WORLD);
    	MPI_Send(&u_e_s,jmax,MPI_DOUBLE,te,0,MPI_COMM_WORLD);
    	/* recevie from the east(left) */
    	MPI_Recv(&v_e_r,jmax+1,MPI_DOUBLE,te,0,MPI_COMM_WORLD,&status);
    	MPI_Recv(&u_e_r,jmax,MPI_DOUBLE,te,0,MPI_COMM_WORLD,&status);
    	for(j=1;j<jmax+2;j++){
    		v[j][imax+1]=v_e_r[j-1];
    		if(j<=jmax){
    			u[j][imax+2]=u_e_r[j-1]; /* is that right ?*/
    		}
    	}
    }
    /* south */
    /* according to west */
    if(rid!=0){
    	/* send to the south(down)*/
    	for(i=1;i<imax+2;i++){
    		u_s_s[i-1]=u[1][i];
    		if(i<=imax){
    			v_s_s[i-1]=v[2][i];
    		}
    	}
    	MPI_Send(&u_s_s,imax+1,MPI_DOUBLE,ts,0,MPI_COMM_WORLD);
    	MPI_Send(&v_s_s,imax,MPI_DOUBLE,ts,0,MPI_COMM_WORLD);
    	/* receive from the north(up)*/
    	MPI_Recv(&u_s_r,imax+1,MPI_DOUBLE,ts,0,MPI_COMM_WORLD,&status);
    	MPI_Recv(&v_s_r,imax,MPI_DOUBLE,ts,0,MPI_COMM_WORLD,&status);
    	for(i=1;i<imax+2;i++){
    		u[0][i]=u_s_r[i-1];
    		if(i<=imax){
    			v[0][i]=v_s_r[i-1];
    		}
    	}
    }
    /* north */
    /* according to east */
    if(rid!=jproc-1){
    	/* send to north(up) */
    	for(i=1;i<imax+2;i++){
    		u_n_s[i-1]=u[jmax][i];
    		if(i<=imax){
    			v_n_s[i-1]=v[jmax][i];
    		}
    	}
    	MPI_Send(&u_n_s,imax+1,MPI_DOUBLE,tn,0,MPI_COMM_WORLD);
    	MPI_Send(&v_n_s,imax,MPI_DOUBLE,tn,0,MPI_COMM_WORLD);
    	/* receive from the north(down)*/
    	MPI_Recv(&u_n_r,imax+1,MPI_DOUBLE,tn,0,MPI_COMM_WORLD,&status);
    	MPI_Recv(&v_n_r,imax,MPI_DOUBLE,tn,0,MPI_COMM_WORLD,&status);
    	for(i=1;i<imax+2;i++){
    		u[jmax+1][i]=u_n_r[i-1];
    		if(i<=imax){
    			v[jmax+2][i]=v_n_r[i-1];
    		}
    	}
    }
	return;
}
