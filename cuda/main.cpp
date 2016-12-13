#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "cuda_util.h"

int imax, jmax;

#define get_index(i,j) ((jmax+2)*i+j)

int read_parameter(char *inputfile){
	return 0;
}

double *RMATRIX(int nrl,int nrh, int ncl,int nch){
/* reserve memory for matrix of size [nrl,nrh]x[ncl,nch]*/
	double* m;
	m = (double*)malloc((unsigned)(nrh-nrl+1)*(nch-ncl+1)*sizeof(double));
	if(m==NULL){
		printf("no more memory \n");
		exit(0);
	}
	return m;
}
void init_uvp(double *u,double *v,double *p,double UI,double VI,double PI){
	int i,j;
	for(j=0;j<jmax+2;j++){
		for(i=0;i<imax+2;i++){
			u[get_index(j,i)] = UI;
			v[get_index(j,i)] = VI;
			p[get_index(j,i)] = PI;
		}
	}
	return;
}

void FREE_RMATRIX(double *m){
	free(m);
    return;
}

void comp_fg(double *u,double *v,double *f,double *g,double delt,double delx,double dely,double gx,double gy,double gamma,double Re){
	int j,i;
	double a,b,c,d,e,ff,gg,h,va,vb,u2x,uvy,u2x2,u2y2;
	double ua,ub,uvx,v2y,v2x2,v2y2;
    for(j=0;j<jmax+2;j++){
        for(i=0;i<imax+2;i++){
        	f[get_index(j,i)] = 0;
        	g[get_index(j,i)] = 0;
        }
    }
	for(j=1;j<jmax+1;j++){
		f[get_index(j,0)] = u[get_index(j,0)];
		f[get_index(j,imax)] = u[get_index(j,imax)];
	}
	for(i=1;i<imax+1;i++){
		g[get_index(0,i)] = v[get_index(0,i)];
		g[get_index(jmax,i)] = v[get_index(jmax,i)];
	}
	for(j=1;j<jmax+1;j++){
		for(i=1;i<imax;i++){
			a = u[get_index(j,i)] + u[get_index(j,i+1)];
			b = u[get_index(j,i-1)] + u[get_index(j,i)];
			c = u[get_index(j,i)] - u[get_index(j,i+1)];
			d = u[get_index(j,i-1)] - u[get_index(j,i)];
			e = u[get_index(j,i)] + u[get_index(j+1,i)];
			ff = u[get_index(j-1,i)] + u[get_index(j,i)];
			gg = u[get_index(j,i)] - u[get_index(j+1,i)];
			h = u[get_index(j-1,i)] - u[get_index(j,i)];

			va = v[get_index(j,i)] + v[get_index(j,i+1)];
			vb = v[get_index(j-1,i)] + v[get_index(j-1,i+1)];

			u2x = 1/delx*((a/2)*(a/2)-(b/2)*(b/2))+gamma*1/delx*(abs(a)/2*c/2-abs(b)/2*d/2);
			uvy = 1/dely*(va/2*e/2-vb/2*ff/2)+gamma*1/dely*(abs(va)/2*gg/2-abs(vb)/2*h/2);
			u2x2 = (u[get_index(j,i+1)] - 2*u[get_index(j,i)] + u[get_index(j,i-1)])/(delx*delx);
			u2y2 = (u[get_index(j+1,i)] - 2*u[get_index(j,i)] + u[get_index(j-1,i)])/(dely*dely);

			f[get_index(j,i)] = u[get_index(j,i)] + delt*(1/Re*(u2x2+u2y2)-u2x-uvy+gx);
		}
	}
	for(j=1;j<jmax;j++){
		for(i=1;i<imax+1;i++){
			a = v[get_index(j,i)] + v[get_index(j,i+1)];
			b = v[get_index(j,i-1)] + v[get_index(j,i)];
			c = v[get_index(j,i)] - v[get_index(j,i+1)];
			d = v[get_index(j,i-1)] - v[get_index(j,i)];
			e = v[get_index(j,i)] + v[get_index(j+1,i)];
			ff = v[get_index(j-1,i)] + v[get_index(j,i)];
			gg = v[get_index(j,i)] - v[get_index(j+1,i)];
			h = v[get_index(j-1,i)] - v[get_index(j,i)];

			ua = u[get_index(j,i)] + u[get_index(j+1,i)];
			ub = u[get_index(j,i-1)] + u[get_index(j+1,i-1)];

			uvx = 1/delx*(ua/2*a/2-ub/2*b/2)+gamma*1/delx*(abs(ua)/2*c/2-abs(ub)/2*d/2);
			v2y = 1/dely*((e/2)*(e/2)-(ff/2)*(ff/2))+gamma*1/dely*(abs(e)/2*gg/2-abs(ff)/2*h/2);
			v2x2 = (v[get_index(j,i+1)] - 2*v[get_index(j,i)] + v[get_index(j,i-1)])/(delx*delx);
			v2y2 = (v[get_index(j+1,i)] - 2*v[get_index(j,i)] + v[get_index(j-1,i)])/(dely*dely);

			g[get_index(j,i)] = v[get_index(j,i)] + delt*(1/Re*(v2x2+v2y2)-uvx-v2y+gy);
		}
	}
	return;
}

void comp_rhs(double *f, double *g,double *rhs,double delt,double delx,double dely){
	int j,i;
    for(j=0;j<jmax+2;j++){
        for(i=0;i<imax+2;i++){
        	rhs[get_index(j,i)] = 0;
        }
    }
	for(j=1;j<jmax+1;j++){
		for(i=1;i<imax+1;i++){
			int tmp = (f[get_index(j,i)]-f[get_index(j,i-1)])/delx + (g[get_index(j,i)]-g[get_index(j-1,i)])/dely;
			rhs[get_index(j,i)] = 1/delt * tmp;
        }
    }
	return;
}

int poisson(double *p,double *rhs,double delx,double dely,double eps,int itermax,double omg){
	int it,j,i,eiw,eie,ejs,ejn;
    double sum;
	double *r;
    double res;
	for(it=0;it<itermax;it++){
		for(j=1;j<jmax+1;j++){
			p[get_index(j,0)] = p[get_index(j,1)];
			p[get_index(j,imax+1)] = p[get_index(j,imax)];
		}
		for(i=1;i<imax+1;i++){
			p[get_index(0,i)] = p[get_index(1,i)];
			p[get_index(jmax+1,i)] = p[get_index(jmax,i)];
		}
		r=RMATRIX(0,jmax+1,0,imax+1);
        for(j=0;j<jmax+2;j++){
            for(i=0;i<imax+2;i++){
                r[get_index(j,i)] = 0;
            }
        }
        sum=0;
		for(j=1;j<jmax+1;j++){
			for(i=1;i<imax+1;i++){
                eiw=1;eie=1;ejs=1;ejn=1;
                p[get_index(j,i)] = (1-omg)*p[get_index(j,i)]
                +
                omg/((eie+eiw)/(delx*delx)+(ejn+ejs)/(dely*dely)) * (
                	(eie*p[get_index(j,i+1)]+eiw*p[get_index(j,i-1)])/(delx*delx)
                	+(ejn*p[get_index(j+1,i)]+ejs*p[get_index(j-1,i)])/(dely*dely)
                	-rhs[get_index(j,i)]
                );

                r[get_index(j,i)] = (
                	eie*(p[get_index(j,i+1)]-p[get_index(j,i)])
                	-eiw*(p[get_index(j,i)]-p[get_index(j,i-1)])
                	)/(delx*delx)
                +	(
                	ejn*(p[get_index(j+1,i)]-p[get_index(j,i)])
                	-ejs*(p[get_index(j,i)]-p[get_index(j-1,i)])
                	)/(dely*dely)
                - rhs[get_index(j,i)];
                sum += r[get_index(j,i)]*r[get_index(j,i)];
			}
		}
        FREE_RMATRIX(r);
        res=sqrt(sum/(imax*jmax));
        if(res<eps){
            printf("Converged...%f\n",res);
        	break;
        }
	}
	return it;
}

void adap_uv(double *u,double *v,double *f,double *g,double *p,double delt,double delx,double dely){
    int i,j;
    for(j=1;j<jmax+1;j++){
    	for(i=1;i<imax;i++){
    		u[get_index(j,i)] = f[get_index(j,i)] - delt/delx*(p[get_index(j,i+1)]-p[get_index(j,i)]);
    	}
    }
    for(j=1;j<jmax;j++){
    	for(i=1;i<imax+1;i++){
    		v[get_index(j,i)] = g[get_index(j,i)] - delt/dely*(p[get_index(j+1,i)]-p[get_index(j,i)]);
    	}
    }
	return;
}


double max(double* u){
	int i,j;
	double max=0;
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
			if(u[get_index(j,i)]>max){
				max = u[get_index(j,i)];
			}
		}
	}
	return max;
}

void comp_delt(double* delt,double delx,double dely,double *u,double *v,double Re,double tau){
    double first,second,third,min;
	double delta = 1/(delx*delx)+1/(dely*dely);
	first = Re/2/delta;
	min=first;
	second = delx/abs(max(u));
	third= dely/abs(max(v));
	if(min>second){
		min=second;
		if(min>third)
			min=third;
	}
	else{
		if(min>third)
			min=third;
	}
	*delt=tau*min;
	return;
}

int main(int argc,char* argv[]){
	int n=0;
	double t=0;
	double xlength,ylength;
	double tend,tau,itermax,eps,omg,gamma;
	double Re,GX,GY,UI,VI,PI;
	int wW,wE,wS,wN;
	double x,y,x1,y1,x2,y2,u1,u2,u3,u4,v1,v2,v3,v4;
	double *u;
	double *v;
	double *p;
	double *f;
	double *g;
	double *rhs;
	double *uu;
	double *vv;
	double *xx;
	double *yy;
	double delx,dely,delt;
	int i,j,ii,jj;
    FILE *input;
    FILE *outputu;
    FILE *outputv;
    FILE *outputu1;
    FILE *outputv1;
    clock_t t1,t2;
    double  total_t;
    input = fopen("input.txt","r");
	if(!input){
		printf("Unable to open the file.\n");
		return -1;
	}

	fscanf(input,"%lf %lf\n",&xlength,&ylength);
	fscanf(input,"%d %d\n",&imax,&jmax);
	fscanf(input,"%lf %lf %lf %lf %lf %lf\n",&Re,&UI,&VI,&PI,&GX,&GY);
	/* tau must be postive*/
	fscanf(input,"%lf %lf %lf %lf %lf %lf \n",&tend,&tau,&itermax,&eps,&omg,&gamma);
	fscanf(input,"%d %d %d %d\n",&wW,&wE,&wN,&wS);
	delx=xlength/imax;
	dely=ylength/jmax;
    printf("xlengh:%f ylength:%f jmax:%d imax:%d Re:%f UI:%f VI:%f PI:%f GX:%f GY:%f tend:%f tau:%f itermax:%f eps:%f omg:%f gamma:%f\n",xlength,ylength,jmax,imax,Re,UI,VI,PI,GX,GY,tend,tau,itermax,eps,omg,gamma);

	/* assign initial values to u,v,p,f,g,rhs,uu,vv*/
	u=RMATRIX(0,imax+1,0,jmax+1);
	v=RMATRIX(0,imax+1,0,jmax+1);
	p=RMATRIX(0,imax+1,0,jmax+1);
	f=RMATRIX(0,imax+1,0,jmax+1);
	g=RMATRIX(0,imax+1,0,jmax+1);
	rhs=RMATRIX(0,imax+1,0,jmax+1);
	uu=RMATRIX(0,imax,0,jmax);
	vv=RMATRIX(0,imax,0,jmax);
	/* allocate memory to xx, yy*/
	xx=(double *)malloc((imax+1)*sizeof(double));
	yy=(double *)malloc((jmax+1)*sizeof(double));
	init_uvp(u,v,p,UI,VI,PI);
    t1=clock();
    cuda_init(imax, jmax);
	while(t<tend){
        if(n==0){
	        delt=0.02;
	        printf("before setbound\n");
			setbound(u,v,imax,jmax,wW,wE,wN,wS);
			comp_fg(u,v,f,g,delt,delx,dely,GX,GY,gamma,Re);
			comp_rhs(f, g,rhs,delt,delx,dely);
			poisson(p,rhs,delx,dely,eps,itermax,omg);
			adap_uv(u,v,f,g,p,delt,delx,dely);
			t=t+delt;
			n++;
			printf("The current t:%f\n",t);
        }else{
            comp_delt(&delt,delx,dely,u,v,Re,tau);
            printf("before setbound\n");
            setbound(u,v,imax,jmax,wW,wE,wN,wS);
            comp_fg(u,v,f,g,delt,delx,dely,GX,GY,gamma,Re);
            comp_rhs(f,g,rhs,delt,delx,dely);
            poisson(p,rhs,delx,dely,eps,itermax,omg);
            adap_uv(u,v,f,g,p,delt,delx,dely);
            t=t+delt;
            n++;
            printf("The current t:%f\n",t);
        }

	}
    t2=clock();
    total_t=(double)(t2-t1)/CLOCKS_PER_SEC;
    printf("Time elapsed:%f\n",total_t);

	outputu = fopen("outputu.txt","w+");
	outputv = fopen("outputv.txt","w+");

	for(j=0;j<jmax+2;j++){
		for(i=0;i<imax+2;i++){
			fprintf(outputu,"%f ", u[get_index(j,i)]);
		}
		fprintf(outputu,"\n");
	}
	for(j=0;j<jmax+2;j++){
		for(i=0;i<imax+2;i++){
            fprintf(outputv,"%f ", v[get_index(j,i)]);
		}
		fprintf(outputv,"\n");
	}
	fclose(outputu);
	fclose(outputv);
	printf("u into file\n");
	printf("v into file\n");

	/* post for visualization*/
	for(i=0;i<imax+1;i++){
		xx[i]=delx*i;
	}
	for(j=0;j<jmax+1;j++){
		yy[j]=dely*j;
	}
	for(j=0;j<jmax+1;j++){/*Is that right?*/
		for(i=0;i<imax+1;i++){/*Is that right?*/
			x=xx[i];/*Is that right?*/
			y=yy[j];/*Is that right?*/
			ii=floor(x/delx)+1;
			jj=floor((y+dely/2)/dely)+1;
			x1=(ii-1)*delx;
			y1=((jj-1)-0.5)*dely;
			x2=ii*delx;
			y2=(jj-1/2)*dely;
			u1 = u[get_index(jj-1,ii-1)];
			u2 = u[get_index(jj-1,ii)];
			u3 = u[get_index(jj,ii-1)];
			u4 = u[get_index(jj,ii)];
			uu[j*(jmax+1)+i] = 1/(delx*dely)*((x2-x)*(y2-y)*u1+(x-x1)*(y2-y)*u2+(x2-x)*(y-y1)*u3+(x-x1)*(y-y1)*u4);
		}
	}
	/* vv*/
	for(j=0;j<jmax+1;j++){/*Is that right?*/
		for(i=0;i<imax+1;i++){
			x=xx[i];
			y=yy[j];
			jj=floor(y/dely)+1;
			ii=floor((x+delx/2)/delx)+1;
			y1=(jj-1)*dely;
			x1=((ii-1)-0.5)*delx;
			y2=jj*dely;
			x2=(ii-0.5)*delx;
			v1 = v[get_index(jj-1,ii-1)];
			v2 = v[get_index(jj-1,ii)];
			v3 = v[get_index(jj,ii-1)];
			v4 = v[get_index(jj,ii)];
			vv[j*(jmax+1)+i] = 1/(delx*dely)*((x2-x)*(y2-y)*v1+(x-x1)*(y2-y)*v2+(x2-x)*(y-y1)*v3+(x-x1)*(y-y1)*v4);
		}
	}

	outputu1 = fopen("post_outputu.txt","w+");
	outputv1 = fopen("post_outputv.txt","w+");
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
			fprintf(outputu1,"%f ", uu[j*(jmax+1)+i]);
		}
		fprintf(outputu1,"\n");
	}
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
            fprintf(outputv1,"%f ", vv[j*(jmax+1)+i]);
		}
		fprintf(outputv1,"\n");
	}
	fclose(outputu1);
	fclose(outputv1);
	printf("uu into file\n");
	printf("vv into file\n");

	FREE_RMATRIX(u);
	FREE_RMATRIX(v);
	FREE_RMATRIX(p);
	FREE_RMATRIX(f);
	FREE_RMATRIX(g);
	FREE_RMATRIX(rhs);
	FREE_RMATRIX(uu);
	FREE_RMATRIX(vv);
	free(xx);
	free(yy);
    return 0;
}
