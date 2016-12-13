#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "cuda_util.h"

int read_parameter(char *inputfile){
	return 0;
}
double **RMATRIX(int nrl,int nrh, int ncl,int nch){
/** reserve memory for matrix of size [nrl,nrh]x[ncl,nch]**/
	int i;
	double **m;
	/** allocate row pointers **/
	if((m=(double**)malloc((unsigned)(nrh-nrl+1)*sizeof(double*)))==NULL){
		printf("no more memory \n");
		exit(0);
	}
	m -= nrl;
	/** allocate rows and set previously allocated row pointers to point to these**/
	for (i=nrl;i<=nrh;i++){
		if((m[i]=(double*)malloc((unsigned)(nch-ncl+1)*sizeof(double)))==NULL){
			printf("no more memory \n");
			exit(0);
		}
		m[i] -= ncl;
	}
	return m;
}
void init_uvp(double **u,double **v,double **p,int imax,int jmax,double UI,double VI,double PI){
	int i,j;
	for(j=0;j<jmax+2;j++){
		for(i=0;i<imax+2;i++){
			u[j][i]=UI;
			v[j][i]=VI;
			p[j][i]=PI;
		}
	}
	return;
}

void FREE_RMATRIX(double **m,int nrl,int nrh,int ncl,int nch){
/** frees memory of matrix allocated by RMATRIX **/
	int i;
	for(i=nrh;i>nrl;i--)
		free((char*)(m[i]+ncl));
	free((char*)(m+nrl));
    return;
}

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
				p[j][i]=(1-omg)*p[j][i]
				+omg/((eie+eiw)/(delx*delx)+(ejn+ejs)/(dely*dely))
				*((eie*p[j][i+1]+eiw*p[j][i-1])/(delx*delx)
				+(ejn*p[j+1][i]+ejs*p[j-1][i])/(dely*dely)-rhs[j][i]);

				r[j][i]=(eie*(p[j][i+1]-p[j][i])-eiw*(p[j][i]-p[j][i-1]))/(delx*delx)
				+(ejn*(p[j+1][i]-p[j][i])-ejs*(p[j][i]-p[j-1][i]))/(dely*dely)-rhs[j][i];
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


double max(double** u,int imax,int jmax){
	int i,j;
	double max=0;
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
			if (u[j][i]>max)
				max=u[j][i];
		}
	}
	return max;
}

void comp_delt(double* delt,int imax,int jmax,double delx,double dely,double **u,double **v,double Re,double tau){
    double first,second,third,min;
	double delta = 1/(delx*delx)+1/(dely*dely);
	first = Re/2/delta;
	min=first;
	second = delx/abs(max(u,imax,jmax));
	third= dely/abs(max(v,imax,jmax));
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
	int imax,jmax;
	int wW,wE,wS,wN;
	double x,y,x1,y1,x2,y2,u1,u2,u3,u4,v1,v2,v3,v4;
	double **u;
	double **v;
	double **p;
	double **f;
	double **g;
	double **rhs;
	double **uu;
	double **vv;
	double *xx;
	double *yy;
	double delx,dely,delt;
	int i,j,ii,jj;
    FILE *input;
    FILE *outputu;
    FILE *outputv;
    FILE *outputu1;
    FILE *outputv1;
    char *outputfilenameu;
    char *outputfilenamev;
    char *outputfilenameu1;
    char *outputfilenamev1;
    char *inputname=NULL;
    clock_t t1,t2;
    double  total_t;
	/*Read command line arguments*/
	inputname = "input.txt";

	/* read inputs*/
	if (inputname==NULL){
		printf("Inputname is incorrect.\n");
		return -1;
	}
    input = fopen(inputname,"r");
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
	init_uvp(u,v,p,imax,jmax,UI,VI,PI);
    t1=clock();
    cuda_init(imax, jmax);
	while(t<tend){
        if(n==0){
	        delt=0.02;
			setbound(u,v,imax,jmax,wW,wE,wN,wS);
			comp_fg(u,v,f,g, imax,jmax,delt,delx,dely,GX,GY,gamma,Re);
			comp_rhs(f, g,rhs,imax,jmax,delt,delx,dely);
			poisson(p,rhs,imax,jmax,delx,dely,eps,itermax,omg);
			adap_uv(u,v,f,g,p,imax,jmax,delt,delx,dely);
			t=t+delt;
			n++;
			printf("The current t:%f\n",t);
        }else{
            comp_delt(&delt,imax,jmax,delx,dely,u,v,Re,tau);
            /*printf("seg fault:%f\n",u[128][128]);*/
            setbound(u,v,imax,jmax,wW,wE,wN,wS);
            comp_fg(u,v,f,g,imax,jmax,delt,delx,dely,GX,GY,gamma,Re);
            comp_rhs(f,g,rhs,imax,jmax,delt,delx,dely);
            poisson(p,rhs,imax,jmax,delx,dely,eps,itermax,omg);
            adap_uv(u,v,f,g,p,imax,jmax,delt,delx,dely);
            t=t+delt;
            n++;
            /*printf("The current delt:%f\n",delt);*/
            printf("The current t:%f\n",t);
        }

	}
    t2=clock();
    total_t=(double)(t2-t1)/CLOCKS_PER_SEC;
    printf("Time elapsed:%f\n",total_t);

	outputfilenameu="outputu.txt";
    outputfilenamev="outputv.txt";
	outputu = fopen(outputfilenameu,"w+");
	outputv = fopen(outputfilenamev,"w+");

	for(j=0;j<jmax+2;j++){
		for(i=0;i<imax+2;i++){
			fprintf(outputu,"%f ",u[j][i]);
		}
		fprintf(outputu,"\n");
	}
	for(j=0;j<jmax+2;j++){
		for(i=0;i<imax+2;i++){
            fprintf(outputv,"%f ",v[j][i]);
		}
		fprintf(outputv,"\n");
	}
	fclose(outputu);
	fclose(outputv);
	printf("u into file:%s\n",outputfilenameu);
	printf("v into file:%s\n",outputfilenamev);

	/** post for visualization*/
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
			u1=u[jj-1][ii-1];
			u2=u[jj-1][ii];
			u3=u[jj][ii-1];
			u4=u[jj][ii];
			uu[j][i]=1/(delx*dely)*((x2-x)*(y2-y)*u1+(x-x1)*(y2-y)*u2+(x2-x)*(y-y1)*u3+(x-x1)*(y-y1)*u4);
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
			v1=v[jj-1][ii-1];
			v2=v[jj-1][ii];
			v3=v[jj][ii-1];
			v4=v[jj][ii];
			vv[j][i]=1/(delx*dely)*((x2-x)*(y2-y)*v1+(x-x1)*(y2-y)*v2+(x2-x)*(y-y1)*v3+(x-x1)*(y-y1)*v4);
		}
	}

	outputfilenameu1="post_outputu.txt";
    outputfilenamev1="post_outputv.txt";
	outputu1 = fopen(outputfilenameu1,"w+");
	outputv1 = fopen(outputfilenamev1,"w+");
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
			fprintf(outputu1,"%f ",uu[j][i]);
		}
		fprintf(outputu1,"\n");
	}
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
            fprintf(outputv1,"%f ",vv[j][i]);
		}
		fprintf(outputv1,"\n");
	}
	fclose(outputu1);
	fclose(outputv1);
	printf("uu into file:%s\n",outputfilenameu1);
	printf("vv into file:%s\n",outputfilenamev1);

	FREE_RMATRIX(u,0,imax+1,0,jmax+1);
	FREE_RMATRIX(v,0,imax+1,0,jmax+1);
	FREE_RMATRIX(p,0,imax+1,0,jmax+1);
	FREE_RMATRIX(f,0,imax+1,0,jmax+1);
	FREE_RMATRIX(g,0,imax+1,0,jmax+1);
	FREE_RMATRIX(rhs,0,imax+1,0,jmax+1);
	FREE_RMATRIX(uu,0,imax,0,jmax);
	FREE_RMATRIX(vv,0,imax,0,jmax);
	free(xx);
	free(yy);
    return 0;
}
