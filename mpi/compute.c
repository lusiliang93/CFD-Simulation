/* Questions: 
Do I need to calculate delt in each processor or just broadcast?
Will the bcast of delt override sub_delta calculated in each processor? 
Whether child processor will wait bcast before calculate res?
Whether child porcessor will wait bcast before return in comp_delt? 
How to make the use of tag? 
Is that right to change calculation method for f,g,rhs,init,adap,setbound??*/
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mpi.h"
#include "init.h"
#include "boundary.h"
#include "uvp.h"


double max(double** u,int imax,int jmax){
	int i,j;
	double max=0;
	for(j=0;j<jmax+1;j++){
		for(i=0;i<imax+1;i++){
			if (u[j][i]>max)
				max=u[j][i];
		}
	}
    /*printf("max u:%f\n",max);*/
	return max;
}
// change  
void comp_delt(double* delt,int imax,int jmax,double delx,double dely,double **u,double **v,double Re,double tau,int procID,int nproc,int iproc,int jproc){
    double first,second,third,min;
    double maxu,maxv,submaxu,submaxv;
	double delta = 1/(delx*delx)+1/(dely*dely);
	//double *sub_delta;
	int i;
	MPI_Status status;
	first = Re/2/delta;
	min=first;
	maxu=max(u,imax,jmax);
	maxv=max(v,imax,jmax);
	if(procID==0){
		for(i=1;i<nproc;i++){
			MPI_Recv(&submaxu,1,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
			MPI_Recv(&submaxv,1,MPI_DOUBLE,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
			if(submaxu>maxu){
				maxu=submaxu;
			}
			if(submaxv>maxv){
				maxv=submaxv;
			}
		}
	    second = delx/abs(maxu);
	    third= dely/abs(maxv);
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
		MPI_Bcast(delt,1,MPI_DOUBLE,0,MPI_COMM_WORLD); /* is that right? Will it override delta in each processor?*/
	} 
	else {
		MPI_Send(&maxu,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		MPI_Send(&maxv,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	}

	return;
}

int compute(int procID,int nproc,char *inputname){
	double t=0,n=0;
	double xlength,ylength;
	double tend,tau,itermax,eps,omg,gamma;
	double Re,GX,GY,UI,VI,PI;
	int imax,jmax;
	int wW,wE,wS,wN;
	int iproc,jproc;
        char procid[10];
	//double x,y,x1,y1,x2,y2,u1,u2,u3,u4,v1,v2,v3,v4;
	double **u;
	double **v;
	double **p;
	double **f;
	double **g;
	double **rhs;
	//double **uu;
	//double **vv;
	//double *xx;
	//double *yy;
	double delx,dely,delt;
	int i,j;//,ii,jj;
	int rid,cid;
	int iw,ie,js,jn,tw,te,ts,tn;
    FILE *input;
    FILE *outputu;
    FILE *outputv;
   // FILE *outputu1;
   // FILE *outputv1;
    char outputfilenameu[128];
    char outputfilenamev[128];
   // char *outputfilenameu1;
   // char *outputfilenamev1;
    clock_t t1,t2;
    double  total_t;
 
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
	/* make sure nproc=iproc*jproc and imax,jmax must be divisible by iproc and jproc */
	fscanf(input,"%d %d\n",&iproc,&jproc);

	delx=xlength/imax;
	dely=ylength/jmax;
    printf("xlengh:%f ylength:%f jmax:%d imax:%d Re:%f UI:%f VI:%f PI:%f GX:%f GY:%f tend:%f tau:%f itermax:%f eps:%f omg:%f gamma:%f\n",xlength,ylength,jmax,imax,Re,UI,VI,PI,GX,GY,tend,tau,itermax,eps,omg,gamma);
    
    /* subdomain index */
    rid=procID/jproc; 
    cid=procID%jproc;
    /* neighboring processors ID */
    tw=(cid>0)? procID-1:-1;
    te=(cid<iproc-1)? procID+1:-1;
    ts=(rid>0)? procID-iproc:-1; //is that right?
    tn=(rid<jproc-1)? procID+iproc:-1;
    //if(tw==-1||te==-1||ts==-1||tn==-1){
    //	printf("Neighboring processor ID is out of bound!");
    //	return -1;
    //}
    /* allocate memory to specified subdomain */
    /* assign initial values to u,v,p,f,g,rhs,uu,vv*/
    /* why do I need iw,ie,js,jn?*/
    iw=imax/iproc*cid+1; /* correct? */
    ie=imax/iproc+iw-1;
    js=jmax/jproc*rid+1;
    jn=jmax/jproc+js-1;
	u=RMATRIX(js-1,jn+1,iw-2,ie+1);
	v=RMATRIX(0,jn+1-(js-2),iw-1,ie+1);
	p=RMATRIX(js-1,jn+1,iw-1,ie+1);
	f=RMATRIX(js-1,jn+1,iw-2,ie+1);
	g=RMATRIX(0,jn+1-(js-2),iw-1,ie+1);
	rhs=RMATRIX(js,jn,iw,ie);
	//uu=RMATRIX(js,jn,iw,ie);
	//vv=RMATRIX(js,jn,iw,ie);
	/* allocate memory to xx, yy*/
	//xx=(double *)malloc((imax/iproc+1)*sizeof(double));
	//yy=(double *)malloc((jmax/jproc+1)*sizeof(double));
	init_uvp(u,v,p,imax/iproc,jmax/jproc,UI,VI,PI);

    t1=clock();
	while(t<tend){
        if(n==0){
		/*comp_delt(&delt,imax,jmax,delx,dely,u,v,Re,tau); */
        delt=0.02;
		setbound(u,v,imax/iproc,jmax/jproc,wW, wE,wN,wS,rid,cid,iproc,jproc);
		comp_fg(u,v,f,g, imax/iproc,jmax/jproc,delt,delx,dely,GX,GY,gamma,Re,rid,cid,iproc,jproc);
		comp_rhs(f, g,rhs,imax/iproc,jmax/jproc,delt,delx,dely);
		/* synchronized */
		poisson(p,rhs,imax/iproc,jmax/jproc,delx,dely,eps,itermax,omg,iw,ie,js,jn,tw,te,ts,tn,rid,cid,iproc,jproc,procID,nproc);
		adap_uv(u,v,f,g,p,imax/iproc,jmax/jproc,delt,delx,dely,tw,te,ts,tn,nproc,rid,cid,iproc,jproc);
		t=t+delt;
		n++;
		printf("The current t:%f\n",t);
        }else{
            comp_delt(&delt,imax/iproc,jmax/jproc,delx,dely,u,v,Re,tau,procID,nproc,iproc,jproc);
            setbound(u,v,imax/iproc,jmax/jproc,wW,wE,wN,wS,rid,cid,iproc,jproc);
            comp_fg(u,v,f,g,imax/iproc,jmax/jproc,delt,delx,dely,GX,GY,gamma,Re,rid,cid,iproc,jproc);
            comp_rhs(f,g,rhs,imax/iproc,jmax/jproc,delt,delx,dely);
            /* synchronized */
            poisson(p,rhs,imax/iproc,jmax/jproc,delx,dely,eps,itermax,omg,iw,ie,js,jn,tw,te,ts,tn,rid,cid,iproc,jproc,procID,nproc);
            adap_uv(u,v,f,g,p,imax/iproc,jmax/jproc,delt,delx,dely,tw,te,ts,tn,nproc,rid,cid,iproc,jproc);
            t=t+delt;
            n++;
            /*printf("The current delt:%f\n",delt);*/
            printf("The current t:%f\n",t);
        }
	}
    t2=clock();
    total_t=(double)(t2-t1)/CLOCKS_PER_SEC;
    printf("Time elapsed:%f\n",total_t);

    /* change outputfilename for each processor!!(to do)*/
    //for(i=0;i<nproc;i++){
	// outputfilenameu="outputu.txt";
 //    outputfilenamev="outputv.txt";
        sprintf(procid,"%d",procID);
    	strcpy(outputfilenameu,"outputu_");
    	strcat(outputfilenameu,procid);
    	strcat(outputfilenameu,".txt");
    	strcpy(outputfilenamev,"outputv_");
    	strcat(outputfilenamev,procid);
    	strcat(outputfilenamev,".txt");


	outputu = fopen(outputfilenameu,"w+");
	outputv = fopen(outputfilenamev,"w+");
    /* may need to change back!!(to do)*/
	for(j=0;j<jmax/jproc+2;j++){
		for(i=0;i<imax/iproc+3;i++){
			fprintf(outputu,"%f ",u[j][i]);
		}
		fprintf(outputu,"\n");
	}
	for(j=0;j<jmax/jproc+3;j++){
		for(i=0;i<imax/iproc+2;i++){
            fprintf(outputv,"%f ",v[j][i]);
		}
		fprintf(outputv,"\n");
	}
	fclose(outputu);
	fclose(outputv);
	printf("u into file:%s\n",outputfilenameu);
	printf("v into file:%s\n",outputfilenamev);

	/** post for visualization*/
	/* change */
	// for(i=0;i<imax+1;i++){
	// 	xx[i]=delx*i;
	// }
	// for(j=0;j<jmax+1;j++){
	// 	yy[j]=dely*j;
	// }
	// for(j=0;j<jmax+1;j++){/*Is that right?*/
	// 	for(i=0;i<imax+1;i++){/*Is that right?*/
	// 		x=xx[i];/*Is that right?*/
	// 		y=yy[j];/*Is that right?*/
	// 		ii=floor(x/delx)+1;
	// 		jj=floor((y+dely/2)/dely)+1;
	// 		x1=(ii-1)*delx;
	// 		y1=((jj-1)-0.5)*dely;
	// 		x2=ii*delx;
	// 		y2=(jj-1/2)*dely;
	// 		u1=u[jj-1][ii-1+1];
	// 		u2=u[jj-1][ii+1];
	// 		u3=u[jj][ii-1+1];
	// 		u4=u[jj][ii+1];
	// 		uu[j][i]=1/(delx*dely)*((x2-x)*(y2-y)*u1+(x-x1)*(y2-y)*u2+(x2-x)*(y-y1)*u3+(x-x1)*(y-y1)*u4);
	// 	}
	// }
	// /* vv*/
	// for(j=0;j<jmax+1;j++){/*Is that right?*/
	// 	for(i=0;i<imax+1;i++){
	// 		x=xx[i];
	// 		y=yy[j];
	// 		jj=floor(y/dely)+1;
	// 		ii=floor((x+delx/2)/delx)+1;
	// 		y1=(jj-1)*dely;
	// 		x1=((ii-1)-0.5)*delx;
	// 		y2=jj*dely;
	// 		x2=(ii-0.5)*delx;
	// 		v1=v[jj-1+1][ii-1];
	// 		v2=v[jj-1+1][ii];
	// 		v3=v[jj+1][ii-1];
	// 		v4=v[jj+1][ii];
	// 		vv[j][i]=1/(delx*dely)*((x2-x)*(y2-y)*v1+(x-x1)*(y2-y)*v2+(x2-x)*(y-y1)*v3+(x-x1)*(y-y1)*v4);
	// 	}
	// }

	// outputfilenameu1="post_outputu.txt";
 //    outputfilenamev1="post_outputv.txt";
	// outputu1 = fopen(outputfilenameu1,"w+");
	// outputv1 = fopen(outputfilenamev1,"w+");
	// for(j=0;j<jmax+1;j++){
	// 	for(i=0;i<imax+1;i++){
	// 		fprintf(outputu1,"%f ",uu[j][i]);
	// 	}
	// 	fprintf(outputu1,"\n");
	// }
	// for(j=0;j<jmax+1;j++){
	// 	for(i=0;i<imax+1;i++){
 //            fprintf(outputv1,"%f ",vv[j][i]);
	// 	}
	// 	fprintf(outputv1,"\n");
	// }
	// fclose(outputu1);
	// fclose(outputv1);
	// printf("uu into file:%s\n",outputfilenameu1);
	// printf("vv into file:%s\n",outputfilenamev1);

	FREE_RMATRIX(u,js-1,jn+1,iw-2,ie+1);
	FREE_RMATRIX(v,js-2,jn+1,iw-1,ie+1);
	FREE_RMATRIX(p,js-1,jn+1,iw-1,ie+1);
	FREE_RMATRIX(f,js-1,jn+1,iw-2,ie+1);
	FREE_RMATRIX(g,js-2,jn+1,iw-1,ie+1);
	FREE_RMATRIX(rhs,js,jn,iw,ie);
	//FREE_RMATRIX(uu,js,jn,iw,ie);
	//FREE_RMATRIX(vv,js,jn,iw,ie);
	//free(xx);
	//free(yy);
	return 0;
}
