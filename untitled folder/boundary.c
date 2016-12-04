#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "boundary.h"
void setbound(double **u,double **v,int imax,int jmax,int wW, int wE,int wN,int wS,int rid,int cid){
	int j,i;
	int us=1;
	/* changed to fit the dimension */
	/* west side */
	if(cid==0){
		for(j=1;j<jmax+1;j++){ /* is that right?*/
			u[j][1]=0;		
			v[j+1][0]=-v[j+1][1]; 
		}
    }
    /* east side */
    if(cid==iproc-1){
    	for(j=1;j<jmax+1;j++){
    		u[j][imax+1]=0; /* is that right?*/
    		v[j+1][imax+1]=-v[j+1][imax]; 
    	}
    }
    /* south side */
    if(rid==0){
    	for(i=1;i<imax+1;i++){
    		v[1][i]=0;		
    		u[0][i+1]=-u[1][i+1];		
    	}
    }
	/* north side */
	if(rid==jproc-1){
		for(i=1;i<imax+1;i++){
			v[jmax+1][i]=0;
			u[jmax+1][i+1]=2*us-u[jmax][i+1];
		}
	}
	return;
}

void setspecbcond(double **u,double **v,int imax,int jmax,char **problem){
	return;
}
