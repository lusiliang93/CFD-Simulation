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

    FREE_RMATRIX(uu);
    FREE_RMATRIX(vv);
    free(xx);
    free(yy);