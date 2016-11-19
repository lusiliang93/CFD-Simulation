X=0:dx:xlength;
Y=0:dy:ylength;
[xx,yy]=meshgrid(X,Y);
uu=zeros(jmax+1,imax+1);
vv=zeros(jmax+1,imax+1);
for m=1:(imax+1)
    for n=1:(jmax+1)
        x=xx(n,m);
        y=yy(n,m);
        i=floor(x/dx)+1;
        j=floor((y+dy/2)/dy)+1;
        x1=(i-1)*dx;
        y1=((j-1)-0.5)*dy;
        x2=i*dx;
        y2=(j-1/2)*dy;
        u1=u(j+1-1,i+1-1);
        u2=u(j+1-1,i+1);
        u3=u(j+1,i+1-1);
        u4=u(j+1,i+1);
        uu(n,m)=1/(dx*dy)*((x2-x)*(y2-y)*u1+(x-x1)*(y2-y)*u2+(x2-x)*(y-y1)*u3+(x-x1)*(y-y1)*u4);
    end
end

for m=1:imax
    for n=1:jmax
        x=xx(n,m);
        y=yy(n,m);
        j=floor(y/dy)+1;
        i=floor((x+dx/2)/dx)+1;
        y1=(j-1)*dy;
        x1=((i-1)-0.5)*dx;
        y2=j*dy;
        x2=(i-1/2)*dx;
        v1=v(j+1-1,i+1-1);
        v2=v(j+1-1,i+1);
        v3=v(j+1,i+1-1);
        v4=v(j+1,i+1);
        vv(n,m)=1/(dx*dy)*((x2-x)*(y2-y)*v1+(x-x1)*(y2-y)*v2+(x2-x)*(y-y1)*v3+(x-x1)*(y-y1)*v4);
    end
end
axis([0,1.1,0,2.5]);
quiver(xx,yy,uu,vv,5);
axis([0,1.1,0,2.5]);
figure;
streamslice(xx,yy,uu,vv);
axis tight
