% Initialize F and G
F=zeros(jmax+2,imax+2); G=zeros(jmax+2,imax+2);
% Calculate F,G interior values
for i=1:imax-1
    for j=1:jmax
        if (i<iB+1)&&(j<jI+1)
            F(j,i)=F(j,i);  % Cell is in solid domain
        else
            a=u(j+1,i+1)+u(j+1,i+1+1);
            b=u(j+1,i-1+1)+u(j+1,i+1);
            c=u(j+1,i+1)-u(j+1,i+1+1);
            d=u(j+1,i-1+1)-u(j+1,i+1);
            e=u(j+1,i+1)+u(j+1+1,i+1);
            f=u(j-1+1,i+1)+u(j+1,i+1);
            g=u(j+1,i+1)-u(j+1+1,i+1);
            h=u(j-1+1,i+1)-u(j+1,i+1);
            va=v(j+1,i+1)+v(j+1,i+1+1);
            vb=v(j-1+1,i+1)+v(j-1+1,i+1+1);
            u2x=1/dx*((a/2)^2-(b/2)^2)...
                +gamma*1/dx*(abs(a)/2*c/2-abs(b)/2*d/2);
            uvy=1/dy*(va/2*e/2-vb/2*f/2)...
                +gamma*1/dy*(abs(va)/2*g/2-abs(vb)/2*h/2);
            u2x2=(u(j+1,i+1+1)-2*u(j+1,i+1)+u(j+1,i-1+1))/(dx^2);
            u2y2=(u(j+1+1,i+1)-2*u(j+1,i+1)+u(j-1+1,i+1))/(dy^2);
            F(j+1,i+1)=u(j+1,i+1)+dt*(1/Re*(u2x2+u2y2)-(u2x)-(uvy)+gx)...
                -beta*dt/2*(T(j,i)+T(j,i+1))*gx;
        end

    end
end

for i=1:imax
    for j=1:jmax-1
       if (i<iB+1)&&(j<jI+1)
           G(j,i)=G(j,i);  % Cell is in solid domain
       else
            a=v(j+1,i+1)+v(j+1,i+1+1);
            b=v(j+1,i-1+1)+v(j+1,i+1);
            c=v(j+1,i+1)-v(j+1,i+1+1);
            d=v(j+1,i-1+1)-v(j+1,i+1);
            e=v(j+1,i+1)+v(j+1+1,i+1);
            f=v(j-1+1,i+1)+v(j+1,i+1);
            g=v(j+1,i+1)-v(j+1+1,i+1);
            h=v(j-1+1,i+1)-v(j+1,i+1);
            ua=u(j+1,i+1)+u(j+1+1,i+1);
            ub=u(j+1,i-1+1)+u(j+1+1,i-1+1);
            uvx=1/dx*(ua/2*a/2-ub/2*b/2)...
                +gamma*1/dx*(abs(ua)/2*c/2-abs(ub)/2*d/2);
            v2y=1/dy*((e/2)^2-(f/2)^2)...
                +gamma*1/dy*(abs(e)/2*g/2-abs(f)/2*h/2);
            v2x2=(v(j+1,i+1+1)-2*v(j+1,i+1)+v(j+1,i-1+1))/(dx^2);
            v2y2=(v(j+1+1,i+1)-2*v(j+1,i+1)+v(j-1+1,i+1))/(dy^2);
            G(j+1,i+1)=v(j+1,i+1)+dt*(1/Re*(v2x2+v2y2)-uvx-v2y+gy)...
                -beta*dt/2*(T(j,i)+T(j+1,i))*gy;
       end
    end
end
% Update F and G at boundaries
F(2:jmax+1,imax+1)=u(2:jmax+1,imax+1);  % Room right, EF 
F(2:jI+1,iB+1)=u(2:jI+1,iB+1);          % Server right, JB
F(jI+2:jmax+1,1)=u(jI+2:jmax+1,1);      % Room left up, HI
G(1,iB+2:imax+1)=v(1,iB+2:imax+1);      % Room bottom, BE
G(jI+1,2:iB+1)=v(jI+1,2:iB+1);          % Server top, IJ
G(jmax+1,2:imax+1)=v(jmax+1,2:imax+1);  % Room top

       

       
       
        