% Compute Temperture for next timestep
Tnew=zeros(jmax+2,imax+2);
for i=2:imax+1
    for j=2:jmax+1
        if (i<iB+1)&&(j<jI+1)
            Tnew(j,i)=T_room;
        else
            AA=1/dx*(u(j,i)*(T(j,i)+T(j,i+1))/2-u(j,i-1)*(T(j,i-1)+T(j,i))/2)...
                +gamma/dx*(abs(u(j,i))*(T(j,i)-T(j,i+1))/2-abs(u(j,i-1))*(T(j,i-1)-T(j,i))/2);
            BB=1/dy*(v(j,i)*(T(j,i)+T(j,i+1))/2-v(j-1,i)*(T(j-1,i)+T(j,i))/2)...
                +gamma/dy*(abs(v(j,i))*(T(j,i)-T(j+1,i))/2-abs(v(j-1,i))*(T(j-1,i)-T(j,i))/2);
            CC=1/(Re*Pr)*((T(j,i+1)-2*T(j,i)+T(j,i-1))/(dx.^2)+(T(j+1,i)-2*T(j,i)+T(j-1,i))/(dy.^2));
            Tnew(j,i)=T(j,i)+dt*(CC-AA-BB);
        end
    end
end
T(2:jmax+1,2:imax+1)=Tnew(2:jmax+1,2:imax+1);
