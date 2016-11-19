RHS=zeros(jmax+2,imax+2);
for i=1:imax
    for j=1:jmax
       RHS(j+1,i+1)=1/dt*((F(j+1,i+1)-F(j+1,i-1+1))/dx+(G(j+1,i+1)-G(j-1+1,i+1))/dy);
    end
end
