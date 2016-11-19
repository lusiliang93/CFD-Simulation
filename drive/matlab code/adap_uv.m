for i=1:imax-1
    for j=1:jmax
        if (i<iB+1)&&(j<jI+1)
            u(j,i)=u(j,i);
        else
            u(j+1,i+1)=F(j+1,i+1)-dt/dx*(p(j+1,i+1+1)-p(j+1,i+1));
        end
    end
end
for i=1:imax
    for j=1:jmax-1
        if (i<iB+1)&&(j<jI+1)
            v(j,i)=v(j,i);
        else
            v(j+1,i+1)=G(j+1,i+1)-dt/dy*(p(j+1+1,i+1)-p(j+1,i+1));
        end
    end
end
