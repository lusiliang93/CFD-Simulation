function p = poisson(RHS,imax,jmax,jI,iB,dx,dy,eps,itermax,omg,p)
for it=1:itermax
    % Copy pressure for boundary cells
    p(2:jI+1,1)=p(2:jI+1,2);                    % Server right, JB
    P(jI+2:jmax+1,1)=p(jI+2:jmax+1,2);          % Room left, HI
    P(2:jmax+1,imax+2)=p(2:jmax+1,imax+1);      % Room right, EF
    P(1,iB+2:imax+1)=P(2,iB+2:imax+1);          % Room bottom, BE
    P(jI+1,2:iB+1)=P(jI+2,2:iB+1);              % Server top, IJ
    P(jmax+1,2:imax+1)=P(jmax,2:imax+1);        % Room top, HF
    P(jI+1,iB+1)=(P(jI+2,iB+1)+P(jI+1,iB+2))/2; % Server NE corner, J
    % Used to compute residual
    rr=0;       % Sum of r(j,i)
    count=0;    % Count number of fluid cells
    for i=1:imax
        for j=1:jmax
            if (i<iB+1)&&(j<jI+1)
                p(j,i)=p(j,i);
            else
                % Compute P(n+1)
                eiw=1;eie=1;ejs=1;ejn=1;
                p(j+1,i+1)=(1-omg)*p(j+1,i+1)+...
                    omg/((eie+eiw)/(dx^2)+(ejn+ejs)/(dy^2))...
                    *((eie*p(j+1,i+1+1)+eiw*p(j+1,i-1+1))/(dx^2)+...
                    (ejn*p(j+1+1,i+1)+ejs*p(j-1+1,i+1))/(dy^2)-RHS(j+1,i+1));
                % Compute residual
                r(j+1,i+1)=(eie*(p(j+1,i+1+1)-p(j+1,i+1))-eiw*(p(j+1,i+1)-p(j+1,i-1+1)))/(dx^2)...
                +(ejn*(p(j+1+1,i+1)-p(j+1,i+1))-ejs*(p(j+1,i+1)-p(j-1+1,i+1)))/(dy^2)-RHS(j+1,i+1);
                count=count+1;
                rr=rr+r(j+1,i+1)^2;
            end
        end
    end
    Ls_norm=(rr/count)^0.5;
    if (Ls_norm<eps)
        disp('Converged...')
        disp(Ls_norm);
        break;
    end      
end
                
            