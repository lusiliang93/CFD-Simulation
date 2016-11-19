%~~~~~~~~~~~~~~~~~~~~~~~
% This is main module
%~~~~~~~~~~~~~~~~~~~~~~~
% Initialize
init_all();
tic;
% Enter main loop
while t<tend
    comp_delt();
    setbcond();
    %comp_temp();
    comp_fg();
    comp_rhs();
    p=poisson(RHS,imax,jmax,jI,iB,dx,dy,eps,itermax,omg,p);
    adap_uv();
    k=k+1;
    t=t+dt;
    %disp([t, round(k), dt]);
end
toc;
printf('The elapsed time is %f s\n',toc);
post



    