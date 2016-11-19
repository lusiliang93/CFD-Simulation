%~~~~~~~~~~~~~~~~~~~~~~~~~~
% Set u,v boundary
%~~~~~~~~~~~~~~~~~~~~~~~~~~
% BC, no-slip
u(1,iB+2:iC+1)=-u(2,iB+2:iC+1);      
v(1,iB+2:iC+1)=0;
% CD, inlet
u(1,iC+2:iD+1)=2*uin-u(2,iC+2:iD+1);
v(1,iC+2:iD+1)=vin;
% DE, no-slip
u(1,iD+2:imax+1)=-u(2,iD+2:imax+1);
v(1,iD+2:imax+1)=0;
% EF, no-slip
u(2:jmax+1,imax+1)=0;
v(2:jmax+1,imax+2)=-v(2:jmax+1,imax+1);
% GF, no-slip
u(jmax+2,iG+2:imax+1)=-u(jmax+1,iG+2:imax+1);
v(jmax+1,iG+2:imax+1)=0;
% HG, outflow
u(jmax+2,2:iG+1)=u(jmax+1,2:iG+1);
v(jmax+1,2:iG+1)=v(jmax,2:iG+1);
% HI, symmerty
u(jI+2:jmax+1,1)=0;
v(jI+2:jmax+1,1)=v(jI+2:jmax+1,2);
% IJ, no-slip (server top)
u(jI+1,2:iB)=u(jI+2,2:iB);
v(jI+1,2:iB)=0;
% Corner J, no-slip on N,E
u(jI+1,iB+1)=0;
v(jI+1,iB+1)=0;
% JB, no-slip (server right)
u(2:jI,iB+1)=0;
v(2:jI,iB+1)=-v(2:jI,iB+2);
%~~~~~~~~~~~~~~~~~~~~~~~~~~
% Set temperature boundary
%~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
% BC, adiabatic
T(1,iB+2:iC+1)=2*T_room-T(2,iB+2:iC+1);      
% CD, inlet
T(1,iC+2:iD+1)=2*T_inlet-T(2,iC+2:iD+1);
% DE, adiabatic
T(1,iD+2:imax+1)=2*T_room-T(2,iD+2:imax+1);
% EF, adiabatic
T(2:jmax+1,imax+2)=2*T_room-T(2:jmax+1,imax+1);
% GF, adiabatic
T(jmax+1,iG+2:imax+1)=2*T_room-T(jmax+1,iG+2:imax+1);
% HG, outflow
T(jmax+1,2:iG+1)=T(jmax,2:iG+1);
% HI, symmerty
T(jI+2:jmax+1,1)=T(jI+2:jmax+1,2);
% IJ, server top
T(jI+1,2:iB)=2*T_svTop-T(jI+2,2:iB);
% JB, server right
T(2:jI,iB+1)=2*T_svRight-T(2:jI,iB+2);
% Corner J
T(jI+1,iB+1)=(T(jI+1,iB+2)+T(jI+2,iB+1))/2;
%}

