tend=5;t=0;k=1;               % End time, start time, # of main loop
imax=88;jmax=200;           % Mesh Size
xlength=1.1;ylength=2.5;    % Base domain
dx=xlength/imax;dy=ylength/jmax;
% Define domain
xB=0.3;xC=0.6;xD=0.8;xG=0.15;yI=1.2;
iB=round(xB/dx);iC=round(xC/dx);iD=round(xD/dx);iG=round(xG/dx);jI=round(yI/dy);
% computation control
tau=0.5;itermax=100;eps=0.001;omg=1.7;
gamma=0.9;Re=17000;gx=0;gy=0;
ui=0; vi=0; pi=0; 
Pr=0.7;beta=0.00034;
% Initialize U,V,P,T
ui=0;vi=0;pi=0;uin=-1;vin=1;
T_room=298;         % Room temperature 25
T_inlet=293;    % Inlet air temperature
T_svTop=318;    % Server top temperature
T_svRight=318;  % Server right temperature
init_uvpt();
% Initialize t,k
t0=0;k=1;
% Store solution
u_solution=zeros(202,90,1000);
v_solution=zeros(202,90,1000);
