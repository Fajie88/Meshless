
%  ========================================================================
%      BKM for the 2D Helmholtz (modified) equation
%          ---- including circle or rectangle boundary (simple_geometry.m)
%          ---- including irregular boundary (BoundaryGeometry.m)
%          ---- Mixed boundary condition
%          ---- test_example (Examples\test_example.m)
% 
%      Author:  Fajie Wang          All Rights Reserved. 
%      Date: 2019-04-18
%  --------------------------------------------------------------------------
%  Qingdao University 
%  National Engineering Research Center for Intelligent Electrical Vehicle Power System  
%  E-mail: wfj1218@126.com
%  ========================================================================
clear; clc; format long;
t0 = clock; % Time Start   
addpath('Boundary_Geometry_data\2D');
addpath('Examples');

% ====== Input
iexample=1;
style='Helmholtz';   
%style='ModifiedHelmholtz';

% % % ===========circle/rectangle==================
% % % ---Boundary nodes and area of influence-------
% iboundary=1; nb=200;
% [xp,yp,n_x,n_y,l,xi,yi,ni] = simple_geometry(iboundary,nb);
% figure(1);plot(xp,yp,'r.',xi,yi,'b.');legend('Bounddary nodes','Test points');

% % % ===========irregular domain==================
% % % ---Boundary nodes and area of influence-------
iboundary=2; idomain=1;
[inp, bp, S, normal_vector, N_total,l]=BoundaryGeometry(iboundary, idomain);
xp=bp(:,1); yp=bp(:,2);
n_x=normal_vector(:,1); n_y=normal_vector(:,2);
nb=length(bp);ni=length(inp);
xi=inp(:,1); yi=inp(:,2);
figure(1);plot(xp,yp,'r.',xi,yi,'b.');legend('Bounddary nodes','Test points');

% ======= Numbers of Dirichlet and Nuemann boundary conditions
nd=round(0.5*nb);
nn=nb-nd;

% ======== Exact solution & Boundary condition
[u_exact,ux_exact,uy_exact,f_source,lamda] = test_example(style, iexample) ; %Helmholtz (Modified) 
b(1:nd) = u_exact(xp(1:nd), yp(1:nd));
b(nd+1:nb) = ux_exact(xp(nd+1:nb), yp(nd+1:nb)).*n_x(nd+1:nb)+uy_exact(xp(nd+1:nb), yp(nd+1:nb)).*n_y(nd+1:nb);
b=b';

% ======== Matrix computation of BKM
r = sqrt((xp-xp').^2+(yp-yp').^2);
ri = sqrt((xi-xp').^2+(yi-yp').^2);
if strcmp(style, 'Helmholtz')==1
G1 = besselj(0,lamda*r(1:nd,:));    % besselh is zero-order Hankel function of the first kind
G2 = lamda*besselj(-1,lamda*r(nd+1:nb,:)).*((xp(nd+1:nb)-xp').*n_x(nd+1:nb)+(yp(nd+1:nb)-yp').*n_y(nd+1:nb))./r(nd+1:nb,:); 
G2(isnan(G2))=0.0;   % besselj(-1,0)==0, (x-xj)/r==1, when r==0.    @@@@@@@
G = [G1; G2];
coeff = G\b;
u_BKM=besselj(0,lamda*ri)*coeff;% Calculate interior points
elseif strcmp(style, 'ModifiedHelmholtz')==1
G1 = besseli(0,lamda*r(1:nd,:));    % besselh is zero-order Hankel function of the first kind
G2 = 1i*lamda*besselj(-1,1i*lamda*r(nd+1:nb,:)).*((xp(nd+1:nb)-xp').*n_x(nd+1:nb)+(yp(nd+1:nb)-yp').*n_y(nd+1:nb))./r(nd+1:nb,:); 
G2(isnan(G2))=0.0;   % besselj(-1,0)==0, (x-xj)/r==1, when r==0.   @@@@@@@
G = [G1; G2];
coeff = G\b;
u_BKM=besseli(0,lamda*ri)*coeff;% Calculate interior points
end

t1 = clock; % Time Start 
% ===== Calculate Errors 
err=abs(u_exact(xi, yi)-u_BKM);
error_max = max(abs(u_exact(xi, yi)-u_BKM)); % Maximum absolue errors
error_global = sqrt(sum((u_exact(xi, yi)-u_BKM).^2)/sum(u_exact(xi, yi).^2)); % Global error

% ===== Figures 
figure(2)
subplot(1,2,1)
plot3(xi, yi, u_exact(xi, yi),'o'); hold on;
plot3(xi, yi, u_BKM,'*');
title('BKM solution');legend('Exact','BKM');xlabel('x');ylabel('y');zlabel('u')
subplot(1,2,2)
plot3(xi, yi, err,'.');title('Absolute error');xlabel('x');ylabel('y');zlabel('Absolute error')

%======================Output========
disp('---- BKM -  Helmholtz (modified) equation ----')
disp(['Eq: ',style])
disp(['Boun. nodes No.: ',num2str(nb)])
disp(['Test nodes No.: ',num2str(ni)])
disp(['Elapsed total_time is: ',num2str(etime(t1,t0)),' seconds.'])
disp(['Max_Error: ',num2str(error_max)]);
disp(['Global_Error: ',num2str(error_global)]);
 