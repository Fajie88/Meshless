
%  ========================================================================
%      SBM for the 2D Helmholtz (modified) equation with Drichlet boundary condition
%          ---- using empirical formula
%          ---- including circle or rectangle boundary 
% 
%          ---- Equation type
%               style=1: Helmholtz  
%               style=2: Modified Helmholtz
%
%          ---- Helmholtz
%          uii = -1/(2*pi)*(log(li/(2*pi))+log(lamda/2)+eulergam)+1i/4; 
%
%          ---- Modified Helmholtz
%          uii = -1/(2*pi)*(log(l(ID)/(2*pi))+log(lamda/2)+eulergam); 
%             
%          li: area of influence,
% 
%
%      Author:  Fajie Wang          All Rights Reserved. 
%      Date: 2019-04-25
%  --------------------------------------------------------------------------
%  Qingdao University 
%  National Engineering Research Center for Intelligent Electrical Vehicle Power System  
%  E-mail: wfj1218@126.com
%  ========================================================================
clear; clc; format long;
t0 = clock; % Time Start   
style=2;   % Equation type---> style=1: Helmholtz;  style=2: Modified Helmholtz
lamda=pi;  % Wave number
iboundary=2;   % boundary geometry: 1-circle; 2-rectangle

% ====== Boundary nodes and area of influence
if iboundary==1   %---circle
nb=800;
r=2.0;thetab=linspace(0,2*pi-2*pi/nb,nb);
xp=r.*cos(thetab);yp=r.*sin(thetab);    
l(1:nb)=r/nb; %area of influence
xp=xp';yp=yp';
elseif iboundary==2   %---rectangle
x_scale=[0 1]; y_scale=[0 1];
nx=200; ny=200; nbx=nx-1;nby=ny-1; nb=2*(nbx+nby);
x=linspace(x_scale(1),x_scale(2),nx);
y=linspace(y_scale(1),y_scale(2),ny);
xp=[(x(1:end-1)+x(2:end))/2 x_scale(2)*ones(1,nby) fliplr((x(1:end-1)+x(2:end))/2) x_scale(1)*ones(1,nby)]; % bottom,right,top,left
yp=[y_scale(1)*ones(1,nbx) (y(1:end-1)+y(2:end))/2 y_scale(2)*ones(1,nbx) fliplr((y(1:end-1)+y(2:end))/2)];
l=[x(2:end)-x(1:end-1) y(2:end)-y(1:end-1) x(2:end)-x(1:end-1) y(2:end)-y(1:end-1)]/(2*pi); %area of influence
xp=xp';yp=yp';
end

% ======= Interior points (for testing)
% XY=textread('innerpoints.txt');
% ni=length(XY);
% xi=XY(:,1);yi=XY(:,2);

ni=20;
theta=linspace(0,2*pi,ni);
xi=0.5+0.2*cos(theta);yi=0.5+0.2*sin(theta); xi=xi'; yi=yi';
figure(1);plot(xp(1:nb),yp(1:nb),'r.',xi,yi,'k.');legend('Bounddary nodes','Test points');

% ======== Exact solution & Boundary condition
if style==1
u_exact = @(x,y) cos(0.5*x+(lamda^2-0.5^2)^(1/2)*y);  % Helmholtz: Computers and Structures 83 (2005) 267¨C278: example 2
elseif style==2
u_exact = @(x,y) exp((lamda/2)*(x+3^(1/2)*y));       % Modified Helmholtz: Engineering Analysis with Boundary Elements 44 (2014) 112¨C119
end
b = u_exact(xp, yp);

% ======== Matrix computation of SBM
eulergam = 0.57721566490153286060651209;

if style==1
G=1i/4*besselh(0,1,lamda*sqrt((xp-xp').^2+(yp-yp').^2));      % besselh is zero-order Hankel function of the first kind
G(logical(eye(size(G))))=-1/(2*pi)*(log(l(:))+log(lamda/2)+eulergam)+1i/4;
coeff=G\b;
u_SBM=1i/4*besselh(0,1,lamda*sqrt((xi-xp').^2+(yi-yp').^2))*coeff; % Calculate interior points
elseif style==2
G=1/(2*pi)*besselk(0,lamda*sqrt((xp-xp').^2+(yp-yp').^2));  % besselk is the 0th order modified Bessel function of the second kind
G(logical(eye(size(G))))=-1.0/(2*pi)*(log(l(:))+log(lamda/2)+eulergam);   % OIF--Uii: Engineering Analysis with Boundary Elements 44 (2014) 112¨C119
% %G(logical(eye(size(G))))=-1/(2*pi)*(log(l(:))+log(1i*lamda/2)+eulergam)+1i/4;
coeff=G\b;
u_SBM=1/(2*pi)*besselk(0,lamda*sqrt((xi-xp').^2+(yi-yp').^2))*coeff; % Calculate interior points
end
t1 = clock; % Time Start 

% ===== Calculate Errors 
err=abs(u_exact(xi, yi)-u_SBM);
error_max = max(abs(u_exact(xi, yi)-u_SBM)); % Maximum absolue errors
error_global = sqrt(sum((u_exact(xi, yi)-u_SBM).^2)/sum(u_exact(xi, yi).^2)); % Global error

% ===== Figures 
figure(2)
subplot(1,2,1)
plot(yi,u_exact(xi, yi),yi,u_SBM,'o');title('SBM solution');legend('Exact','SBM');xlabel('\theta');ylabel('u')
subplot(1,2,2)
plot(yi,err,'-.');title('Absolute error');xlabel('\theta');ylabel('Absolute error')

% figure(2)
% subplot(1,2,1)
% plot(theta,u_exact(xi, yi),theta,u_SBM,'o');title('SBM solution');legend('Exact','SBM');xlabel('\theta');ylabel('u')
% subplot(1,2,2)
% plot(theta,err,'-.');title('Absolute error');xlabel('\theta');ylabel('Absolute error')

%======================Output========
disp('---- SBM -  Helmholtz equation ----')
disp(['Boun. nodes No.: ',num2str(nb)])
disp(['Test nodes No.: ',num2str(ni)])
disp(['Elapsed total_time is: ',num2str(etime(t1,t0)),' seconds.'])
disp(['Max_Error: ',num2str(error_max)]);
disp(['Global_Error: ',num2str(error_global)]);
 