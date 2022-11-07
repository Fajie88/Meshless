
%  ========================================================================
%      SBM for the 2D Helmholtz (modified) equation with mixed boundary condition
%          ---- using empirical formula
%          ---- including circle or rectangle boundary 
%          ---- Equation type
%               style=1: Helmholtz  
%               style=2: Modified Helmholtz

%          ---- Helmholtz
%          uii = -1/(2*pi)*(log(li/(2*pi))+log(lamda/2)+eulergam)+1i/4; 
%          qii = -sum(G0(i~=j));   %  li: area of influence,
%                                     G0: fundamental solution of Laplace equation
%
%          ---- Modified Helmholtz
%          uii = -1/(2*pi)*(log(l(ID)/(2*pi))+log(lamda/2)+eulergam); 
%          qii = -sum(G0(i~=j));   %  li: area of influence,
%                                     G0: fundamental solution of Laplace equation                   
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

iboundary=2;   % boundary geometry: 1-circle; 2-rectangle
% style='Helmholtz';   % iexample=1,2,3,4
style='ModifiedHelmholtz';  % iexample=1,2
iexample=2;

% ====== Boundary nodes and area of influence
if iboundary==1   %---circle
nb=200;
r=2.0; thetab=linspace(0,2*pi-2*pi/nb,nb);
xp=r.*cos(thetab);yp=r.*sin(thetab); 
n_x=cos(thetab); n_y=sin(thetab); n_x=n_x'; n_y=n_y'; %Normal vector
l(1:nb)=2*pi*r/nb; %area of influence
xp=xp';yp=yp';
elseif iboundary==2   %---rectangle
x_scale=[-1 1]; y_scale=[-1 1];
nx=201; ny=nx; nbx=nx-1;nby=ny-1; nb=2*(nbx+nby);
x=linspace(x_scale(1),x_scale(2),nx);
y=linspace(y_scale(1),y_scale(2),ny);
xp=[(x(1:end-1)+x(2:end))/2 x_scale(2)*ones(1,nby) fliplr((x(1:end-1)+x(2:end))/2) x_scale(1)*ones(1,nby)]; % bottom,right,top,left
yp=[y_scale(1)*ones(1,nbx) (y(1:end-1)+y(2:end))/2 y_scale(2)*ones(1,nbx) fliplr((y(1:end-1)+y(2:end))/2)];
n_x=[zeros(nbx,1);ones(nby,1);zeros(nbx,1);-ones(nby,1)];  %Normal vector1
n_y=[-ones(nbx,1);zeros(nby,1);ones(nbx,1);zeros(nby,1)];  %Normal vector2
l=[x(2:end)-x(1:end-1) y(2:end)-y(1:end-1) x(2:end)-x(1:end-1) y(2:end)-y(1:end-1)]; %area of influence
xp=xp';yp=yp';
end
ID=1:floor(nb/2);
IN=setdiff((1:nb),ID);
% ======= Interior points (for testing)
% XY=textread('innerpoints.txt');
% ni=length(XY);
% xi=XY(:,1);yi=XY(:,2);

ni=100;
theta=linspace(0,2*pi,ni);
xi=0.5*cos(theta);yi=0.5*sin(theta); xi=xi'; yi=yi';
figure(1);plot(xp(1:nb),yp(1:nb),'r.',xi,yi,'k.');legend('Bounddary nodes','Test points');

% ======== Exact solution & Boundary condition
[u_exact,ux_exact,uy_exact,f_source,lamda] = test_example(style, iexample) ; %ModifiedHelmholtz


b1 = u_exact(xp(ID), yp(ID));
b2 = ux_exact(xp(IN), yp(IN)).*n_x(IN)+uy_exact(xp(IN), yp(IN)).*n_y(IN);
b=[b1; b2];

% ======== Matrix computation of SBM
eulergam = 0.57721566490153286060651209;

if strcmp(style, 'Helmholtz')
dr=sqrt((xp-xp').^2+(yp-yp').^2); rdotn = (xp(IN)-xp').*n_x(IN)+(yp(IN)-yp').*n_y(IN);
G1=1i/4*besselh(0,1,lamda*dr(ID,:));      % besselh is zero-order Hankel function of the first kind
G2=-1i/4*lamda*besselh(1,1,lamda*dr(IN,:)).*rdotn./dr(IN,:); 
G=[G1; G2];
%-----uii & qii---------
[~,E0] = GE0(dr,xp,yp,xp',yp',n_x',n_y');
E0(logical(eye(nb)))=0; 
Asum = sum(E0(IN,:),2);
uii=-1/(2*pi)*(log(l(ID)/(2*pi))+log(lamda/2)+eulergam)+1i/4;
% qii=0.5./l(IN);
qii=-Asum';
%-----------------------
GG=[uii qii];
G(logical(eye(size(G))))=GG; 
coeff=G\b;
u_SBM=1i/4*besselh(0,1,lamda*sqrt((xi-xp').^2+(yi-yp').^2))*coeff; % Calculate interior points
elseif strcmp(style, 'ModifiedHelmholtz')
dr=sqrt((xp-xp').^2+(yp-yp').^2); rdotn = (xp(IN)-xp').*n_x(IN)+(yp(IN)-yp').*n_y(IN);
G1=1/(2*pi)*besselk(0,lamda*dr(ID,:));      % besselk is the 0th order modified Bessel function of the second kind
G2=-lamda/(2*pi)*besselk(1,lamda*dr(IN,:)).*rdotn./dr(IN,:); 
G=[G1; G2];
%---------uii & qii-------------
[~,E0] = GE0(dr,xp,yp,xp',yp',n_x',n_y');
E0(logical(eye(nb)))=0; 
Asum = sum(E0(IN,:),2);
uii=-1/(2*pi)*(log(l(ID)/(2*pi))+log(lamda/2)+eulergam);
qii=-Asum';
%---------------------------------
GG=[uii qii];
G(logical(eye(size(G))))=GG; 
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
plot(theta,u_exact(xi, yi),theta,u_SBM,'o');title('SBM solution');legend('Exact','SBM');xlabel('\theta');ylabel('u')
subplot(1,2,2)
plot(theta,err,'-.');title('Absolute error');xlabel('\theta');ylabel('Absolute error')

% figure(2)
% subplot(1,2,1)
% plot(theta,u_exact(xi, yi),theta,u_SBM,'o');title('SBM solution');legend('Exact','SBM');xlabel('\theta');ylabel('u')
% subplot(1,2,2)
% plot(theta,err,'-.');title('Absolute error');xlabel('\theta');ylabel('Absolute error')

%======================Output========
disp('---- SBM -  Helmholtz / modified equation ----')
% disp(['Total number of points: ',num2str(nb)])
disp(['Boun. nodes No.: ',num2str(nb)])
disp(['Test nodes No.: ',num2str(ni)])
disp(['Elapsed total_time is: ',num2str(etime(t1,t0)),' seconds.'])
disp(['Max_Error: ',num2str(error_max)]);
disp(['Global_Error: ',num2str(error_global)]);
 