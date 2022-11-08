
%  ========================================================================
%      MFS for the 2D Helmholtz /Modified Helmholtz equation
%          ---- including circle, rectangle, irregular boundary
%           ---- Mixed Boundary Conditions 
%           ---- test_example (Examples\test_example.m)
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
addpath('Examples');

style='Helmholtz';   % iexample=1,2,3,4,5
% style='ModifiedHelmholtz';  % iexample=1,2,3
iexample=5;
iboundary=2;   % boundary geometry: 1-circle; 2-rectangle
D=2;        % distance of the real and virtual boundary

% ====== Boundary nodes and source points
if iboundary==1   %---circle
nb=50;
r=2.0;thetab=linspace(0,2*pi-2*pi/nb,nb);
xp=r.*cos(thetab);yp=r.*sin(thetab);      % Boundary nodes
xpv=(r+D).*cos(thetab);ypv=(r+D).*sin(thetab);  % Source points
xp=xp';yp=yp'; xpv=xpv';ypv=ypv'; 
nx=cos(thetab); ny=sin(thetab);nx=nx'; ny=ny';  %Normal vector
elseif iboundary==2   %---irregular
    nb=50;rr=5;
    theta=linspace(0,2*pi-2*pi/nb,nb);
    n=8;
    r=(n^2+2*n+2-2*(n+1)*cos(n.*theta))./n^2;
    rt=(sin(n.*theta).*(2*n + 2))/n;
    xp=r.*cos(theta); yp=r.*sin(theta); % Boundary nodes
    xpv=(r+D).*cos(theta); ypv=(r+D).*sin(theta); % Source points
    nx=rt.*sin(theta)+r.*cos(theta);
    ny=-(rt.*cos(theta)-r.*sin(theta));
    xp=xp';yp=yp'; nx=nx'; ny=ny'; xpv=xpv';ypv=ypv'; 
%------------------边界节点所在影响域的弧长----------------------  
    ac=@(thea) sqrt((((sin(n.*thea).*(2*n + 2))/n).*cos(thea)-((n^2+2*n+2-2*(n+1)*cos(n.*thea))./n^2).*sin(thea)).^2+...
    (((sin(n.*thea).*(2*n + 2))/n).*sin(thea)+((n^2+2*n+2-2*(n+1)*cos(n.*thea))./n^2).*cos(thea)).^2);
% for k11=1:nb
%     l(k11)=quadl(ac,theta(k11)-2*pi/nb,theta(k11)+2*pi/nb)/(2*pi);
% end
       
elseif iboundary==3   %---rectangle
x_scale=[-1 1]; y_scale=[-1 1];
nx=22; ny=nx; nbx=nx-1;nby=ny-1; nb=2*(nbx+nby);
x=linspace(x_scale(1),x_scale(2),nx);
y=linspace(y_scale(1),y_scale(2),ny);
xp=[(x(1:end-1)+x(2:end))/2 x_scale(2)*ones(1,nby) fliplr((x(1:end-1)+x(2:end))/2) x_scale(1)*ones(1,nby)]; % bottom,right,top,left
yp=[y_scale(1)*ones(1,nbx) (y(1:end-1)+y(2:end))/2 y_scale(2)*ones(1,nbx) fliplr((y(1:end-1)+y(2:end))/2)];
nx=[zeros(nbx,1);ones(nby,1);zeros(nbx,1);-ones(nby,1)];  %Normal vector1
ny=[-ones(nbx,1);zeros(nby,1);ones(nbx,1);zeros(nby,1)];  %Normal vector2
% xpv=xp;ypv=yp;
% ypv(1:nbx)=yp(1:nbx)-D;xpv(nbx+1:nbx+nby)=xp(nbx+1:nbx+nby)+D; ypv(nbx+nby+1:2*nbx+nby)=yp(nbx+nby+1:2*nbx+nby)+D;xpv(2*nbx+nby+1:2*nbx+2*nby)=xp(2*nbx+nby+1:2*nbx+2*nby)-D; 
% xp=xp';yp=yp';xpv=xpv';ypv=ypv';

thetab=linspace(0,2*pi-2*pi/nb,nb); r=sqrt((x_scale(2)-x_scale(1))^2+(y_scale(2)-y_scale(1))^2);
xpv=(r+D)/2.*cos(thetab);ypv=(r+D)/2.*sin(thetab); 
xp=xp';yp=yp';xpv=xpv';ypv=ypv';
end

% ======= Numbers of Dirichlet and Nuemann boundary conditions
nd=round(0.5*nb);
nn=nb-nd;


% ======= Interior points (for testing)
ni=20;
theta=linspace(0,2*pi,ni);
xi=0.5*cos(theta);yi=0.5*sin(theta); xi=xi'; yi=yi';
figure(1);h=fill(xp,yp,'k'); hold on; plot(xp,yp,'r.',xi,yi,'b.',xpv,ypv,'k.');legend('Domain', 'Bounddary nodes','Test points','Virtual points');
set(h,'edgealpha',0,'facealpha',0.3)

% ======== Exact solution & Boundary condition
[u_exact,ux_exact,uy_exact,f_source,lamda] = test_example(style, iexample) ; %ModifiedHelmholtz
 
b(1:nd) = u_exact(xp(1:nd), yp(1:nd));
b(nd+1:nb) = ux_exact(xp(nd+1:nb), yp(nd+1:nb)).*nx(nd+1:nb)+uy_exact(xp(nd+1:nb), yp(nd+1:nb)).*ny(nd+1:nb);
b=b';

% ======== Matrix computation of MFS
H01x = @(x) -besselj(1,x)+1i*0.5*(bessely(-1,x)-bessely(1,x));  % Derivative of zero-order Hankel function of the first kind
K0x = @(x) pi/2*(besselj(1,1i*x)+1i*0.5*(bessely(1,1i*x)-(bessely(-1,1i*x))));  % Derivative of zero-order Hankel function of the first kind

r = sqrt((xp-xpv').^2+(yp-ypv').^2);
ri = sqrt((xi-xpv').^2+(yi-ypv').^2);
if strcmp(style, 'Helmholtz')
G1=1i/4*besselh(0,1,lamda*r(1:nd,:));    % besselh is zero-order Hankel function of the first kind
G2=1i/4*lamda*H01x(lamda*r(nd+1:nb,:)).*((xp(nd+1:nb)-xpv').*nx(nd+1:nb)+(yp(nd+1:nb)-ypv').*ny(nd+1:nb))./r(nd+1:nb,:); 
G=[G1;G2];
coeff=G\b;
u_MFS=1i/4*besselh(0,1,lamda*ri)*coeff; % Calculate interior points

elseif strcmp(style, 'ModifiedHelmholtz')
G1=1/(2*pi)*besselk(0,lamda*r(1:nd,:));  % besselk is the 0th order modified Bessel function of the second kind
G2=1/(2*pi)*lamda*K0x(lamda*r(nd+1:nb,:)).*((xp(nd+1:nb)-xpv').*nx(nd+1:nb)+(yp(nd+1:nb)-ypv').*ny(nd+1:nb))./r(nd+1:nb,:); 
G=[G1;G2];
coeff=G\b;
u_MFS=1/(2*pi)*besselk(0,lamda*ri)*coeff; % Calculate interior points
end

t1 = clock; % Time Start 
% ===== Calculate Errors 
err=abs(u_exact(xi, yi)-u_MFS);
error_max = max(abs(u_exact(xi, yi)-u_MFS)); % Maximum absolue errors
error_global = sqrt(sum((u_exact(xi, yi)-u_MFS).^2)/sum(u_exact(xi, yi).^2)); % Global error

% ===== Figures 
figure(2)
subplot(1,2,1)
plot(theta,u_exact(xi, yi),theta,u_MFS,'o');title('MFS solution');legend('Exact','MFS');xlabel('\theta');ylabel('u')
subplot(1,2,2)
plot(theta,err,'-.');title('Absolute error');xlabel('\theta');ylabel('Absolute error')

%======================Output========
disp('---- MFS -  Helmholtz /Modified Helmholtz equation ----')
% disp(['Total number of points: ',num2str(nb)])
disp(['Boun. nodes No.: ',num2str(nb)])
disp(['Test nodes No.: ',num2str(ni)])
disp(['Elapsed total_time is: ',num2str(etime(t1,t0)),' seconds.'])
disp(['Max_Error: ',num2str(error_max)]);
disp(['Global_Error: ',num2str(error_global)]);
 
