
%  ========================================================================
%      MFS for the 2D Laplace equation
%          ---- including circle, rectangle, irregular boundary
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
clear; clc; format short e;
t0 = clock; % Time Start   

addpath('Examples');
iexample=2;
style='Laplace'; 
iboundary=2;   % boundary geometry: 1-circle; 2-rectangle
D=2;           % distance of the real and virtual boundary

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
ni=50; 
theta=linspace(0,2*pi,ni);
xi=0.5*cos(theta);yi=0.5*sin(theta); xi=xi'; yi=yi';
% xi=linspace(-1,1,ni);yi=0.5*ones(1,ni); xi=xi'; yi=yi';
%xi=0.0*ones(1,ni);yi=linspace(-1,1,ni); xi=xi'; yi=yi';

% nni=20; ni=nni*nni;
% [xxi,yyi] = meshgrid(linspace(-1,1,nni),linspace(-1,1,nni));  
% xi = reshape(xxi,ni,1);
% yi = reshape(yyi,ni,1);

figure(1);h=fill(xp,yp,'k'); hold on; plot(xp,yp,'r.',xi,yi,'b.',xpv,ypv,'k.');
legend('Domain', 'Bounddary nodes','Test points','Virtual points');
xlabel('x');ylabel('y')
set(h,'edgealpha',0,'facealpha',0.3)

% ======== Exact solution & Boundary condition
[u_exact,ux_exact,uy_exact,f_source,lamda] = test_example(style, iexample);

b = u_exact(xp, yp);
b(1:nd) = u_exact(xp(1:nd), yp(1:nd));
b(nd+1:nb) = ux_exact(xp(nd+1:nb), yp(nd+1:nb)).*nx(nd+1:nb)+uy_exact(xp(nd+1:nb), yp(nd+1:nb)).*ny(nd+1:nb);

% ======== Matrix computation of MFS
r = sqrt((xp-xpv').^2+(yp-ypv').^2);
ri = sqrt((xi-xpv').^2+(yi-ypv').^2);
G1=-1/2/pi*log(r(1:nd,:));
G2=-1/2/pi*((xp(nd+1:nb)-xpv').*nx(nd+1:nb)+(yp(nd+1:nb)-ypv').*ny(nd+1:nb))./r(nd+1:nb,:).^2;
G=[G1;G2];
coeff=G\b;
u_MFS=-1/2/pi*log(ri)*coeff;% Calculate interior points

% ===== Calculate Errors 
err=abs(u_exact(xi, yi)-u_MFS);
error_max = max(abs(u_exact(xi, yi)-u_MFS)); % Maximum absolue errors
error_global = sqrt(sum((u_exact(xi, yi)-u_MFS).^2)/sum(u_exact(xi, yi).^2)); % Global error

% ===== Figures 
figure(2)
subplot(1,2,1)
plot(1:ni,u_exact(xi, yi),1:ni,u_MFS,'o');title('MFS solution');
legend('Exact','MFS');xlabel('number of nodes');ylabel('u')
subplot(1,2,2)
plot(1:ni,err,'-.');title('Absolute error');xlabel('number of nodes');ylabel('Absolute error')

%======================Output========
disp('---- MFS -  Laplace equation ----')
disp(['Boun. nodes No.: ',num2str(nb)])
disp(['Test nodes No.: ',num2str(ni)])
disp(['Elapsed total_time is: ',num2str(etime(clock,t0)),' seconds.'])
disp(['Max_Error: ',num2str(error_max)]);
disp(['Global_Error: ',num2str(error_global)]);
 