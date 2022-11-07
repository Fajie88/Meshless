
%  ========================================================================
%      SBM for the 2D Laplace equation with Drichlet boundary condition
%                  ---- using empirical formula
%                  ---- including circle or rectangle boundary 
% 
%                  uii = -1/(2*pi)*log(li)/(2*pi)) 
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

iboundary=2;   % boundary geometry: 1-circle; 2-rectangle
style='Laplace';   % iexample=1,2,3,...,9
iexample=9;

% ====== Boundary nodes and area of influence
if iboundary==1   %---circle
nb=100;
r=2.0;thetab=linspace(0,2*pi-2*pi/nb,nb);
xp=r.*cos(thetab);yp=r.*sin(thetab);    
l(1:nb)=2*pi*r/nb; %area of influence
xp=xp';yp=yp';
elseif iboundary==2   %---rectangle
x_scale=[-1 1]; y_scale=[-1 1];
nx=102; ny=nx; nbx=nx-1;nby=ny-1; nb=2*(nbx+nby);
x=linspace(x_scale(1),x_scale(2),nx);
y=linspace(y_scale(1),y_scale(2),ny);
xp=[(x(1:end-1)+x(2:end))/2 x_scale(2)*ones(1,nby) fliplr((x(1:end-1)+x(2:end))/2) x_scale(1)*ones(1,nby)]; % bottom,right,top,left
yp=[y_scale(1)*ones(1,nbx) (y(1:end-1)+y(2:end))/2 y_scale(2)*ones(1,nbx) fliplr((y(1:end-1)+y(2:end))/2)];
l=[x(2:end)-x(1:end-1) y(2:end)-y(1:end-1) x(2:end)-x(1:end-1) y(2:end)-y(1:end-1)]; %area of influence
xp=xp';yp=yp';
end

% ======= Interior points (for testing)
% XY=textread('innerpoints.txt');
% ni=length(XY);
% xi=XY(:,1);yi=XY(:,2);
% 
% yi=linspace(0.0001,0.9999,ni); yi=yi';
% xi=0.5*ones(1,ni); xi=xi';

% yi=linspace(0.01,0.99,ni); yi=yi';
% xi=linspace(0.01,0.99,ni); xi=xi';

ni=40;
theta=linspace(0,2*pi,ni);
xi=0.5*cos(theta);yi=0.5*sin(theta); xi=xi'; yi=yi';
figure(1);plot(xp(1:nb),yp(1:nb),'r.',xi,yi,'k.');legend('Bounddary nodes','Test points');

% ======== Exact solution & Boundary condition
[u_exact,ux_exact,uy_exact,f_source,lamda] = test_example(style, iexample) ; % Laplace
b = u_exact(xp, yp);

% ======== Matrix computation of SBM
G=-1/2/pi*log(sqrt((xp-xp').^2+(yp-yp').^2));
G(logical(eye(size(G))))=-1/2/pi*log(l(:)/2/pi);
coeff=G\b;
u_SBM=-1/2/pi*log(sqrt((xi-xp').^2+(yi-yp').^2))*coeff;% Calculate interior points

% ===== Calculate Errors 
err=abs((u_exact(xi, yi)-u_SBM)./u_exact(xi, yi));
error_max = max(abs(u_exact(xi, yi)-u_SBM)); % Maximum absolue errors
error_global = sqrt(sum((u_exact(xi, yi)-u_SBM).^2)/sum(u_exact(xi, yi).^2)); % Global error

% ===== Figures 
% figure(2)
% subplot(1,2,1)
% plot(yi,u_exact(xi, yi),yi,u_SBM,'o');title('SBM solution');legend('Exact','SBM');xlabel('\theta');ylabel('u')
% subplot(1,2,2)
% plot(yi,err,'-.');title('Absolute error');xlabel('\theta');ylabel('Absolute error')

figure(2)
subplot(1,2,1)
plot(theta,u_exact(xi, yi),theta,u_SBM,'o');title('SBM solution');legend('Exact','SBM');xlabel('\theta');ylabel('u')
subplot(1,2,2)
plot(theta,err,'-.');title('Absolute error');xlabel('\theta');ylabel('Absolute error')
%======================Output========
disp('---- SBM -  Laplace equation ----')
% disp(['Total number of points: ',num2str(nb)])
disp(['Boun. nodes No.: ',num2str(nb)])
disp(['Test nodes No.: ',num2str(ni)])
disp(['Elapsed total_time is: ',num2str(etime(clock,t0)),' seconds.'])
disp(['Max_Error: ',num2str(error_max)]);
disp(['Global_Error: ',num2str(error_global)]);

fid1=fopen('sbm.txt','w');
fprintf(fid1,'%g\n', err');
fclose(fid1);
 