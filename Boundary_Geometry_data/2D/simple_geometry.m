function [xp,yp,n_x,n_y,l,xi,yi,ni] = simple_geometry(iboundary,nb)

%including circle or rectangle boundary 

if iboundary==1   %---circle
nb=200;
r=2.0; thetab=linspace(0,2*pi-2*pi/nb,nb);
xp=r.*cos(thetab);yp=r.*sin(thetab); 
n_x=cos(thetab); n_y=sin(thetab); n_x=n_x'; n_y=n_y'; %Normal vector
l(1:nb)=2*pi*r/nb; %area of influence
xp=xp';yp=yp';
elseif iboundary==2   %---rectangle
x_scale=[0 1]; y_scale=[0 1];
nx=nb/4+1; ny=nx; nbx=nx-1;nby=ny-1; %nb=2*(nbx+nby);
x=linspace(x_scale(1),x_scale(2),nx);
y=linspace(y_scale(1),y_scale(2),ny);
xp=[(x(1:end-1)+x(2:end))/2 x_scale(2)*ones(1,nby) fliplr((x(1:end-1)+x(2:end))/2) x_scale(1)*ones(1,nby)]; % bottom,right,top,left
yp=[y_scale(1)*ones(1,nbx) (y(1:end-1)+y(2:end))/2 y_scale(2)*ones(1,nbx) fliplr((y(1:end-1)+y(2:end))/2)];
n_x=[zeros(nbx,1);ones(nby,1);zeros(nbx,1);-ones(nby,1)];  %Normal vector1
n_y=[-ones(nbx,1);zeros(nby,1);ones(nbx,1);zeros(nby,1)];  %Normal vector2
l=[x(2:end)-x(1:end-1) y(2:end)-y(1:end-1) x(2:end)-x(1:end-1) y(2:end)-y(1:end-1)]; %area of influence
xp=xp';yp=yp';
end
% ======= Interior points (for testing)
ni=50;
theta=linspace(0,2*pi,ni);
xi=0.5+0.1*cos(theta);yi=0.5+0.1*sin(theta); 
xi=xi'; yi=yi';
end

