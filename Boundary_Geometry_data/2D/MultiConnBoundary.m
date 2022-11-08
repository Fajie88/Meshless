function [xb,yb,xi,yi,nb,ni,pn1,pn2]=MultiConnBoundary(nb1,nb2,nb3,nb4,irand,delta,x0,x1,y0,y1,nx,ny)
   nbtest=200;
   nb=nb1+nb2+nb3+nb4;
%Boundary points
    theta1=linspace(0,2*pi-2*pi/(nb1-1),nb1);
    theta2=linspace(0,2*pi-2*pi/(nb2-1),nb2);
    theta3=linspace(0,2*pi-2*pi/(nb3-1),nb3);
    theta4=linspace(0,2*pi-2*pi/(nb4-1),nb4);
    
    r1=4*sqrt(cos(2.*theta1)+sqrt(1.1-sin(2.*theta1).*sin(2.*theta1)));
    rt1=-(2.*(2.*sin(2.*theta1) + (2.*cos(2.*theta1).*sin(2.*theta1))./(11/10 -...
        sin(2.*theta1).^2).^(1/2)))./(cos(2.*theta1) + (11/10 - sin(2.*theta1).^2).^(1/2)).^(1/2);
    xb1=r1.*cos(theta1);
    yb1=r1.*sin(theta1);
    pn11=rt1.*sin(theta1)+r1.*cos(theta1);
    pn12=-(rt1.*cos(theta1)-r1.*sin(theta1));
    
    
    r2=0.5*(cos(3.*theta2)+sqrt(2-sin(3.*theta2).*sin(3.*theta2))).^(1/3);
    rt2=-0.5*(3.*sin(3.*theta2) + (3.*cos(3.*theta2).*sin(3.*theta2))./(2 - ...
        sin(3.*theta2).^2).^(1/2))./(3.*(cos(3.*theta2) + (2 - sin(3.*theta2).^2).^(1/2)).^(2/3));
    xb2=r2.*cos(theta2);
    yb2=r2.*sin(theta2);
    pn21=-(rt2.*sin(theta2)+r2.*cos(theta2));
    pn22=(rt2.*cos(theta2)-r2.*sin(theta2));
    
    r3=0.5+cos(4*theta3).^2;
    rt3=-8*cos(4*theta3).*sin(4*theta3);
    xb3=-3+r3.*cos(theta3);
    yb3=r3.*sin(theta3); 
    pn31=-(rt3.*sin(theta3)+r3.*cos(theta3));
    pn32=(rt3.*cos(theta3)-r3.*sin(theta3));
    
    r4=0.7*(exp(sin(theta4)).*sin(2.*theta4).*sin(2.*theta4)+exp(cos(theta4)).*cos(2.*theta4).*cos(2.*theta4));
    rt4=0.7*(4.*cos(2.*theta4).*sin(2.*theta4).*exp(sin(theta4)) - cos(2.*theta4).^2.*exp(cos(theta4)).*sin(theta4) - ...
    4.*cos(2.*theta4).*sin(2.*theta4).*exp(cos(theta4)) + sin(2.*theta4).^2.*exp(sin(theta4)).*cos(theta4));
    xb4=3+r4.*cos(theta4);
    yb4=r4.*sin(theta4);
    pn41=-(rt4.*sin(theta4)+r4.*cos(theta4));
    pn42=(rt4.*cos(theta4)-r4.*sin(theta4));
    
    xb(1:nb1)=xb1(1:nb1);
    yb(1:nb1)=yb1(1:nb1);
    pn1(1:nb1)=pn11(1:nb1);
    pn2(1:nb1)=pn12(1:nb1);
    xb(nb1+1:nb1+nb2)=xb2(1:nb2);
    yb(nb1+1:nb1+nb2)=yb2(1:nb2);
    pn1(nb1+1:nb1+nb2)=pn21(1:nb2);
    pn2(nb1+1:nb1+nb2)=pn22(1:nb2);
    xb(nb1+nb2+1:nb1+nb2+nb3)=xb3(1:nb3);
    yb(nb1+nb2+1:nb1+nb2+nb3)=yb3(1:nb3);
    pn1(nb1+nb2+1:nb1+nb2+nb3)=pn31(1:nb3);
    pn2(nb1+nb2+1:nb1+nb2+nb3)=pn32(1:nb3);
    xb(nb1+nb2+nb3+1:nb1+nb2+nb3+nb4)=xb4(1:nb4);
    yb(nb1+nb2+nb3+1:nb1+nb2+nb3+nb4)=yb4(1:nb4);
    pn1(nb1+nb2+nb3+1:nb1+nb2+nb3+nb4)=pn41(1:nb4);
    pn2(nb1+nb2+nb3+1:nb1+nb2+nb3+nb4)=pn42(1:nb4);
    
    
pn1=pn1./sqrt(pn1.^2+pn2.^2);
pn2=pn2./sqrt(pn1.^2+pn2.^2);
pn1=pn1';pn2=pn2';

% Inner points
    theta1=linspace(0,2*pi,nbtest);
    theta2=linspace(0,2*pi,nbtest);
    theta3=linspace(0,2*pi,nbtest);
    theta4=linspace(0,2*pi,nbtest);
    
    r1=4*sqrt(cos(2.*theta1)+sqrt(1.1-sin(2.*theta1).*sin(2.*theta1)));
    xx1=r1.*cos(theta1);
    yy1=r1.*sin(theta1);
    
    
    r2=0.5*(cos(3.*theta2)+sqrt(2-sin(3.*theta2).*sin(3.*theta2))).^(1/3);
    xx2=r2.*cos(theta2);
    yy2=r2.*sin(theta2);

    r3=0.5+cos(4*theta3).^2;
    xx3=-3+r3.*cos(theta3);
    yy3=r3.*sin(theta3); 
    
    r4=0.7*(exp(sin(theta4)).*sin(2.*theta4).*sin(2.*theta4)+exp(cos(theta4)).*cos(2.*theta4).*cos(2.*theta4));
    xx4=3+r4.*cos(theta4);
    yy4=r4.*sin(theta4);
    
    x(1:nbtest)=xx1(1:nbtest);
    y(1:nbtest)=yy1(1:nbtest);
    x(nbtest+1:2*nbtest)=xx2(1:nbtest);
    y(nbtest+1:2*nbtest)=yy2(1:nbtest);
    x(2*nbtest+1:3*nbtest)=xx3(1:nbtest);
    y(2*nbtest+1:3*nbtest)=yy3(1:nbtest);
    x(3*nbtest+1:4*nbtest)=xx4(1:nbtest);
    y(3*nbtest+1:4*nbtest)=yy4(1:nbtest);

ti=linspace(x0,x1,nx);
tj=linspace(y0,y1,ny);
NN=size(ti,2)*size(tj,2);
dx=(x1-x0)/nx;
dy=(y1-y0)/ny;
[xt,yt] = meshgrid(ti,tj);
B=[xt(:)';yt(:)'];
xx=B(1,:);
yy=B(2,:);
if irand==1
xx=xx;
yy=yy;
elseif irand==2
    for i=1:NN
xx(i)=xx(i)+dx*delta*rand;
yy(i)=yy(i)+dy*delta*rand;
    end
elseif irand==3
    for i=1:NN
xx(i)=x0+(x1-x0)*rand;
yy(i)=y0+(y1-y0)*rand;
    end
end
       
in1=inpolygon(xx,yy,xx1,yy1);
in2=inpolygon(xx,yy,xx2,yy2);
in3=inpolygon(xx,yy,xx3,yy3);
in4=inpolygon(xx,yy,xx4,yy4);
k=0;
for j=1:NN
    if in1(j)==1 && in2(j)==0 && in3(j)==0 && in4(j)==0
        k=k+1;
        xi(k)=xx(j);
        yi(k)=yy(j);
    else
    end
end
ni=k;
% 
% xb=xb+6;
% yb=yb+3;
% xi=xi+6;
% yi=yi+3;
% xx1=xx1+6;
% xx2=xx2+6;
% xx3=xx3+6;
% yy1=yy1+3;
% yy2=yy2+3;
% yy3=yy3+3;

% figure(1);
% plot(xi,yi,'b.')%,xx(~in),yy(~in),'.b')
% hold on;
% plot(xb,yb,'r.')
% hold on;
plot(xx1,yy1,'k-')
hold on;
plot(xx2,yy2,'k-')
hold on;
plot(xx3,yy3,'k-')
hold on;
plot(xx4,yy4,'k-')
hold on;
% fill(xx1,yy1,'m')
% fill(xx2,yy2,'w')
% fill(xx3,yy3,'w')
% fill(xx4,yy4,'w')
end