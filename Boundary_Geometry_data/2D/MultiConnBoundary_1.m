function [xb,yb,xi,yi,nb,ni,pn1,pn2,lii]=MultiConnBoundary_1(nb1,nb2,irand,delta,x0,x1,y0,y1,nx,ny)
   nbtest=500;
   nb=nb1+nb2;
   n=9;
%Boundary points
    theta1=linspace(0,2*pi-2*pi/(nb1-1),nb1);
    theta2=linspace(0,2*pi-2*pi/(nb2-1),nb2);

%------------------边界节点所在影响域的弧长----------------------  
ac1=@(thea) sqrt((((n.*cos(n.*thea))./2).*cos(thea)-(2+0.5.*sin(n.*thea)).*sin(thea)).^2+...
    (((n.*cos(n.*thea))./2).*sin(thea)+(2+0.5.*sin(n.*thea)).*cos(thea)).^2);
for kk=1:nb1
    l1(kk)=quadl(ac1,theta1(kk)-2*pi/nb1,theta1(kk)+2*pi/nb1)/(2*pi);
end
ac2=@(thea) sqrt(((-8*cos(4*thea).*sin(4*thea)).*cos(thea)-(1+cos(4*thea).^2).*sin(thea)).^2+...
    ((-8*cos(4*thea).*sin(4*thea)).*sin(thea)+(1+cos(4*thea).^2).*cos(thea)).^2);
for kk=1:nb2
    l2(kk)=quadl(ac2,theta2(kk)-2*pi/nb2,theta2(kk)+2*pi/nb2)/(2*pi);
end


    r1=2+0.5.*sin(n.*theta1);
    rt1=(n.*cos(n.*theta1))./2;
    xb1=r1.*cos(theta1+0.5.*sin(n.*theta1));
    yb1=r1.*sin(theta1+0.5.*sin(n.*theta1));  
    pn11=rt1.*sin(theta1+0.5.*sin(n.*theta1))+r1.*(cos(theta1 + ...
        sin(n.*theta1)./2).*((n.*cos(n.*theta1))./2 + 1));
    pn12=-(rt1.*cos(theta1+0.5.*sin(n.*theta1))+r1.*(-sin(theta1 + ...
        sin(n.*theta1)./2).*((n.*cos(n.*theta1))./2 + 1)));
    
    
    r2=0.3+cos(4*theta2).^2;
    rt2=-8*cos(4*theta2).*sin(4*theta2);
    xb2=r2.*cos(theta2);
    yb2=r2.*sin(theta2);
    pn21=-(rt2.*sin(theta2)+r2.*cos(theta2));
    pn22=(rt2.*cos(theta2)-r2.*sin(theta2));
    
    
    xb(1:nb1)=xb1(1:nb1);
    yb(1:nb1)=yb1(1:nb1);
    pn1(1:nb1)=pn11(1:nb1);
    pn2(1:nb1)=pn12(1:nb1);
    xb(nb1+1:nb1+nb2)=xb2(1:nb2);
    yb(nb1+1:nb1+nb2)=yb2(1:nb2);
    pn1(nb1+1:nb1+nb2)=pn21(1:nb2);
    pn2(nb1+1:nb1+nb2)=pn22(1:nb2);
    
    lii(1:nb1)=l1(1:nb1);
    lii(nb1+1:nb1+nb2)=l2(1:nb2); 
    
pn1=pn1./sqrt(pn1.^2+pn2.^2);
pn2=pn2./sqrt(pn1.^2+pn2.^2);
pn1=pn1';pn2=pn2';

% Inner points
    theta1=linspace(0,2*pi,nbtest);
    theta2=linspace(0,2*pi,nbtest);
    
    r1=2+0.5.*sin(n.*theta1); 
    xx1=r1.*cos(theta1+0.5.*sin(n.*theta1));
    yy1=r1.*sin(theta1+0.5.*sin(n.*theta1));  
    
    r2=0.3+cos(4*theta2).^2;
    xx2=r2.*cos(theta2);
    yy2=r2.*sin(theta2);

    x(1:nbtest)=xx1(1:nbtest);
    y(1:nbtest)=yy1(1:nbtest);
    x(nbtest+1:2*nbtest)=xx2(1:nbtest);
    y(nbtest+1:2*nbtest)=yy2(1:nbtest);


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

k=0;
for j=1:NN
    if in1(j)==1 && in2(j)==0
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
fill(xx1,yy1,'m')
fill(xx2,yy2,'w')
% fill(xx3,yy3,'w')
% fill(xx4,yy4,'w')
end