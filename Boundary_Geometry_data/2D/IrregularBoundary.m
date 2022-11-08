function [xb,yb,xi,yi,ni,pn1,pn2,lii]=IrregularBoundary(nb,idomain,irand,delta,x0,x1,y0,y1,nx,ny)
nbtest=200;

%Boundary points
theta=linspace(0,2*pi-2*pi/nb,nb);
if idomain==1
    r=1;
    xb=r.*cos(theta);
    yb=r.*sin(theta);
    pn1=cos(theta);
    pn2=sin(theta);
    
    lii=r/nb*ones(1,nb);
elseif idomain==2
    n=8;
    r=(n^2+2*n+2-2*(n+1)*cos(n.*theta))./n^2;
    rt=(sin(n.*theta).*(2*n + 2))/n;
    xb=r.*cos(theta);
    yb=r.*sin(theta);
    pn1=rt.*sin(theta)+r.*cos(theta);
    pn2=-(rt.*cos(theta)-r.*sin(theta));
%------------------边界节点所在影响域的弧长----------------------  
    ac=@(thea) sqrt((((sin(n.*thea).*(2*n + 2))/n).*cos(thea)-((n^2+2*n+2-2*(n+1)*cos(n.*thea))./n^2).*sin(thea)).^2+...
    (((sin(n.*thea).*(2*n + 2))/n).*sin(thea)+((n^2+2*n+2-2*(n+1)*cos(n.*thea))./n^2).*cos(thea)).^2);
for k11=1:nb
    lii(k11)=quadl(ac,theta(k11)-2*pi/nb,theta(k11)+2*pi/nb)/(2*pi);
end
%----------------------------------------------------------------

elseif idomain==3
    r=exp(sin(theta)).*sin(2.*theta).*sin(2.*theta)+exp(cos(theta)).*cos(2.*theta).*cos(2.*theta);
    rt=4.*cos(2.*theta).*sin(2.*theta).*exp(sin(theta)) - cos(2.*theta).^2.*exp(cos(theta)).*sin(theta) - ...
        4.*cos(2.*theta).*sin(2.*theta).*exp(cos(theta)) + sin(2.*theta).^2.*exp(sin(theta)).*cos(theta);
    xb=r.*cos(theta);
    yb=r.*sin(theta);
    pn1=rt.*sin(theta)+r.*cos(theta);
    pn2=-(rt.*cos(theta)-r.*sin(theta));   
%------------------边界节点所在影响域的弧长----------------------  
ac=@(thea) sqrt(((4.*cos(2.*thea).*sin(2.*thea).*exp(sin(thea)) - cos(2.*thea).^2.*exp(cos(thea)).*sin(thea) - ...
        4.*cos(2.*thea).*sin(2.*thea).*exp(cos(thea)) + sin(2.*thea).^2.*exp(sin(thea)).*cos(thea)).*cos(thea)-...
        (exp(sin(thea)).*sin(2.*thea).*sin(2.*thea)+exp(cos(thea)).*cos(2.*thea).*cos(2.*thea)).*sin(thea)).^2+...
        ((4.*cos(2.*thea).*sin(2.*thea).*exp(sin(thea)) - cos(2.*thea).^2.*exp(cos(thea)).*sin(thea) - ...
        4.*cos(2.*thea).*sin(2.*thea).*exp(cos(thea)) + sin(2.*thea).^2.*exp(sin(thea)).*cos(thea)).*sin(thea)+...
        (exp(sin(thea)).*sin(2.*thea).*sin(2.*thea)+exp(cos(thea)).*cos(2.*thea).*cos(2.*thea)).*cos(thea)).^2);
for k11=1:nb
    lii(k11)=quadl(ac,theta(k11)-2*pi/nb,theta(k11)+2*pi/nb)/(2*pi);
end
%----------------------------------------------------------------

elseif idomain==4
    r=sqrt(4-2.*cos(2.*theta));
    rt=(2^(1/2).*sin(2.*theta))./(2 - cos(2.*theta)).^(1/2);
    xb=r.*cos(theta);
    yb=r.*sin(theta);
    pn1=rt.*sin(theta)+r.*cos(theta);
    pn2=-(rt.*cos(theta)-r.*sin(theta));
    %------------------边界节点所在影响域的弧长----------------------  
    ac=@(thea) sqrt((((2^(1/2).*sin(2.*thea))./(2 - cos(2.*thea)).^(1/2)).*cos(thea)-...
        (sqrt(4-2.*cos(2.*thea))).*sin(thea)).^2+(((2^(1/2).*sin(2.*thea))./(2 -...
        cos(2.*thea)).^(1/2)).*sin(thea)+(sqrt(4-2.*cos(2.*thea))).*cos(thea)).^2);
for k11=1:nb
    lii(k11)=quadl(ac,theta(k11)-2*pi/nb,theta(k11)+2*pi/nb)/(2*pi);
end    
%----------------------------------------------------------------

elseif idomain==5
    r=4*sqrt(cos(2.*theta)+sqrt(1.1-sin(2.*theta).*sin(2.*theta)));
    rt=-(2.*(2.*sin(2.*theta) + (2.*cos(2.*theta).*sin(2.*theta))./(11/10 ...
        - sin(2.*theta).^2).^(1/2)))./(cos(2.*theta) + (11/10 - sin(2.*theta).^2).^(1/2)).^(1/2);
    xb=r.*cos(theta);
    yb=r.*sin(theta);
    pn1=rt.*sin(theta)+r.*cos(theta);
    pn2=-(rt.*cos(theta)-r.*sin(theta));
    %------------------边界节点所在影响域的弧长----------------------  
    ac=@(thea) sqrt(((-(2.*(2.*sin(2.*thea) + (2.*cos(2.*thea).*sin(2.*thea))./(11/10 ...
        - sin(2.*thea).^2).^(1/2)))./(cos(2.*thea) + (11/10 - sin(2.*thea).^2).^(1/2)).^(1/2)).*cos(thea)-...
        (4*sqrt(cos(2.*thea)+sqrt(1.1-sin(2.*thea).*sin(2.*thea)))).*sin(thea)).^2+((-(2.*(2.*sin(2.*thea) +...
        (2.*cos(2.*thea).*sin(2.*thea))./(11/10- sin(2.*thea).^2).^(1/2)))./(cos(2.*thea) +(11/10 -...
        sin(2.*thea).^2).^(1/2)).^(1/2)).*sin(thea)+(4*sqrt(cos(2.*thea)+sqrt(1.1-sin(2.*thea).*sin(2.*thea)))).*cos(thea)).^2);
for k11=1:nb
    lii(k11)=quadl(ac,theta(k11)-2*pi/nb,theta(k11)+2*pi/nb)/(2*pi);
end   
%----------------------------------------------------------------

elseif idomain==6
    r=(cos(3.*theta)+sqrt(2-sin(3.*theta).*sin(3.*theta))).^(1/3);
    rt=-(3.*sin(3.*theta) + (3.*cos(3.*theta).*sin(3.*theta))./(2 - ...
        sin(3.*theta).^2).^(1/2))./(3.*(cos(3.*theta) + (2 - sin(3.*theta).^2).^(1/2)).^(2/3));
    xb=r.*cos(theta);
    yb=r.*sin(theta);
    pn1=rt.*sin(theta)+r.*cos(theta);
    pn2=-(rt.*cos(theta)-r.*sin(theta));
    %------------------边界节点所在影响域的弧长----------------------  
    ac=@(thea) sqrt(((-(3.*sin(3.*thea) + (3.*cos(3.*thea).*sin(3.*thea))./(2 - ...
        sin(3.*thea).^2).^(1/2))./(3.*(cos(3.*thea) + (2 - sin(3.*thea).^2).^(1/2)).^(2/3))).*cos(thea)-...
        ((cos(3.*thea)+sqrt(2-sin(3.*thea).*sin(3.*thea))).^(1/3)).*sin(thea)).^2+((-(3.*sin(3.*thea) +...
        (3.*cos(3.*thea).*sin(3.*thea))./(2 - ...
        sin(3.*thea).^2).^(1/2))./(3.*(cos(3.*thea) + (2 - sin(3.*thea).^2).^(1/2)).^(2/3))).*sin(thea)+...
        ((cos(3.*thea)+sqrt(2-sin(3.*thea).*sin(3.*thea))).^(1/3)).*cos(thea)).^2);
for k11=1:nb
    lii(k11)=quadl(ac,theta(k11)-2*pi/nb,theta(k11)+2*pi/nb)/(2*pi);
    end
%----------------------------------------------------------------

elseif idomain==7
    r=1;
    yb=r.*(2.*cos(theta)-cos(2.*theta));
    xb=r.*(2.*sin(theta)-sin(2.*theta)); 
    pn1=r.*(2.*sin(2.*theta) - 2.*sin(theta));
    pn2=r.*(2.*cos(2.*theta) - 2.*cos(theta));
    %------------------边界节点所在影响域的弧长----------------------  
%     dyt=r.*(-2.*sin(theta)+sin(2.*theta));
%     dxt=r.*(2.*cos(theta)-cos(2.*theta)); 
    ac=@(thea) sqrt((r.*(-2.*sin(thea)+sin(2.*thea))).^2+(r.*(2.*cos(thea)-cos(2.*thea))).^2);
for k11=1:nb
    lii(k11)=quadl(ac,theta(k11)-2*pi/nb,theta(k11)+2*pi/nb)/(2*pi);
    end
%----------------------------------------------------------------

elseif idomain==8
    n=9;
    r=2+0.5.*sin(n.*theta);
    rt=(n.*cos(n.*theta))./2;
    xb=r.*cos(theta+0.5.*sin(n.*theta));
    yb=r.*sin(theta+0.5.*sin(n.*theta));  
    pn1=rt.*sin(theta+0.5.*sin(n.*theta))+r.*(cos(theta + ...
        sin(n.*theta)./2).*((n.*cos(n.*theta))./2 + 1));
    pn2=-(rt.*cos(theta+0.5.*sin(n.*theta))+r.*(-sin(theta + ...
        sin(n.*theta)./2).*((n.*cos(n.*theta))./2 + 1)));
    %------------------边界节点所在影响域的弧长----------------------  
ac=@(thea) sqrt((((n.*cos(n.*thea))./2).*cos(thea)-(2+0.5.*sin(n.*thea)).*sin(thea)).^2+...
    (((n.*cos(n.*thea))./2).*sin(thea)+(2+0.5.*sin(n.*thea)).*cos(thea)).^2);
for k11=1:nb
    lii(k11)=quadl(ac,theta(k11)-2*pi/nb,theta(k11)+2*pi/nb)/(2*pi);
    end
%----------------------------------------------------------------

elseif idomain==9
    a=4; b=1;
    r=sqrt((a+b)^2+1-2*(a+b)*cos(a*theta/b));
    rt=(a*sin((a*theta)/b)*(2*a + 2*b))./(2*b*((a + b)^2 - cos((a*theta)/b)*(2*a + 2*b) + 1).^(1/2));
    xb=r.*cos(theta);
    yb=r.*sin(theta);  
    pn1=rt.*sin(theta)+r.*cos(theta);
    pn2=-(rt.*cos(theta)-r.*sin(theta));
    %------------------边界节点所在影响域的弧长----------------------  
ac=@(thea) sqrt((((a*sin((a*thea)/b)*(2*a + 2*b))./(2*b*((a + b)^2 - ...
    cos((a*thea)/b)*(2*a + 2*b) + 1).^(1/2))).*cos(thea)-(sqrt((a+b)^2+1-2*(a+b)*cos(a*thea/b))).*sin(thea)).^2+...
    (((a*sin((a*thea)/b)*(2*a + 2*b))./(2*b*((a + b)^2 -...
    cos((a*thea)/b)*(2*a + 2*b) + 1).^(1/2))).*sin(thea)+(sqrt((a+b)^2+1-2*(a+b)*cos(a*thea/b))).*cos(thea)).^2);
for k11=1:nb
    lii(k11)=quadl(ac,theta(k11)-2*pi/nb,theta(k11)+2*pi/nb)/(2*pi);
    end
%----------------------------------------------------------------

elseif idomain==10    
    r=1+cos(4*theta).^2;
    rt=-8*cos(4*theta).*sin(4*theta);
    xb=r.*cos(theta);
    yb=r.*sin(theta); 
    pn1=(rt.*sin(theta)+r.*cos(theta));
    pn2=-(rt.*cos(theta)-r.*sin(theta));
    %------------------边界节点所在影响域的弧长----------------------  
ac=@(thea) sqrt(((-8*cos(4*thea).*sin(4*thea)).*cos(thea)-(1+cos(4*thea).^2).*sin(thea)).^2+...
    ((-8*cos(4*thea).*sin(4*thea)).*sin(thea)+(1+cos(4*thea).^2).*cos(thea)).^2);
for k11=1:nb
    lii(k11)=quadl(ac,theta(k11)-2*pi/nb,theta(k11)+2*pi/nb)/(2*pi);
    end
%----------------------------------------------------------------

elseif idomain==11    
    r=(3+cos(theta-pi/7).*sin(4*theta))./(5+sin(2*theta));
    rt=(3*cos(3*theta + pi/7) + 5*sin(5*theta + ...
        (5*pi)/14))./(2*(sin(2*theta) + 5)) - ...
        (2*cos(2*theta).*(sin(4*theta).*sin(theta + ...
        (5*pi)/14) + 3))./(sin(2*theta) + 5).^2;
    xb=r.*cos(theta);
    yb=r.*sin(theta); 
    pn1=(rt.*sin(theta)+r.*cos(theta));
    pn2=-(rt.*cos(theta)-r.*sin(theta));
    %------------------边界节点所在影响域的弧长----------------------  
ac=@(thea) sqrt((((3*cos(3*thea + pi/7) + 5*sin(5*thea + ...
        (5*pi)/14))./(2*(sin(2*thea) + 5)) - ...
        (2*cos(2*thea).*(sin(4*thea).*sin(thea + ...
        (5*pi)/14) + 3))./(sin(2*thea) + 5).^2).*cos(thea)-...
        ((3+cos(thea-pi/7).*sin(4*thea))./(5+sin(2*thea))).*sin(thea)).^2+...
        (((3*cos(3*thea + pi/7) + 5*sin(5*thea + ...
        (5*pi)/14))./(2*(sin(2*thea) + 5)) - ...
        (2*cos(2*thea).*(sin(4*thea).*sin(thea + ...
        (5*pi)/14) + 3))./(sin(2*thea) + 5).^2).*sin(thea)+...
        ((3+cos(thea-pi/7).*sin(4*thea))./(5+sin(2*thea))).*cos(thea)).^2);
for k11=1:nb
    lii(k11)=quadl(ac,theta(k11)-2*pi/nb,theta(k11)+2*pi/nb)/(2*pi);
    end
%----------------------------------------------------------------
elseif idomain==12
    nbb=nb/2;
    theta=linspace(0,2*pi-2*pi/nbb,nbb);
    n=9;
    r=2+0.5.*sin(n.*theta);
    rt=(n.*cos(n.*theta))./2;
    xb=r.*cos(theta+0.5.*sin(n.*theta));
    yb=r.*sin(theta+0.5.*sin(n.*theta));  
    pn1=rt.*sin(theta+0.5.*sin(n.*theta))+r.*(cos(theta + ...
        sin(n.*theta)./2).*((n.*cos(n.*theta))./2 + 1));
    pn2=-(rt.*cos(theta+0.5.*sin(n.*theta))+r.*(-sin(theta + ...
        sin(n.*theta)./2).*((n.*cos(n.*theta))./2 + 1)));
    %------------------边界节点所在影响域的弧长----------------------  
ac=@(thea) sqrt((((n.*cos(n.*thea))./2).*cos(thea)-(2+0.5.*sin(n.*thea)).*sin(thea)).^2+...
    (((n.*cos(n.*thea))./2).*sin(thea)+(2+0.5.*sin(n.*thea)).*cos(thea)).^2);
for k11=1:nbb
    lii(k11)=quadl(ac,theta(k11)-2*pi/nbb,theta(k11)+2*pi/nbb)/(2*pi);
end

    r=1+cos(4*theta).^2;
    rt=-8*cos(4*theta).*sin(4*theta);
    xb(nbb+1:nb)=r.*cos(theta);
    yb(nbb+1:nb)=r.*sin(theta); 
    pn1(nbb+1:nb)=-(rt.*sin(theta)+r.*cos(theta));
    pn2(nbb+1:nb)=(rt.*cos(theta)-r.*sin(theta));
    %------------------边界节点所在影响域的弧长----------------------  
ac=@(thea) sqrt(((-8*cos(4*thea).*sin(4*thea)).*cos(thea)-(1+cos(4*thea).^2).*sin(thea)).^2+...
    ((-8*cos(4*thea).*sin(4*thea)).*sin(thea)+(1+cos(4*thea).^2).*cos(thea)).^2);
for k11=1:nbb
    lii(k11+nbb)=quadl(ac,theta(k11)-2*pi/nbb,theta(k11)+2*pi/nbb)/(2*pi);
end
    
    
end
pn1=pn1./sqrt(pn1.^2+pn2.^2);
pn2=pn2./sqrt(pn1.^2+pn2.^2);
pn1=pn1';pn2=pn2';

% Inner points
theta=linspace(0,2*pi,nbtest); 
if idomain==1
    r=1;
    x=r.*cos(theta);
    y=r.*sin(theta);
elseif idomain==2
    n=8;
    r=(n^2+2*n+2-2*(n+1)*cos(n.*theta))./n^2;
    x=r.*cos(theta);
    y=r.*sin(theta);
elseif idomain==3
    r=exp(sin(theta)).*sin(2.*theta).*sin(2.*theta)+exp(cos(theta)).*cos(2.*theta).*cos(2.*theta);
    x=r.*cos(theta);
    y=r.*sin(theta);
elseif idomain==4
    r=sqrt(4-2.*cos(2.*theta));
    x=r.*cos(theta);
    y=r.*sin(theta);
elseif idomain==5
    r=4*sqrt(cos(2.*theta)+sqrt(1.1-sin(2.*theta).*sin(2.*theta)));
    x=r.*cos(theta);
    y=r.*sin(theta);
elseif idomain==6
    r=(cos(3.*theta)+sqrt(2-sin(3.*theta).*sin(3.*theta))).^(1/3);
    x=r.*cos(theta);
    y=r.*sin(theta);
elseif idomain==7
    r=1;
    y=r.*(2.*cos(theta)-cos(2.*theta));
    x=r.*(2.*sin(theta)-sin(2.*theta));    
elseif idomain==8
    n=9;
    r=2+0.5.*sin(n.*theta); 
    x=r.*cos(theta+0.5.*sin(n.*theta));
    y=r.*sin(theta+0.5.*sin(n.*theta));  
elseif idomain==9
    a=4; b=1;
    r=sqrt((a+b)^2+1-2*(a+b)*cos(a*theta/b));
    x=r.*cos(theta);
    y=r.*sin(theta);   
elseif idomain==10    
    r=1+cos(4*theta).^2;
    x=r.*cos(theta);
    y=r.*sin(theta); 
 elseif idomain==11    
    r=(3+cos(theta-pi/7).*sin(4*theta))./(5+sin(2*theta));
    x=r.*cos(theta);
    y=r.*sin(theta); 
 elseif idomain==12  
    nbb=nb/2;
    theta=linspace(0,2*pi-2*pi/nbb,nbb);
    n=9;
    r=2+0.5.*sin(n.*theta);
    x=r.*cos(theta+0.5.*sin(n.*theta));
    y=r.*sin(theta+0.5.*sin(n.*theta));  

    r=1+cos(4*theta).^2;
    x(nbb+1:nb)=r.*cos(theta);
    y(nbb+1:nb)=r.*sin(theta); 
end


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
xx=x0+(x1-x0)*rand(nx,1);
yy=y0+(y1-y0)*rand(ny,1);
end
       
in=inpolygon(xx,yy,x,y);

k=0;
for j=1:NN
    if in(j)==1
        k=k+1;
        xi(k)=xx(j);
        yi(k)=yy(j);
    else
    end
end
ni=k;

%  plot(xi,yi,'b.') %,xx(~in),yy(~in),'.b')
%  hold on;
%  plot(xb(1:nb/2),yb(1:nb/2),'r.')
%  hold on;
%  plot(xb(nb/2+1:nb),yb(nb/2+1:nb),'b.')
%  hold on;
%  plot(x(1:nbtest/2),y(1:nbtest/2),'k-')
%  hold on;
%  plot(x(nbtest/2:nbtest),y(nbtest/2:nbtest),'r-')

end