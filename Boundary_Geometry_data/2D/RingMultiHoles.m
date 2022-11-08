function [xb,yb,xi,yi,nb,ni,pn1,pn2]=RingMultiHoles(nc,nb1,nb2,nb3,irand,delta,x0,x1,y0,y1,nx,ny)
   r1=1;
   r2=0.4;
   rd=r1-r2;
   r3=0.1;%1/5*rd;
   nbtest=200;
   nb=nb1+nb2+nc*nb3;
%Boundary points
    theta1=linspace(0,2*pi-2*pi/(nb1-1),nb1);
    theta2=linspace(0,2*pi-2*pi/(nb2-1),nb2);
    theta3=linspace(0,2*pi-2*pi/(nb3-1),nb3);

    xb(1:nb1)=r1.*cos(theta1(1:nb1));
    yb(1:nb1)=r1.*sin(theta1(1:nb1));
    pn1(1:nb1)=cos(theta1(1:nb1));
    pn2(1:nb1)=sin(theta1(1:nb1));
    
    xb(nb1+1:nb1+nb2)=r2.*cos(theta2(1:nb2));
    yb(nb1+1:nb1+nb2)=r2.*sin(theta2(1:nb2));
    pn1(nb1+1:nb1+nb2)=-cos(theta2(1:nb2));
    pn2(nb1+1:nb1+nb2)=-sin(theta2(1:nb2));   
    
    k=nb1+nb2;
    for i=1:nc
        thb(i)=2*pi*(i-1)/nc;
        xo(i)=(r2+rd/2)*cos(thb(i));
        yo(i)=(r2+rd/2)*sin(thb(i));
        for j=1:nb3
            k=k+1;
            thc(j)=2*pi*(j-1)/nb3;
            xb(k)=xo(i)+r3*cos(thc(j));
            yb(k)=yo(i)+r3*sin(thc(j));
            pn1(k)=-cos(thc(j));
            pn2(k)=-sin(thc(j));
        end
    end
pn1=pn1'; pn2=pn2';
 plot(xb,yb,'*')
    
% Inner points
    theta1=linspace(0,2*pi,nbtest);
    theta2=linspace(0,2*pi,nbtest);
    theta3=linspace(0,2*pi,nbtest);
    
     
 thb=2*pi.*(0:nc-1)/nc;
 xo=(r2+rd/2).*cos(thb);
 yo=(r2+rd/2).*sin(thb);
 xx=zeros(nc,nbtest);
 yy=zeros(nc,nbtest); 
 xx(1,:)=r1*cos(theta1);
 yy(1,:)=r1*sin(theta1);
 xx(2,:)=r2*cos(theta2);
 yy(2,:)=r2*sin(theta2); 
 for kkk=1:nc
     xx(kkk+2,:)=xo(kkk)+r3*cos(theta3);
     yy(kkk+2,:)=yo(kkk)+r3*sin(theta3);
 end
    
%     kk=2*nbtest;
%     for i=1:nc
%         thb(i)=2*pi*(i-1)/nc;
%         xo(i)=(r2+rd/2)*cos(thb(i));
%         yo(i)=(r2+rd/2)*sin(thb(i));
%         xx(2+i)=xo(i)+r3*cos(theta3);
%         yy(2+i)=yo(i)+r3*sin(theta3);
%  
%     end
    
    
    x(1:nbtest)=r1.*cos(theta1(1:nbtest));
    y(1:nbtest)=r1.*sin(theta1(1:nbtest));
    x(nbtest+1:2*nbtest)=r2.*cos(theta2(1:nbtest));
    y(nbtest+1:2*nbtest)=r2.*sin(theta2(1:nbtest));
    kk=2*nbtest;
    for i=1:nc
        thb(i)=2*pi*(i-1)/nc;
        xo(i)=(r2+rd/2)*cos(thb(i));
        yo(i)=(r2+rd/2)*sin(thb(i));
        for j=1:nbtest
            kk=kk+1;
            x(kk)=xo(i)+r3*cos(theta3(j));
            y(kk)=yo(i)+r3*sin(theta3(j));
        end
    end


ti=linspace(x0,x1,nx);
tj=linspace(y0,y1,ny);
NN=size(ti,2)*size(tj,2);
dx=(x1-x0)/nx;
dy=(y1-y0)/ny;
[xt,yt] = meshgrid(ti,tj);
B=[xt(:)';yt(:)'];
xxx=B(1,:);
yyy=B(2,:);
if irand==1
xxx=xxx;
yyy=yyy;
elseif irand==2
    for i=1:NN
xxx(i)=xxx(i)+dx*delta*(2*rand-1);
yyy(i)=yyy(i)+dy*delta*(2*rand-1);
    end
elseif irand==3
    for i=1:NN
xxx(i)=x0+(x1-x0)*rand;
yyy(i)=y0+(y1-y0)*rand;
    end
end

for k1=1:nc+2
in(k1,:)=inpolygon(xxx,yyy,xx(k1,:),yy(k1,:));
end


k=0;
%for i=1:nc
for j=1:NN
    if in(1,j)==1 & in(2:nc+2,j)==0
        k=k+1;
        xi(k)=xxx(j);
        yi(k)=yyy(j);
    else
    end

end
%end
ni=k;

plot(xi,yi,'b.')%,xx(~in),yy(~in),'.b')
hold on;
plot(xb,yb,'r.')


for i=1:2+nc
k=nbtest*(i-1);
hold on;
plot(x(k+1:k+nbtest),y(k+1:k+nbtest),'k-')
end

end