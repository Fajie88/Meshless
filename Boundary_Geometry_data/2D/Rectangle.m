function [xb,yb,xi,yi,nb,ni,pn1,pn2,lii] = Rectangle(x_scale,y_scale,irand,delta)
nix=100;
niy=nix;
nbx=100;
nby=nbx;
ni=nix*niy;
nb=2*(nbx+nby);
N=ni+nb;

x=linspace(x_scale(1),x_scale(2),nbx+1);
y=linspace(y_scale(1),y_scale(2),nby+1);
dx=(x_scale(2)-x_scale(1))/(nbx);
dy=(y_scale(2)-y_scale(1))/(nby);

lii=[dx.*ones(1,nbx) dy.*ones(1,nby) dx.*ones(1,nbx) dy.*ones(1,nby)];

ti=linspace(x_scale(1)+dx,x_scale(2)-dx,nbx);
tii=linspace(x_scale(2)-dx,x_scale(1)+dx,nbx);
tj=linspace(y_scale(1)+dy,y_scale(2)-dy,nby);
tjj=linspace(y_scale(2)-dy,y_scale(1)+dy,nby);

%±ß½çµã×ø±ê
A(:,1:nb)=[ti x_scale(2).*ones(1,nby) tii x_scale(1).*ones(1,nby);y_scale(1).*ones(1,nbx) tj y_scale(2).*ones(1,nbx) tjj];
pn1(1:nb)=[zeros(1,nbx) ones(1,nby) zeros(1,nbx) -ones(1,nby)];
pn2(1:nb)=[-ones(1,nbx) zeros(1,nby) ones(1,nbx) zeros(1,nby)];
pn1=pn1'; pn2=pn2';

NN=size(ti,2)*size(tj,2);
[xi,yi] = meshgrid(ti,tj);

B=[xi(:)';yi(:)'];

xp(1:ni)=B(1,1:ni);
yp(1:ni)=B(2,1:ni);

xp(ni+1:N)=A(1,1:nb);
yp(ni+1:N)=A(2,1:nb);

if irand==1

 xp = x_scale(1)+(x_scale(2)-x_scale(1))*rand(1,N);
 yp = y_scale(1)+(y_scale(2)-y_scale(1))*rand(1,N);
 
 yp(ni+1:ni+nb/4)=y_scale(1)*ones(1,nb/4);
 xp(ni+nb/4+1:ni+nb/2)=x_scale(2)*ones(1,nb/4);
 yp(ni+nb/2+1:ni+3*nb/4)=y_scale(2)*ones(1,nb/4);
 xp(ni+3*nb/4+1:N)=x_scale(1)*ones(1,nb/4);
%  for i=ni+1:N
%      if i<=ni+nb/4
%        yp(i)=y_scale(1);
%      elseif i<=ni+nb/2
%        xp(i)=x_scale(2);
%      elseif i<=ni+3*nb/4
%        yp(i)=y_scale(2);
%      else
%        xp(i)=x_scale(1);
%      end
%  end
 
elseif irand==2

 %nodePosition1(1,:)=[1,1];
 %nodePosition1(Nx,:)=[3,1];
 %nodePosition1((Ny-1)*Nx+1,:)=[1,3];
 %nodePosition1(Nx*Ny,:)=[3,3];
 
xp(1:ni)=xp(1:ni)+dx*delta*(2*rand(1,ni)-1);
yp(1:ni)=yp(1:ni)+dy*delta*(2*rand(1,ni)-1);

xp(ni+1:ni+nb/4)=xp((ni+1:ni+nb/4))+dx*delta*rand(1,nb/4);
yp(ni+nb/4+1:ni+nb/2)=yp(ni+nb/4+1:ni+nb/2)+dy*delta*rand(1,nb/4);
xp(ni+nb/2+1:ni+3*nb/4)=xp(ni+nb/2+1:ni+3*nb/4)+dx*delta*rand(1,nb/4);
yp(ni+3*nb/4+1:N)=yp(ni+3*nb/4+1:N)+dy*delta*rand(1,nb/4);

%  for j=1:N
%      if j<=ni
%        xp(j)=xp(j)+dx*delta*rand;
%        yp(j)=yp(j)+dy*delta*rand;
%      elseif j<=ni+nb/4
%        xp(j)=xp(j)+dx*delta*rand;
%      elseif j<=ni+nb/2
%        yp(j)=yp(j)+dy*delta*rand;
%      elseif j<=ni+3*nb/4
%        xp(j)=xp(j)+dx*delta*rand;
%      else
%        yp(j)=yp(j)+dy*delta*rand;
%      end           
%  end
 
end
xi=xp(1:ni);yi=yp(1:ni);
xb=xp(ni+1:N);yb=yp(ni+1:N);

xxbb=[x_scale(1) x_scale(2) x_scale(2) x_scale(1) x_scale(1)];
yybb=[y_scale(1) y_scale(1) y_scale(2) y_scale(2) y_scale(1)];
plot(xxbb,yybb);

end