% IrregularBall
function [bp, inp, S, N_flux] = ThreeD_IrregularBall(r, n1, n2, ni)
% n1=50; n2=50; ni=20;
a=1/10; b=6; c=5;
st=linspace(0,2*pi,n1);
fa=linspace(0,pi,n2);
[s, f] = meshgrid(st, fa);
rho=r+a*sin(b*s).*sin(c*f);
x=rho.*cos(s).*sin(f);    
y=rho.*sin(s).*sin(f); 
z=rho.*cos(f);
bp = [x(:), y(:), z(:)];
bp = unique(bp,'rows'); %怎样去掉其中的重复行呢？

% surf(x, y, z);

xs=a*b*cos(b*s).*sin(c*f).*cos(s).*sin(f) - sin(f).*sin(s).*(r + a*sin(c*f).*sin(b*s));
ys=cos(s).*sin(f).*(r + a*sin(c*f).*sin(b*s)) + a*b*cos(b*s).*sin(c*f).*sin(f).*sin(s); 
zs=a*b*cos(b*s).*sin(c*f).*cos(f);
xf=cos(f).*cos(s).*(r + a*sin(c*f).*sin(b*s)) + a*c*cos(c*f).*sin(b*s).*cos(s).*sin(f); 
yf=cos(f).*sin(s).*(r + a*sin(c*f).*sin(b*s)) + a*c*cos(c*f).*sin(b*s).*sin(f).*sin(s); 
zf=a*c*cos(c*f).*sin(b*s).*cos(f) - sin(f).*(r + a*sin(c*f).*sin(b*s));
N_flux = [ys(:).*zf(:)-yf(:).*zs(:), zs(:).*xf(:)-zf(:).*xs(:), xs(:).*yf(:)-xf(:).*ys(:)];
N_flux = unique(N_flux,'rows'); %怎样去掉其中的重复行呢？

t=linspace(-r-a,r+a,ni);
[xx, yy, zz] = meshgrid(t,t,t);
X = reshape(xx,ni^3,1);Y = reshape(yy,ni^3,1);Z = reshape(zz,ni^3,1);
distP2PMat = sqrt(X(:).^2+Y(:).^2+Z(:).^2);
ptsInSphere = find(distP2PMat < r-2*a);
xi = X(ptsInSphere,:);yi = Y(ptsInSphere,:);zi = Z(ptsInSphere,:);
inp=[xi(:),yi(:),zi(:)];
inp = unique(inp,'rows'); %怎样去掉其中的重复行呢？

S = [inp; bp];
% plot3(bp(:,1),bp(:,2),bp(:,3),'.')
% hold on; plot3(inp(:,1),inp(:,2),inp(:,3),'r.')
end