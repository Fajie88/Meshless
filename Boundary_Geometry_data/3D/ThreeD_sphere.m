% sphere
function [bp, inp, S, N_flux] = ThreeD_sphere(r,n)
ni=10;
[X,Y,Z] = sphere(n);
bpp = [X(:),Y(:),Z(:)]; 
bp = unique(r*bpp,'rows'); %怎样去掉其中的重复行呢？

% length(bp)
figure;
surf(X, Y, Z);

t=linspace(-r,r,ni);
[xx, yy, zz] = meshgrid(t,t,t);
X = reshape(xx,ni^3,1);Y = reshape(yy,ni^3,1);Z = reshape(zz,ni^3,1);
distP2PMat = sqrt(X(:).^2+Y(:).^2+Z(:).^2);
ptsInSphere = find(distP2PMat <= r);
xi = X(ptsInSphere,:);yi = Y(ptsInSphere,:);zi = Z(ptsInSphere,:);
inp=[xi(:),yi(:),zi(:)];

S = [inp; bp];
% plot3(bp(:,1),bp(:,2),bp(:,3),'.')
% hold on; plot3(xi(:),yi(:),zi(:),'r.')
N_flux =bp./r;
end