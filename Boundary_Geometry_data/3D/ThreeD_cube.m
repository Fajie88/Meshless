%  Bottom up left right back front
function [bp, inp, S, N_flux] = ThreeD_cube(xscal,yscal,zscal,n1,n2,n3)
nxy=n1*n2; nxz=n1*n3; nyz=n2*n3;
dx=(xscal(2)-xscal(1))/(n1+1);
dy=(yscal(2)-yscal(1))/(n2+1);
dz=(zscal(2)-zscal(1))/(n3+1);
t1=linspace(xscal(1)+dx,xscal(2)-dx,n1);
t2=linspace(yscal(1)+dy,yscal(2)-dy,n2);
t3=linspace(zscal(1)+dz,zscal(2)-dz,n3);
[tx1,ty1] = meshgrid(t1,t2); 
[tx2,tz2] = meshgrid(t1,t3); 
[ty3,tz3] = meshgrid(t2,t3); 
bp1=[tx1(:), ty1(:), zscal(1)*ones(nxy,1)];
bp2=[tx1(:), ty1(:), zscal(2)*ones(nxy,1)];
bp3=[tx2(:), yscal(1)*ones(nxz,1),  tz2(:)];
bp4=[tx2(:), yscal(2)*ones(nxz,1),  tz2(:)];
bp5=[xscal(1)*ones(nyz,1), ty3(:), tz3(:)];
bp6=[xscal(2)*ones(nyz,1), ty3(:), tz3(:)];
bp=[bp1;bp2;bp3;bp4;bp5;bp6];

[xi,yi,zi] = meshgrid(t1,t2,t3);
inp=[xi(:),yi(:),zi(:)];

S = [inp; bp];
% plot3(bp(:,1),bp(:,2),bp(:,3),'.')
% hold on; plot3(xi(:),yi(:),zi(:),'r.')
% N_flux = [[zeros(nxy,1); zeros(nxy,1); zeros(nxz,1); zeros(nxz,1); -ones(nyz,1); ones(nyz,1)],...
%[zeros(nxy,1); zeros(nxy,1); -ones(nxz,1); ones(nxz,1); zeros(nyz,1); zeros(nyz,1)],...
%    [-ones(nxy,1); ones(nxy,1); zeros(nxz,1); zeros(nxz,1); zeros(nyz,1); zeros(nyz,1)]];

N_flux = zeros(2*(nxy+nxz+nyz),3);
N_flux(1:nxy,3)=-1; N_flux(nxy+1:2*nxy,3)=1;
N_flux(2*nxy+1:2*nxy+nxz,2)=-1; N_flux(2*nxy+nxz+1:2*(nxy+nxz),2)=1;
N_flux(2*(nxy+nxz)+1:2*(nxy+nxz)+nyz,1)=-1; N_flux(2*(nxy+nxz)+nyz+1:2*(nxy+nxz+nyz),1)=1;
end