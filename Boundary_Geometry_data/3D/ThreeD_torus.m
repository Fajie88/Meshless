
% torus
% n1: 大圆角度划分; n2: 小圆角度划分; n3: 小圆半径划分（内点）
function [bp, inp, S, N_flux] = ThreeD_torus(r1,r2,n1,n2,n3)
s1=linspace(0,2*pi,n1);
s2=linspace(0,2*pi,n2);
[u, v] = meshgrid(s1, s2);
x=(r1+r2*cos(v)).*cos(u);    
y=(r1+r2*cos(v)).*sin(u);  
z=r2*sin(v);
bp = [x(:), y(:), z(:)];
bp = unique(bp,'rows'); %怎样去掉其中的重复行呢？

surf(x, y, z);

xu=-(r1+r2*cos(v)).*sin(u); yu=(r1+r2*cos(v)).*cos(u); zu=0;
xv=-r2*sin(v).*cos(u); yv=-r2*sin(v).*sin(u); zv=r2*cos(v);
N_flux = [yu(:).*zv(:)-yv(:).*zu(:), zu(:).*xv(:)-zv(:).*xu(:), xu(:).*yv(:)-xv(:).*yu(:)];
N_flux = unique(N_flux,'rows'); %怎样去掉其中的重复行呢？

s1=linspace(0,2*pi,n1);
s2=linspace(0,2*pi,n2);
r22=linspace(0,r2-r2/(n3-1),n3);

[u, v, h] = meshgrid(s1, s2, r22);
x=(r1+h.*cos(v)).*cos(u);    
y=(r1+h.*cos(v)).*sin(u);  
z=h.*sin(v);
inp = [x(:), y(:), z(:)];
inp = unique(inp,'rows'); %怎样去掉其中的重复行呢？

S = [inp; bp];
% plot3(bp(:,1),bp(:,2),bp(:,3),'.')
% hold on; plot3(inp(:,1),inp(:,2),inp(:,3),'r.')
end