% cylinder
function [bp, inp, S, N_flux] = ThreeD_cylinder(za,zb,r,n1,n2,n3)

stt=linspace(0,2*pi,n1);
roo=linspace(0,r,n2);   %------------为了棱角上的法向量
zoo=linspace(za,zb,n3);
[theta, zzz] = meshgrid(stt, zoo);
x1=r.*cos(theta);    
y1=r.*sin(theta);   
z1=zzz;

[rho, theta] = meshgrid(roo, stt);
x2=rho.*cos(theta);    
y2=rho.*sin(theta);   
z2=za*ones(size(theta));

x3=rho.*cos(theta);    
y3=rho.*sin(theta);   
z3=zb*ones(size(theta));

Cbp = [x1(:), y1(:), z1(:)]; 
Dbp = [x2(:), y2(:), z2(:)]; 
Ubp = [x3(:), y3(:), z3(:)]; 
bpx = [Cbp; Dbp; Ubp];
[bp,cc,~] = unique(bpx,'rows'); %怎样去掉其中的重复行呢？

Cn = Cbp/r;
Dn = zeros(length(Dbp),3);
Un = zeros(length(Ubp),3);
Dn(:,3) = -1;
Un(:,3) = 1;
N_N = [Cn; Dn; Un];
N_flux = N_N(cc,:);
% N_flux = unique(N_flux,'rows'); %怎样去掉其中的重复行呢？

% plot3(bp(:,1),bp(:,2),bp(:,3),'.');
% figure; surf(x1, y1, z1); hold on; surf(x2, y2, z2); hold on; surf(x3, y3, z3);

stt=linspace(0,2*pi,n1);
roo=linspace(0,r-r/(n2-1),n2);  
zoo=linspace(za+(zb-za)/(n3-1),zb-(zb-za)/(n3-1),n3);
[u, v, h] = meshgrid(stt, zoo, roo);
x=h.*cos(u);    
y=h.*sin(u);  
z=v;
inp = [x(:), y(:), z(:)];
inp = unique(inp,'rows'); %怎样去掉其中的重复行呢？

S = [inp; bp];

% plot3(bp(:,1),bp(:,2),bp(:,3),'.')
% hold on; plot3(inp(:,1),inp(:,2),inp(:,3),'r.')
end

