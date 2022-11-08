function [Fbp, Bbp, Ubp, Dbp, Lbp, Rbp, inp, N_N, bp, S] = geodata_rectangle(x, y, z, nx, ny, nz) 
% ========================  interior points ========
theta = linspace(0,2*pi,nt); 
phai = linspace(0,pi,np); 
rho=1+1/8*sin(10*theta).*sin(9*phai);
x=rho.*cos(theta).*sin(phai);
y=rho.*sin(theta).*sin(phai);
z=rho.*cos(phai);

[A, B, C] = meshgrid(x_Temp, y_Temp, z_Temp);

XX = A(:); 
YY = B(:);
ZZ = C(:);

inp = [XX YY ZZ];
%======================= boundary points ========
[AA, BB] = meshgrid(x_Temp, y_Temp);
BX = AA(:);
BY = BB(:);
BZ = z(1)*ones(length(BX),1);
Dbp = [BX BY BZ];    %   Down surface  下面

BZ = z(2)*ones(length(BX),1);
Ubp = [BX BY BZ];   %  Upper surface   上面
% ----------------------------  前后
[AA, BB] = meshgrid(y_Temp, z_Temp);
BY = AA(:);
BZ = BB(:);
BX = x(2)*ones(length(BY), 1);
Fbp = [BX BY BZ];    %   Front surface  前面

BX = x(1)*ones(length(BY), 1);
Bbp = [BX BY BZ];    %   Back surface  后面
% ----------------------------  左右
[AA, BB] = meshgrid(x_Temp, z_Temp);
BX = AA(:);
BZ = BB(:);
BY = y(1)*ones(length(BX), 1);
Lbp = [BX BY BZ];    %   Lift surface  左面

BY = y(2)*ones(length(BX), 1);
Rbp = [BX BY BZ];    %   Right surface  右面
% ===========================
Fn = zeros(length(Fbp),3);
Bn = zeros(length(Bbp),3);
Un = zeros(length(Ubp),3);
Dn = zeros(length(Dbp),3);
Ln = zeros(length(Lbp),3);
Rn = zeros(length(Rbp),3);

Fn(:,1) = 1;
Bn(:,1) = -1;

Un(:,3) = 1;
Dn(:,3) = -1;

Ln(:,2) = -1;
Rn(:,2) = 1;
N_N = [Fn; Bn; Un; Dn; Ln; Rn];

bp = [Fbp; Bbp; Ubp; Dbp; Lbp; Rbp];  %  前后上下左右
S = [inp; bp];
end





