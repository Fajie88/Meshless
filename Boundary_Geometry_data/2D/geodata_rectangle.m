function [bp1,bp2,bp3,bp4,inp,N, bp, S] = geodata_rectangle(x, y, nx, ny)
bp1 = zeros(ny-2, 2);     % left boundary
bp2 = zeros(nx-2, 2);     % lower boundary
bp3 = zeros(ny-2, 2);     % right boundary 
bp4 = zeros(nx-2, 2);     % upper boundary

x_Temp = linspace(x(1),x(2),nx); 
y_Temp = linspace(y(1),y(2),ny);

bp1(:,1) = x(1);  
bp1(:,2) = y_Temp(2:ny-1);      

bp2(:,1) = x_Temp(2:nx-1);      
bp2(:,2) = y(1);

bp3(:,1) = x(2);                        
bp3(:,2) = y_Temp(2:ny-1);

bp4(:,1) = x_Temp(2:nx-1);      
bp4(:,2) = y(2);

n1 = zeros(ny-2,2);
n2 = zeros(nx-2,2);
n3 = zeros(ny-2,2);
n4 = zeros(nx-2,2);

n1(:,1) = -1;
n2(:,2) = -1;
n3(:,1) = 1;
n4(:,2) = 1;

N = [n1;n2;n3;n4];
% =============  interior points =============
[Allx,Ally] = meshgrid(x_Temp,y_Temp);   

XX = reshape(Allx(2:end-1,2:end-1),(nx-2)*(ny-2),1);
YY = reshape(Ally(2:end-1,2:end-1),(nx-2)*(ny-2),1);

inp = [XX YY];
% ====================== 内部随机布点=============
% dx = (x(2)-x(1))/(nx-1);
% dy = (y(2)-y(1))/(ny-1);
% rr = 2*rand(length(inp),1)-1;
% inp(:,1) = inp(:,1)+rr*dx/3;
% inp(:,2) = inp(:,2)+rr*dy/3;
%====================================
bp = [bp1;bp2;bp3;bp4];    %  whole boundary
S = [bp;inp];  %  the total points = boundary points + interior points
end