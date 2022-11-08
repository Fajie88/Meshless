%  
%
clear all;
clc;
addpath('E:\�ҵĳ���\BoundaryGeometry\2D');
addpath('E:\�ҵĳ���\BoundaryGeometry\3D');
%
% % --------  ��ά�߽���״���߽�ڵ����ꡢ������ --------------
% iboundary=1; %--Types of boundary geometry: retangular,irregular,simply connected multiple connected domain
% idomain=1;   %--Computational domian
% %----------  HOW TO USE IT ?  -----------
% [inp, bp, S, normal_vector, N] = BoundaryGeometry(iboundary, idomain);


% --------  ��ά�߽���״���߽�ڵ����ꡢ������ --------------
xscal = [0 1];  yscal = [0 1];  zscal = [0 1]; nx = 100; ny = nx; nz = nx; 
tic;
% [bp1, bp2, bp3, bp4, bp5, bp6, inp, N_flux, bp, S] = geodata_cube(xscal, yscal, zscal, nx, ny, nz); 
%[bp, inp, S, N_flux] = ThreeD_cube(xscal, yscal, zscal, nx, ny, nz); 
%r=2; n=10; [bp, inp, S, N_flux] = ThreeD_sphere(r,n);
%r1=4; r2=2; n1=50; n2=50; n3=20; [bp, inp, S, N_flux] = ThreeD_torus(r1,r2,n1,n2,n3);
%za=0; zb=5; r=2; n1=20; n2=10; n3=20; [bp, inp, S, N_flux] = ThreeD_cylinder(za,zb,r,n1,n2,n3);
r=1; n1=50; n2=50; n3=10; [bp, inp, S, N_flux] = ThreeD_IrregularBall(r, n1, n2, n3);
t=toc
plot3(bp(:,1),bp(:,2),bp(:,3),'r.');
hold on; plot3(inp(:,1),inp(:,2),inp(:,3),'k.');
