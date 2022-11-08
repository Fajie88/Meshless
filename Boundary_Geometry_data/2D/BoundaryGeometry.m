function [inp, bp, S, normal_vector, N_total,lii]=BoundaryGeometry(iboundary, idomain)
%----------Boundary Geometry------------------
%       --iboundary=1: Rectangle
%       --iboundary=2: SimplyConnDomain
%             idomain=1: Circle
%             idomain=2: Gear wheel shape
%             idomain=3: Amoeba-like domain
%             idomain=4: Peanut1
%             idomain=5: Peanut2
%             idomain=6: Triangle 
%             idomain=7: Heart
%             idomain=8: Gear-shaped domain
%             idomain=9: Epitrochoid boundary shape
%             idomain=10: Petal
%             idomain=11: Irregular
%       --iboundary=3: MultiConnDomain
%             idomain=1: Peanut with two petals 
%       --iboundary=4: MultiConnDomain
%             idomain=1: Ring with nc holes (nc=10\15\20)
%
% addpath('F:\Matlab_Subroutine\Boundary_Geometry_data');
%===================== geometry=====
irand=2;       % Style of collocation nodes
delta=0.0;     % Perturbation for Irregular I
x0=-6 ;  x1=-x0;  y0=x0;  y1=-x0;
nx=150;  ny=nx;
if iboundary==1
    x_scale=[0 1]; y_scale=[0 1];
   [xb,yb,xi,yi,nb,ni,pn1,pn2,lii] = Rectangle(x_scale,y_scale,irand,delta);   
%    xscal = [0 1]; yscal = [0 1]; nx = 52;    ny = nx;
%    [bp1,bp2,bp3,bp4,inp,N,bp,S] = geodata_rectangle(xscal, yscal, nx, ny); 
elseif iboundary==2
   nb=200;
   [xb,yb,xi,yi,ni,pn1,pn2,lii] = IrregularBoundary(nb,idomain,irand,delta,x0,x1,y0,y1,nx,ny);
elseif iboundary==3
   nb1=600; nb2=50; nb3=200; nb4=200;
   [xb,yb,xi,yi,nb,ni,pn1,pn2] = MultiConnBoundary(nb1,nb2,nb3,nb4,irand,delta,x0,x1,y0,y1,nx,ny);
   lii=[];
elseif iboundary==4
   nc=15;
   nb1=100; nb2=50; nb3=10;
   [xb,yb,xi,yi,nb,ni,pn1,pn2] = RingMultiHoles(nc,nb1,nb2,nb3,irand,delta,x0,x1,y0,y1,nx,ny);
   lii=[];
elseif iboundary==5
   nb1=400; nb2=400;  
   [xb,yb,xi,yi,nb,ni,pn1,pn2,lii]=MultiConnBoundary_1(nb1,nb2,irand,delta,x0,x1,y0,y1,nx,ny);
end
bp=[xb' yb']; inp=[xi' yi'];
normal_vector=[pn1 pn2];
S = [bp;inp];
N_total = length(S);

% figure(1)
% plot(xi,yi,'b.')
% hold on;
% plot(xb,yb,'r.')
end