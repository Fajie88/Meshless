function [u_exact, q_exact, f_source, ux_exact, uy_exact, uz_exact, eq_style] = example3D_Helmholtz(iexample)
%  -----------Exact solution----------
%         iexample=1: u_exact=sin(x).*cos(y).*cos(z);
%         iexample=2: u_exact=sin(x)+sin(y)+sin(z);
%         iexample=3:cos(x+2*y+2*z);
%         iexample=4:sin(x).*sin(y).*sin(z); 
%         iexample=5:exp(a1*x+a2*y+a3*z);
%         iexample=6: u_exact=cos(x+2*y+lamda*z); 
%         iexample=7: u_exact=exp(x+y+lamda*z)
global lamda
if iexample==1
   eq_style=1;
   lamda=3^0.5;  % wave number
   u_exact = @(x,y,z) sin(x).*cos(y).*cos(z);       %  u exact solutions
   q_exact = @(x,y,z,n_x,n_y,n_z) (cos(x).*cos(y).*cos(z)).*n_x+(-sin(x).*sin(y).*cos(z)).*n_y+(-sin(x).*cos(y).*sin(z)).*n_z;
   f_source = @(x,y,z) 0*x;    %  source term
   ux_exact = @(x,y,z) cos(x).*cos(y).*cos(z);        %  ux exact solutions
   uy_exact = @(x,y,z) -sin(x).*sin(y).*cos(z);       %  uy exact solutions
   uz_exact = @(x,y,z) -sin(x).*cos(y).*sin(z);       %  uz exact solutions
elseif iexample==2
   eq_style=1;
   lamda=1.0;  % wave number
   u_exact = @(x,y,z) sin(x)+sin(y)+sin(z);       %  exact solutions
   q_exact = @(x,y,z,n_x,n_y,n_z) cos(x).*n_x+cos(y).*n_y+cos(z).*n_z;
   f_source = @(x,y,z) 0*x;    %  source term
   ux_exact = @(x,y,z) cos(x);       %  ux exact solutions
   uy_exact = @(x,y,z) cos(y);       %  ux exact solutions
   uz_exact = @(x,y,z) cos(z);       %  ux exact solutions
elseif iexample==3    % Applied Mathematics and Computation, 2005, 165(2): 355-374.
   eq_style=1;       
   lamda=3.0;  % wave number
   u_exact = @(x,y,z) cos(x+2*y+2*z);       %  exact solutions
   q_exact = @(x,y,z,n_x,n_y,n_z) -sin(x+2*y+2*z).*n_x-2*sin(x+2*y+2*z).*n_y-2*sin(x+2*y+2*z).*n_z;
   f_source = @(x,y,z) 0*x;    %  source term
   ux_exact = @(x,y,z) -sin(x+2*y+2*z);       %  ux exact solutions
   uy_exact = @(x,y,z) -2*sin(x+2*y+2*z);       %  ux exact solutions
   uz_exact = @(x,y,z) -2*sin(x+2*y+2*z);       %  ux exact solutions
elseif iexample==4
   eq_style=1;
   lamda=3^0.5;  % wave number
   u_exact = @(x,y,z) sin(x).*sin(y).*sin(z);       %  u exact solutions
   q_exact = @(x,y,z,n_x,n_y,n_z) (cos(x).*sin(y).*sin(z)).*n_x+(sin(x).*cos(y).*sin(z)).*n_y+(sin(x).*sin(y).*cos(z)).*n_z;
   f_source = @(x,y,z) 0*x;    %  source term
   ux_exact = @(x,y,z) cos(x).*sin(y).*sin(z);        %  ux exact solutions
   uy_exact = @(x,y,z) sin(x).*cos(y).*sin(z);       %  uy exact solutions
   uz_exact = @(x,y,z) sin(x).*sin(y).*cos(z);       %  uz exact solutions
elseif iexample==5  % Applied Mathematics and Computation, 2005, 165(2): 355-374.
   eq_style=2;
   lamda=2.0;  % wave number
   a1=1.0; a2=0.5; a3=sqrt(lamda^2-a1^2-a2^2);
   u_exact = @(x,y,z) exp(a1*x+a2*y+a3*z);       %  u exact solutions
   q_exact = @(x,y,z,n_x,n_y,n_z) (a1*exp(a1*x+a2*y+a3*z)).*n_x+(a2*exp(a1*x+a2*y+a3*z)).*n_y+(a3*exp(a1*x+a2*y+a3*z)).*n_z;
   f_source = @(x,y,z) 0*x;    %  source term
   ux_exact = @(x,y,z) a1*exp(a1*x+a2*y+a3*z);        %  ux exact solutions
   uy_exact = @(x,y,z) a2*exp(a1*x+a2*y+a3*z);       %  uy exact solutions
   uz_exact = @(x,y,z) a3*exp(a1*x+a2*y+a3*z);       %  uz exact solutions
elseif iexample==6        %--- Note !!!!!!!  Nonhomogeneous,   MFS,BKM can not be directly used
   eq_style=2;
   lamda=2.0;  % wave number
   u_exact = @(x,y,z) cos(x+2*y+lamda*z);       %  exact solutions
   q_exact = @(x,y,z,n_x,n_y,n_z) -sin(x+2*y+lamda*z).*n_x-2*sin(x+2*y+lamda*z).*n_y-lamda*sin(x+2*y+lamda*z).*n_z;
   f_source = @(x,y,z) -5*cos(x+2*y+lamda*z);    %  source term
   ux_exact = @(x,y,z) -sin(x+2*y+lamda*z);       %  ux exact solutions
   uy_exact = @(x,y,z) -2*sin(x+2*y+lamda*z);       %  ux exact solutions
   uz_exact = @(x,y,z) -lamda*sin(x+2*y+lamda*z);       %  ux exact solutions
elseif iexample==7        %--- Note !!!!!!! Nonhomogeneous,   MFS,BKM can not be directly used
    eq_style=2;
   lamda=3^0.5;  % wave number
   u_exact = @(x,y,z) exp(x+y+lamda*z);       %  exact solutions
   q_exact = @(x,y,z,n_x,n_y,n_z) exp(x+y+lamda*z).*n_x+exp(x+y+lamda*z).*n_y+lamda*exp(x+y+lamda*z).*n_z;
   f_source = @(x,y,z) 2*exp(x+y+lamda*z);    %  source term
   ux_exact = @(x,y,z) exp(x+y+lamda*z);       %  ux exact solutions
   uy_exact = @(x,y,z) exp(x+y+lamda*z);       %  ux exact solutions
   uz_exact = @(x,y,z) lamda*exp(x+y+lamda*z);       %  ux exact solutions 
end
end