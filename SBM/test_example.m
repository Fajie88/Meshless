function [u_exact,ux_exact,uy_exact,f_source,lamda] = test_example(type_eq, iexample) 

% type_eq: 'Laplace', 'Poisson', 'Helmholtz', 'General'

%--------------Exact solution----------Laplace Equations------
%         iexample=1: u_exact=exp(x).*sin(y);
%         iexample=2: u_exact=exp(x).*sin(y)+exp(y).*sin(x);
%         iexample=3: u_exact=x.^2-y.^2;
%         iexample=4: u_exact=x.^2-y.^2+(x-0.5).*(y-0.5);
%         iexample=5: u_exact=cos(x+y)+x.^3+y.^3+1;
%         iexample=6: u_exact=cos(x).*sinh(y)+x.^2-y.^2+x.*y+1;
%         iexample=7: u_exact=cos(x).*cosh(y)+sin(x).*sinh(y);
%         iexample=8: u_exact=exp(x).*cos(y); 
%--------------Exact solution---------Helmholtz(modified) Equations------
%         iexample=1: u_exact=sin(x).*sin(y);
%         iexample=2: u_exact=x.*sin(sqrt(2)*y);
%         iexample=3: u_exact=sin(10*x)+sin(10*y);
%         iexample=4: u_exact=sin(sqrt(2)*x).*sinh(y)+cos(y); 
%         iexample=1: u_exact=sin(x).*cosh(sqrt(3)*y)+cos(x).*sinh(sqrt(3)*y); %modified Helmholtz equation
%         iexample=2: u_exact=exp(2*x-y); %modified Helmholtz equation
%--------------Exact solution---------Poisson Equations------
%         iexample=1: u_exact=x.^2+y.^2+exp(x+y);  
%         iexample=2: u_exact=x.^3.*y-2*x.*y.^3+x+10;   Eng. Anal. Bound. Elem. 29 (2005) 756¨C760
%         iexample=3: u_exact=exp(2*x+2*y);             Eng. Anal. Bound. Elem. 28 (2004) 1417¨C1425
%         iexample=4: u_exact=sin(pi*x/6).*sin(7*pi*x/4).*sin(3*pi*y/4).*sin(5*pi*y/4)   Eng. Andy. Boundary Elements, 303-311, 1994.
%--------------Exact solution---------General Equations------
%         iexample=1: u_exact=x.^2+y.^2-5*x.*y;     Eng. Anal. Bound. Elem. 24 (2000) 549-557.
%         iexample=2: u_exact=exp(x).*sin(y);

if strcmp(type_eq,'Laplace')
    f_source = @(x,y) 0.0;
    if iexample==1
        u_exact = @(x,y) exp(x).*sin(y);
        ux_exact = @(x,y) exp(x).*sin(y);
        uy_exact = @(x,y) exp(x).*cos(y);
    elseif iexample==2
        u_exact = @(x,y) exp(x).*sin(y)+exp(y).*sin(x);
        ux_exact = @(x,y) exp(x).*sin(y)+exp(y).*cos(x);
        uy_exact = @(x,y) exp(x).*cos(y)+exp(y).*sin(x);
    elseif iexample==3
        u_exact = @(x,y) x.^2-y.^2;
        ux_exact = @(x,y) 2*x;
        uy_exact = @(x,y) -2*y;
    elseif iexample==4
        u_exact = @(x,y) x.^2-y.^2+(x-0.5).*(y-0.5);
        ux_exact = @(x,y) 2*x+(y-0.5);
        uy_exact = @(x,y) -2*y+(x-0.5);
    elseif iexample==5  % non-harmonic boundary condition
        u_exact = @(x,y) cos(x+y)+x.^3+y.^3+1;
        ux_exact = @(x,y) -sin(x+y)+3*x.^2;
        uy_exact = @(x,y) -sin(x+y)+3*y.^2;
    elseif iexample==6
        u_exact = @(x,y) cos(x).*sinh(y)+x.^2-y.^2+x.*y+1;
        ux_exact = @(x,y) -sin(x).*sinh(y)+2*x+y;
        uy_exact = @(x,y) cos(x).*cosh(y)-2*y+x;
    elseif iexample==7
        u_exact = @(x,y) cos(x).*cosh(y)+sin(x).*sinh(y);
        ux_exact = @(x,y) -sin(x).*cosh(y)+cos(x).*sinh(y);
        uy_exact = @(x,y) cos(x).*sinh(y)+sin(x).*cosh(y);
    elseif iexample==8
        u_exact = @(x,y) exp(x).*cos(y);
        ux_exact = @(x,y) exp(x).*cos(y);
        uy_exact = @(x,y) -exp(x).*sin(y);
    end
    
elseif strcmp(type_eq,'Helmholtz')
    f_source = @(x,y) 0.0;
    if iexample==1
        u_exact = @(x,y) sin(x).*sin(y);
        ux_exact = @(x,y) cos(x).*sin(y);
        uy_exact = @(x,y) sin(x).*cos(y); 
    elseif iexample==2
        u_exact = @(x,y) x.*sin(sqrt(2)*y);
        ux_exact = @(x,y) sin(sqrt(2)*y);
        uy_exact = @(x,y) sqrt(2)*x.*cos(sqrt(2)*y);
    elseif iexample==3
        u_exact = @(x,y) sin(10*x)+sin(10*y);
        ux_exact = @(x,y) 10*cos(10*x);
        uy_exact = @(x,y) 10*cos(10*y);
    elseif iexample==4 
        u_exact = @(x,y) sin(sqrt(2)*x).*sinh(y)+cos(y);
        ux_exact = @(x,y) sqrt(2)*cos(sqrt(2)*x).*sinh(y);
        uy_exact = @(x,y) sin(sqrt(2)*x).*cosh(y)-sin(y); 
    end
    
elseif strcmp(type_eq,'ModifiedHelmholtz')
    f_source = @(x,y) 0.0;
    if iexample==1   %modified Helmholtz equation
        u_exact = @(x,y) sin(x).*cosh(sqrt(3)*y)+cos(x).*sinh(sqrt(3)*y);
        ux_exact = @(x,y) cos(x).*cosh(sqrt(3)*y)-sin(x).*sinh(sqrt(3)*y);
        uy_exact = @(x,y) sqrt(3)*sin(x).*sinh(sqrt(3)*y)+sqrt(3)*cos(x).*cosh(sqrt(3)*y);
    elseif iexample==2  %modified Helmholtz equation
        u_exact = @(x,y) exp(2*x-y);
        ux_exact = @(x,y) 2*exp(2*x-y);
        uy_exact = @(x,y) -exp(2*x-y);
    end
    
elseif strcmp(type_eq,'Poisson')   
    if iexample==1
        u_exact = @(x,y) x.^2+y.^2+exp(x+y);
        ux_exact = @(x,y) 2*x+exp(x+y);
        uy_exact = @(x,y) 2*y+exp(x+y);
        f_source = @(x,y) 4+2*exp(x+y);
    elseif iexample==2
        u_exact = @(x,y) x.^3.*y-2*x.*y.^3+x+10;
        ux_exact = @(x,y) 3*x.^2.*y-2*y.^3+1;
        uy_exact = @(x,y) x.^3-6*x.*y.^2;
        f_source = @(x,y) -6*x.*y;
    elseif iexample==3
        u_exact = @(x,y) exp(2*x+2*y);
        ux_exact = @(x,y) 2*exp(2*x+2*y);
        uy_exact = @(x,y) 2*exp(2*x+2*y);
        f_source = @(x,y) 8*exp(2*x+2*y);
    elseif iexample==4
        u_exact = @(x,y)   sin(pi*x/6).*sin(7*pi*x/4).*sin(3*pi*y/4).*sin(5*pi*y/4); 
        ux_exact = @(x,y) (pi*cos((pi*x)/6).*sin((7*pi*x)/4).*sin((3*pi*y)/4).*sin((5*pi*y)/4))/6 +...
                (7*pi*cos((7*pi*x)/4).*sin((pi*x)/6).*sin((3*pi*y)/4).*sin((5*pi*y)/4))/4;
        uy_exact = @(x,y) (3*pi*cos((3*pi*y)/4).*sin((pi*x)/6).*sin((7*pi*x)/4).*sin((5*pi*y)/4))/4 +...
                 (5*pi*cos((5*pi*y)/4).*sin((pi*x)/6).*sin((7*pi*x)/4).*sin((3*pi*y)/4))/4;
        f_source = @(x,y) (-751*pi^2/144).*sin(pi*x/6).*sin(7*pi*x/4).*sin(3*pi*y/4).*sin(5*pi*y/4) +...
                 (7*pi^2/12).*cos(pi*x/6).*cos(7*pi*x/4).*sin(3*pi*y/4).*sin(5*pi*y/4) +...
                 (15*pi^2/8).*sin(pi*x/6).*sin(7*pi*x/4).*cos(3*pi*y/4).*cos(5*pi*y/4);
    end
    
elseif strcmp(type_eq,'General')   
    if iexample==1
        u_exact = @(x,y) x.^2+y.^2-5*x.*y;
        ux_exact = @(x,y) 2*x-5*y;
        uy_exact = @(x,y) 2*y-5*x;
        f_source = @(x,y) -18*(x.^2+y.^2)+16;
    elseif iexample==2
        u_exact = @(x,y) exp(x).*sin(y);
        ux_exact = @(x,y) exp(x).*sin(y);
        uy_exact = @(x,y) exp(x).*cos(y);
        f_source = @(x,y) exp(x).*sinh(x).*cosh(y).*cos(y);
    end
end

%====wave number for Helmholtz or ModifiedHelmholtz
if strcmp(type_eq, 'Helmholtz')==1  && iexample==1
    lamda=sqrt(2);
elseif strcmp(type_eq, 'Helmholtz')==1 && iexample==2
    lamda=sqrt(2);
elseif strcmp(type_eq, 'Helmholtz')==1 && iexample==3
    lamda=10;   
elseif strcmp(type_eq, 'Helmholtz')==1 && iexample==4
    lamda=1;
elseif strcmp(type_eq, 'ModifiedHelmholtz')==1 && iexample==1
    lamda=sqrt(2);
elseif strcmp(type_eq, 'ModifiedHelmholtz')==1 && iexample==2
    lamda=sqrt(5);
else
    lamda=[];  
end

end