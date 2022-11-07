function [G0,E0] = GE0(r,x1,x2,s1,s2,n_s1,n_s2)
%GE0 此处显示有关此函数的摘要
%   此处显示详细说明
G0 =-1/(2*pi)*log(r);
G0_ds1 = 1/(2*pi)*(x1-s1)./r.^2;
G0_ds2 = 1/(2*pi)*(x2-s2)./r.^2;
E0 = G0_ds1.*n_s1+G0_ds2.*n_s2;  %G0_ns
end

