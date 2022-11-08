function [S, bp, N_N, b_xyz] = geometry_data_complex(point_all, boundary_index, cs)
%  ----------Boundary Geometry--------
%  cs = 1 (无人机  m = 10 w);  cs = 14 (4000);  
%  cs = 2 (教练机 m = 15 w); cs = 15 (2w); 
%  cs = 3 (核潜艇  12 w (close to singular)); cs = 18 (8 w); cs = 19 (4 w)
%  cs = 4 (客机  m = 13w);
%  cs = 8 (简单潜艇  m = 1w);  
%  cs = 99 (战斗机 m = 21w); cs = 999 (m = 3w); cs = 9999 (m = 1.5w)
%  cs = 51 (简易 船 m = 5w); cs = 551 (m = 7000)

%  cs = 5 (齿轮  m = 7w);  cs = 22 (m = 2w);  cs = 222 (m = 1w)
%  cs = 9 (汽车轮胎  m = 9w); cs = 13 (汽车轮胎  m = 3w); cs = 20 (m=2w)
%  cs = 10 (机械零件 [5通]  23w); cs = 11 (4w); cs = 12 (2w);  
%  cs = 17 (简易发动机  m = 6700)
%  cs = 33 (涡轮 4叶 m = 2w);  cs = 333 (6000);
%  cs = 123 (涡轮 5叶 m = 1w)
%  cs = 66 (涡轮 另一种 m = 1.2w)
%  cs = 77 (风车 m = 3.5w);
%  cs = 71 (蜂窝结构 m = 8000);

%  cs = 6 (海豚  m = 20w ); cs = 16 (6000)
%  cs = 7 (人头  m = 1w); 
%  cs = 55 (手  m = 5w); cs = 555 (m = 5000)
%  cs = 88 (整个人体 m = 1w)

%  cs = 10001 (CL201-1铰链  m = 3w)
%  cs = 10002 (愤怒的小鸟   m = 2w)
%  cs = 10003 (Bloc Motor de Aeromodelism 航空发动机   m = 5w)
%  cs = 10004 (Helical Gear 螺旋伞齿轮   m = 5w)
%  cs = 10005 (sprocket 链轮   m = 2.3w)
%  cs = 100061(房子 m = 2.56w); 10062 (m = 10w)
%  cs = 10007(车子 m = 7000); cs = 100071(车子 m = 1.2w); cs = 100072(m = 1.7w); 100073(m = 4w);cs = 100074(m = 8.3w);

point_all(:, 1) = [];
boundary_index(:, 1) = [];

if cs ==1 || cs == 14   %  无人机
    point_all = point_all/3;   %  图形缩小比例
elseif cs == 2 || cs == 15   %  教练机
    point_all = point_all/600;
elseif cs == 3 || cs == 18 || cs == 19 %  核潜艇
    point_all = point_all/800;
elseif cs == 4   %    客机
    point_all = point_all/160;
elseif cs == 5 || cs == 22 || cs == 222  %    
    point_all = point_all/160;
elseif cs == 6 || cs == 16   %    
    point_all = point_all/20;
elseif cs == 7   %    
    point_all = point_all/600;
elseif cs == 8   %    
    point_all = point_all/100;
elseif cs == 9 || cs == 13 || cs == 20   %    
    point_all = point_all/300;
elseif cs == 10 || cs == 11 || cs == 12  %    
    point_all = point_all/80;
elseif cs == 17   %    
    point_all = point_all/260;
elseif cs == 33 || cs == 333   %    
    point_all = point_all/3;
elseif cs == 123   %    
    point_all = point_all/30;
elseif cs == 55 || cs == 555   %    
    point_all = point_all/4;
elseif cs == 77 || cs == 777   %    
    point_all = point_all/800;
elseif cs == 88   %    
    point_all = point_all/30;
elseif cs == 99 || cs == 999  || cs == 9999  %    
    point_all = point_all/100;
elseif cs == 51 || cs == 551 %    
    point_all = point_all/30;
elseif cs == 10003   %    
    point_all = point_all/10;
elseif cs == 10007 || cs == 100071 || cs == 100072  || cs == 100073 || cs == 100074  %    
    point_all = point_all/100;
end

%----------  boundary nodes
A1 = [point_all(boundary_index(:, 1), 1), point_all(boundary_index(:, 1), 2), point_all(boundary_index(:, 1), 3)];
A2 = [point_all(boundary_index(:, 2), 1), point_all(boundary_index(:, 2), 2), point_all(boundary_index(:, 2), 3)];
A3 = [point_all(boundary_index(:, 3), 1), point_all(boundary_index(:, 3), 2), point_all(boundary_index(:, 3), 3)];

x1 = (A1(: ,1) + A2(: ,1) + A3(: ,1))/3;
x2 = (A1(: ,2) + A2(: ,2) + A3(: ,2))/3;
x3 = (A1(: ,3) + A2(: ,3) + A3(: ,3))/3;
bp = [x1, x2, x3];
%----------  N_N   Normal Vectors
vec_1 = A2 - A1;
vec_2 = A3 - A1;
N_N = cross(vec_1, vec_2);
N_N = N_N./sqrt(sum(N_N.^2, 2));
%----------  interior nodes
boundary_index = boundary_index(:);
boundary_index = unique(boundary_index);    %  去掉相同的值
%------- 边界所有三角形顶点的坐标
bounary_xyz = point_all(boundary_index, :);   %  tecplot 使用
%---------------------wfj
b_xyz=bounary_xyz(:);   
%-------------------
%-------
point_all_index = (1:length(point_all))';
inp_index = setdiff(point_all_index, boundary_index);   %  前一个向量中排除后一个向量的值
inp = point_all(inp_index, :);

bp=sortrows(bp,2);
S = [bp; inp];
end