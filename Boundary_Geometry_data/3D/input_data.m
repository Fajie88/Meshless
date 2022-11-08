function [point_all, boundary_index] = input_data(cs)
%  ----------Boundary Geometry--------
%  cs = 1 (���˻�  m = 10 w);  cs = 14 (4000);  
%  cs = 2 (������ m = 15 w); cs = 15 (2w); 
%  cs = 3 (��Ǳͧ  12 w (close to singular)); cs = 18 (8 w); cs = 19 (4 w)
%  cs = 4 (�ͻ�  m = 13w);
%  cs = 8 (��Ǳͧ  m = 1w);  
%  cs = 99 (ս���� m = 21w); cs = 999 (m = 3w); cs = 9999 (m = 1.5w)
%  cs = 51 (���� �� m = 5w); cs = 551 (m = 7000)

%  cs = 5 (����  m = 7w);  cs = 22 (m = 2w);  cs = 222 (m = 1w)
%  cs = 9 (������̥  m = 9w); cs = 13 (������̥  m = 3w); cs = 20 (m=2w)
%  cs = 10 (��е��� [5ͨ]  23w); cs = 11 (4w); cs = 12 (2w);  
%  cs = 17 (���׷�����  m = 6700)
%  cs = 33 (���� 4Ҷ m = 2w);  cs = 333 (6000);
%  cs = 123 (���� 5Ҷ m = 1w)
%  cs = 66 (���� ��һ�� m = 1.2w)
%  cs = 77 (�糵 m = 3.5w);
%  cs = 71 (���ѽṹ m = 8000);

%  cs = 6 (����  m = 20w ); cs = 16 (6000)
%  cs = 7 (��ͷ  m = 1w); 
%  cs = 55 (��  m = 5w); cs = 555 (m = 5000)
%  cs = 88 (�������� m = 1w)

%  cs = 10001 (CL201-1����  m = 3w)
%  cs = 10002 (��ŭ��С��   m = 2w)
%  cs = 10003 (Bloc Motor de Aeromodelism ���շ�����   m = 5w)
%  cs = 10004 (Helical Gear ����ɡ����   m = 5w)
%  cs = 10005 (sprocket ����   m = 2.3w)
%  cs = 100061(���� m = 2.56w); 10062 (m = 10w)
%  cs = 10007(���� m = 7000); cs = 100071(���� m = 1.2w); cs = 100072(m = 1.7w); 100073(m = 4w);cs = 100074(m = 8.3w);
addpath('F:\Matlab_Subroutine\Model_Yan_Wang');
if cs == 1   %  ���˻�  m = 10 ��
    point_all = load('wurenji.txt');  %  ���е������(�߽�+�ڵ�)
    boundary_index = load('wurenji_boundary_index.txt'); %   �߽�С��Ԫ����������ţ����ظ���ţ�
elseif cs == 2  %   ������ m = 15 ��
    point_all = load('jiaolianji.txt');  
    boundary_index = load('jiaolianji_boundary_index.txt'); 
elseif cs == 3    %   ��Ǳͧ  m = 12 �� (close to singular)
    point_all = load('heqianting.txt');  
    boundary_index = load('heqianting_boundary_index.txt'); 
elseif cs == 4     %    �ͻ�  m = 13��
    point_all = load('air_plane.txt');  
    boundary_index = load('air_plane_boundary_index.txt'); 
elseif cs == 5     %    ����  m = 7��
    point_all = load('chilun.txt');  
    boundary_index = load('chilun_boundary_index.txt'); 
elseif cs == 6     %    ����  m = 20�� �����²��㣬���ٵ�10�����ڣ�
    point_all = load('dolphin.txt');  
    boundary_index = load('dolphin_boundary.txt'); 
elseif cs == 7     %    ��ͷ  m = 1��
    point_all = load('humanhead_7349.txt');  
    boundary_index = load('humanhead_7349_boundary.txt'); 
elseif cs == 8     %    ��Ǳͧ  m = 1��
    point_all = load('submarine.txt');  
    boundary_index = load('submarine_boundary_index.txt'); 
elseif cs == 9     %    ������̥  m = 9��
    point_all = load('rim18.txt');  
    boundary_index = load('rim18_bounary_index.txt'); 
elseif cs == 10     %    ��е��� 01  m = 23��
    point_all = load('tool.txt');  
    boundary_index = load('tool_boundary_index.txt'); 
elseif cs == 11     %    ��е��� 01  m = 4��
    point_all = load('tool_1_30000.txt');  
    boundary_index = load('tool_1_30000_bounary_index.txt'); 
elseif cs == 12     %    ��е��� 01  m = 2��
    point_all = load('tool_1_6000.txt');  
    boundary_index = load('tool_1_6000_bounary_index.txt'); 
elseif cs == 13     %    ������̥  m = 3��
    point_all = load('tool_2.txt');  
    boundary_index = load('tool_2_boundary_index.txt'); 
elseif cs == 14     %    ���˻�  m = 4000
    point_all = load('wurenji_2.txt');  
    boundary_index = load('wurenji_2_boundary_index.txt'); 
elseif cs == 15     %    ������  m = 2��
    point_all = load('jiaolianji_2.txt');  
    boundary_index = load('jiaolianji_2_boundary_index.txt'); 
elseif cs == 16     %    ����  m = 6000
    point_all = load('haidun.txt');  
    boundary_index = load('haidun_boundary_inex.txt'); 
elseif cs == 17     %    ���׷�����  m = 6700
    point_all = load('jianyifadongji.txt');  
    boundary_index = load('jianyifadongji_bounary_index.txt'); 
elseif cs == 18
    point_all = load('heqianting_80000.txt');  
    boundary_index = load('heqianting_80000_bounary_index.txt'); 
elseif cs == 19
    point_all = load('heqianting_40000.txt');  
    boundary_index = load('heqianting_40000_boundary_index.txt'); 
elseif cs == 20
    point_all = load('luntai.txt');  
    boundary_index = load('luntai_boundary.txt'); 
elseif cs == 22
    point_all = load('chilun_2.txt');  
    boundary_index = load('chilun_boundary_2.txt'); 
elseif cs == 222
    point_all = load('chilun_3.txt');  
    boundary_index = load('chilun_3_boundary.txt'); 
elseif cs == 33
    point_all = load('fengshan.txt');  
    boundary_index = load('fengshan_boundary.txt');
elseif cs == 333
    point_all = load('fengshan_2.txt');  
    boundary_index = load('fengshan_boundary_2.txt');
elseif cs == 123
    point_all = load('fengshan_6ye.txt');  
    boundary_index = load('fengshan_6ye_boundary.txt');
elseif cs == 55
    point_all = load('hand.txt');  
    boundary_index = load('hand_boundary.txt');
elseif cs == 555
    point_all = load('hand_2.txt');  
    boundary_index = load('hand_2_boundary.txt');
elseif cs == 66
    point_all = load('wolun.txt');  
    boundary_index = load('wolun_boundary.txt');
elseif cs == 77
    point_all = load('fengche.txt');  
    boundary_index = load('fengche_boundary.txt');
elseif cs == 88
    point_all = load('whole_body.txt');  
    boundary_index = load('whole_body_boundary.txt');
elseif cs == 99
    point_all = load('zhandouji.txt');  
    boundary_index = load('zhandouji_boundary.txt');
elseif cs == 999
    point_all = load('zhandouji_2.txt');  
    boundary_index = load('zhandouji_2_boundary.txt');
elseif cs == 9999
    point_all = load('zhandouji_3.txt');  
    boundary_index = load('zhandouji_3_boundary.txt');
elseif cs == 51
    point_all = load('boat.txt');  
    boundary_index = load('boat_boundary.txt');
elseif cs == 551
    point_all = load('boat_2.txt');  
    boundary_index = load('boat_2_boundary.txt');
elseif cs == 71
    point_all = load('fengwo.txt');  
    boundary_index = load('fengwo_boundary.txt');
    
elseif cs == 10001
    point_all = load('CL201_1_hinge.txt');  
    boundary_index = load('CL201_1_hinge_boundary.txt');
elseif cs == 10002
    point_all = load('bird.txt');  
    boundary_index = load('bird_boundary.txt');
elseif cs == 10003
    point_all = load('BlocMotor.txt');  
    boundary_index = load('BlocMotor_boundary.txt');
elseif cs == 10004
    point_all = load('HelicalGear.txt');  
    boundary_index = load('HelicalGear_boundary.txt');
elseif cs == 10005
    point_all = load('sprocket.txt');  
    boundary_index = load('sprocket_boundary.txt');
elseif cs == 100061
    point_all = load('room25600.txt');  
    boundary_index = load('room25600_boundary.txt');
elseif cs == 100062
    point_all = load('room109014.txt');  
    boundary_index = load('room109014_boundary.txt');
elseif cs == 10007
    point_all = load('car4458.txt');  
    boundary_index = load('car4458_boundary.txt');    
elseif cs == 100071
    point_all = load('car8003.txt');  
    boundary_index = load('car8003_boundary.txt');
elseif cs == 100072
    point_all = load('car11419.txt');  
    boundary_index = load('car11419_boundary.txt');
elseif cs == 100073
    point_all = load('car29780.txt');  
    boundary_index = load('car29780_boundary.txt');
    elseif cs == 100074
    point_all = load('car59628.txt');  
    boundary_index = load('car59628_boundary.txt');
end

end