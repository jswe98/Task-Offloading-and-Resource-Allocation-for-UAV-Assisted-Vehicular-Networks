HU=[-125,200;170,200];%合法用户位置
HE=[-100,50;-350,-380];%窃听者估计位置
HG=[0;0];
H_NFZ=[-100,50;-200,-250];%飞行禁区中心位置


r=45;%飞行禁区1
a=-100;
b=-200;

r1=30;%飞行禁区2
a1=50;
b1=-250;

r2=10;%窃听者1
a2=-100;
b2=-350;

r3=40;%窃听者2
a3=50;
b3=-380;

para=[a-r,b-r,2*r,2*r];
para1=[a1-r1,b1-r1,2*r1,2*r1];


HJ_zuiyou_70=xlsread('zuiyou_T70.xlsx');
% HJ_zuiyou_46=xlsread('zuiyou_T46.xlsx');
HJ_zuiyou_40=xlsread('zuiyou_T40.xlsx');

HJ_wulubang_70=xlsread('wulubang_T70.xlsx');
% HJ_wulubang_46=xlsread('wulubang_T46.xlsx');
HJ_wulubang_40=xlsread('wulubang_T40.xlsx');

HJ_gudingPJ_70=xlsread('gudinPJ_T70.xlsx');
% HJ_gudingPJ_46=xlsread('gudinPJ_T46.xlsx');
HJ_gudingPJ_40=xlsread('gudinPJ_T40.xlsx');

HJ_gudingguiji_70=xlsread('gudinguiji_T70.xlsx');
% HJ_gudingguiji_46=xlsread('gudinguiji_T46.xlsx');
HJ_gudingguiji_40=xlsread('gudinguiji_T40.xlsx');


figure(1)
plot(HJ_zuiyou_70(1,:),HJ_zuiyou_70(2,:),'-rd', ...%绘制最优的轨迹_70
'LineWidth',1, ...
'MarkerSize',4, ...
'MarkerEdgeColor','r')
hold on
% plot(HJ_zuiyou_50(1,:),HJ_zuiyou_50(2,:),'-ro', ...%绘制最优的轨迹_50
% 'LineWidth',1, ...
% 'MarkerSize',7, ...
% 'MarkerEdgeColor','b')
% hold on
% plot(HJ_zuiyou_40(1,:),HJ_zuiyou_40(2,:),'-yo', ...%绘制最优的轨迹_9=40
% 'LineWidth',1, ...
% 'MarkerSize',7, ...
% 'MarkerEdgeColor','b')
% hold on


plot(HJ_wulubang_70(1,:),HJ_wulubang_70(2,:),'-bs', ...%绘制无鲁棒的轨迹_70
'LineWidth',1, ...
'MarkerSize',4, ...
'MarkerEdgeColor','b')
hold on
% plot(HJ_wulubang_50(1,:),HJ_wulubang_50(2,:),'-r^', ...%绘制无鲁棒的轨迹_50
% 'LineWidth',1, ...
% 'MarkerSize',7, ...
% 'MarkerEdgeColor','b')
% hold on
% plot(HJ_wulubang_40(1,:),HJ_wulubang_40(2,:),'-y^', ...%绘制无鲁棒的轨迹_40
% 'LineWidth',1, ...
% 'MarkerSize',7, ...
% 'MarkerEdgeColor','b')
% hold on


plot(HJ_gudingPJ_70(1,:),HJ_gudingPJ_70(2,:),'-go', ...%绘制固定PJ的轨迹_70
'LineWidth',1, ...
'MarkerSize',4, ...
'MarkerEdgeColor','g')
hold on
% plot(HJ_gudingPJ_46(1,:),HJ_gudingPJ_46(2,:),'-rv', ...%绘制固定PJ的轨迹_50
% 'LineWidth',1, ...
% 'MarkerSize',7, ...
% 'MarkerEdgeColor','b')
% hold on
% plot(HJ_gudingPJ_40(1,:),HJ_gudingPJ_40(2,:),'-yv', ...%绘制固定PJ的轨迹_40
% 'LineWidth',1, ...
% 'MarkerSize',7, ...
% 'MarkerEdgeColor','b')
% hold on


plot(HJ_gudingguiji_70(1,:),HJ_gudingguiji_70(2,:),'-m^', ...%绘制固定轨迹的轨迹_70
'LineWidth',1, ...
'MarkerSize',4, ...
'MarkerEdgeColor','m')
hold on
% plot(HJ_gudingguiji_46(1,:),HJ_gudingguiji_46(2,:),'-rs', ...%绘制固定轨迹的轨迹_50
% 'LineWidth',1, ...
% 'MarkerSize',7, ...
% 'MarkerEdgeColor','b')
% hold on
% plot(HJ_gudingguiji_40(1,:),HJ_gudingguiji_40(2,:),'-ys', ...%绘制固定轨迹的轨迹_40
% 'LineWidth',1, ...
% 'MarkerSize',7, ...
% 'MarkerEdgeColor','b')
% hold on


HU1=[HU(1,1);HU(2,1)];
plot(HU1(1,1),HU1(2,1),'o', ...%用户一（红色o）
'MarkerSize',7, ...
'MarkerEdgeColor','g', ...
'MarkerFaceColor','g')
hold on

HE1=[HE(1,1);HE(2,1)];
plot(HE1(1,1),HE1(2,1),'v', ...%窃听者一（黑色下三角）
'MarkerSize',7, ...
'MarkerEdgeColor','k', ...
'MarkerFaceColor','k')
hold on

plot(HG(1,1),HG(2,1),'h', ...%基站（绿色六角星）
'MarkerSize',7, ...
'MarkerEdgeColor','r', ...
'MarkerFaceColor','r')
hold on

% H_NFZ1=[H_NFZ(1,1);H_NFZ(2,1)];
% plot(H_NFZ1(1,1),H_NFZ1(2,1),'x', ...%禁飞区1（红色上三角）
% 'MarkerSize',7, ...
% 'MarkerEdgeColor','r', ...
% 'MarkerFaceColor','r')
% hold on

HU2=[HU(1,2);HU(2,2)];
plot(HU2(1,1),HU2(2,1),'o', ...%用户二（红色o）
'MarkerSize',7, ...
'MarkerEdgeColor','g', ...
'MarkerFaceColor','g')
hold on

HE2=[HE(1,2);HE(2,2)];
plot(HE2(1,1),HE2(2,1),'v', ...%窃听者二（黑色下三角）
'MarkerSize',7, ...
'MarkerEdgeColor','k', ...
'MarkerFaceColor','k')
hold on 

rectangle('Position',para,'Curvature',[1,1],'EdgeColor','b');
rectangle('Position',para1,'Curvature',[1,1],'EdgeColor','b');

% axis equal

% hold on
% H_NFZ2=[H_NFZ(1,2);H_NFZ(2,2)];
% plot(H_NFZ2(1,1),H_NFZ2(2,1),'x', ...%禁飞区2（红色上三角）
% 'MarkerSize',7, ...
% 'MarkerEdgeColor','r', ...
% 'MarkerFaceColor','r')
grid on
xlabel('x(m)')
ylabel('y(m)')
legend('Proposed scheme','Non-robust ','Without-JPC','Fixed trajectory')