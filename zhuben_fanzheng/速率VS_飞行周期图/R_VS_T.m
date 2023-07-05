T_40=xlsread('R_avg_zuiyou_40.xlsx');
T_50=xlsread('R_avg_zuiyou_50.xlsx');
T_60=xlsread('R_avg_zuiyou_60.xlsx');
T_70=xlsread('R_avg_zuiyou_70.xlsx');
T_80=xlsread('R_avg_zuiyou_80.xlsx');
T_90=xlsread('R_avg_zuiyou_90.xlsx');


a=[T_40,T_50,T_60,T_70,T_80,T_90];
b=[];
for i=1:6
    b(1,i)=30+10*i;
end
R_avg_zuiyou_T=[b;a];

T2_40=xlsread('R_avg_wulubang_40.xlsx');
T2_50=xlsread('R_avg_wulubang_50.xlsx');
T2_60=xlsread('R_avg_wulubang_60.xlsx');
T2_70=xlsread('R_avg_wulubang_70.xlsx');
T2_80=xlsread('R_avg_wulubang_80.xlsx');
T2_90=xlsread('R_avg_wulubang_90.xlsx');


a2=[T2_40,T2_50,T2_60,T2_70,T2_80,T2_90];
b2=[];
for i=1:6
    b2(1,i)=30+10*i;
end
R_avg_wulubang_T=[b2;a2];

T3_40=xlsread('gudinPJ_R_avg_40.xlsx');
T3_50=xlsread('gudinPJ_R_avg_50.xlsx');
T3_60=xlsread('gudinPJ_R_avg_60.xlsx');
T3_70=xlsread('gudinPJ_R_avg_70.xlsx');
T3_80=xlsread('gudinPJ_R_avg_80.xlsx');
T3_90=xlsread('gudinPJ_R_avg_90.xlsx');


a3=[T3_40,T3_50,T3_60,T3_70,T3_80,T3_90];
b3=[];
for i=1:6
    b3(1,i)=30+10*i;
end
R_avg_gudinPJ_T=[b3;a3];

T4_40=xlsread('gudinguiji_40.xlsx');
T4_50=xlsread('gudinguiji_50.xlsx');
T4_60=xlsread('gudinguiji_60.xlsx');
T4_70=xlsread('gudinguiji_70.xlsx');
T4_80=xlsread('gudinguiji_80.xlsx');
T4_90=xlsread('gudinguiji_90.xlsx');


a4=[T4_40,T4_50,T4_60,T4_70,T4_80,T4_90];
b4=[];
for i=1:6
    b4(1,i)=30+10*i;
end
R_avg_gudinguiji_T=[b4;a4];

figure(1)
plot(R_avg_zuiyou_T(1,:),R_avg_zuiyou_T(2,:),'-rd', ...%最优方案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','r')
hold on
plot(R_avg_wulubang_T(1,:),R_avg_wulubang_T(2,:),'-bs', ...%无鲁棒案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','b')
hold on
plot(R_avg_gudinPJ_T(1,:),R_avg_gudinPJ_T(2,:),'-go', ...%固定PJ案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','g')
hold on
plot(R_avg_gudinguiji_T(1,:),R_avg_gudinguiji_T(2,:),'-m^', ...%固定轨迹案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','m')
grid on
xlabel('T(s)')
ylabel('Max-min average secrecy rate (bps/Hz)')
legend('Proposed scheme','No-robust','Without-PJC','Fixed trajectory')