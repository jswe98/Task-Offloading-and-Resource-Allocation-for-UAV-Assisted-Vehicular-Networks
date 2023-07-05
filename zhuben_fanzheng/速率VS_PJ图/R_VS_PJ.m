PJ_0_0=xlsread('R_avg_zuiyou_PJ_0.xlsx');
PJ_0_0003=xlsread('R_avg_zuiyou_PJ_0.0003.xlsx');
PJ_0_0006=xlsread('R_avg_zuiyou_PJ_0.0006.xlsx');
PJ_0_0009=xlsread('R_avg_zuiyou_PJ_0.0009.xlsx');
PJ_0_0012=xlsread('R_avg_zuiyou_PJ_0.0012.xlsx');
PJ_0_0015=xlsread('R_avg_zuiyou_PJ_0.0015.xlsx');
PJ_0_0018=xlsread('R_avg_zuiyou_PJ_0.0018.xlsx');
PJ_0_0021=xlsread('R_avg_zuiyou_PJ_0.0021.xlsx');
PJ_0_0024=xlsread('R_avg_zuiyou_PJ_0.0024.xlsx');
PJ_0_0027=xlsread('R_avg_zuiyou_PJ_0.0027.xlsx');

a=[PJ_0_0,PJ_0_0003,PJ_0_0006,PJ_0_0009,PJ_0_0012,PJ_0_0015,PJ_0_0018,PJ_0_0021,PJ_0_0024,PJ_0_0027];
for i=1:10
    b(1,i)=-0.0003+0.0003*i;
end
R_avg_zuiyou_PJ=[b;a];

PJ2_0_0=xlsread('R_avg_wulubang_PJ_0.xlsx');
PJ2_0_0003=xlsread('R_avg_wulubang_PJ_0.0003.xlsx');
PJ2_0_0006=xlsread('R_avg_wulubang_PJ_0.0006.xlsx');
PJ2_0_0009=xlsread('R_avg_wulubang_PJ_0.0009.xlsx');
PJ2_0_0012=xlsread('R_avg_wulubang_PJ_0.0012.xlsx');
PJ2_0_0015=xlsread('R_avg_wulubang_PJ_0.0015.xlsx');
PJ2_0_0018=xlsread('R_avg_wulubang_PJ_0.0018.xlsx');
PJ2_0_0021=xlsread('R_avg_wulubang_PJ_0.0021.xlsx');
PJ2_0_0024=xlsread('R_avg_wulubang_PJ_0.0024.xlsx');
PJ2_0_0027=xlsread('R_avg_wulubang_PJ_0.0027.xlsx');


a2=[PJ2_0_0,PJ2_0_0003,PJ2_0_0006,PJ2_0_0009,PJ2_0_0012,PJ2_0_0015,PJ2_0_0018,PJ2_0_0021,PJ2_0_0024,PJ2_0_0027];
for i=1:10
    b2(1,i)=-0.0003+0.0003*i;
end
R_avg_wulubang_PJ=[b2;a2];

PJ3_0_0=xlsread('gudinPJ_R_avg_PJ_0.xlsx');
PJ3_0_0003=xlsread('gudinPJ_R_avg_PJ_0.0003.xlsx');
PJ3_0_0006=xlsread('gudinPJ_R_avg_PJ_0.0006.xlsx');
PJ3_0_0009=xlsread('gudinPJ_R_avg_PJ_0.0009.xlsx');
PJ3_0_0012=xlsread('gudinPJ_R_avg_PJ_0.0012.xlsx');
PJ3_0_0015=xlsread('gudinPJ_R_avg_PJ_0.0015.xlsx');
PJ3_0_0018=xlsread('gudinPJ_R_avg_PJ_0.0018.xlsx');
PJ3_0_0021=xlsread('gudinPJ_R_avg_PJ_0.0021.xlsx');
PJ3_0_0024=xlsread('gudinPJ_R_avg_PJ_0.0024.xlsx');
PJ3_0_0027=xlsread('gudinPJ_R_avg_PJ_0.0027.xlsx');

a3=[PJ3_0_0,PJ3_0_0003,PJ3_0_0006,PJ3_0_0009,PJ3_0_0012,PJ3_0_0015,PJ3_0_0018,PJ3_0_0021,PJ3_0_0024,PJ3_0_0027];
for i=1:10
    b3(1,i)=-0.0003+0.0003*i;
end
R_avg_gudinPJ_PJ=[b3;a3];

PJ4_0_0=xlsread('gudinguiji_PJ_0.xlsx');
PJ4_0_0003=xlsread('gudinguiji_PJ_0.0003.xlsx');
PJ4_0_0006=xlsread('gudinguiji_PJ_0.0006.xlsx');
PJ4_0_0009=xlsread('gudinguiji_PJ_0.0009.xlsx');
PJ4_0_0012=xlsread('gudinguiji_PJ_0.0012.xlsx');
PJ4_0_0015=xlsread('gudinguiji_PJ_0.0015.xlsx');
PJ4_0_0018=xlsread('gudinguiji_PJ_0.0018.xlsx');
PJ4_0_0021=xlsread('gudinguiji_PJ_0.0021.xlsx');
PJ4_0_0024=xlsread('gudinguiji_PJ_0.0024.xlsx');
PJ4_0_0027=xlsread('gudinguiji_PJ_0.0027.xlsx');

a4=[PJ4_0_0,PJ4_0_0003,PJ4_0_0006,PJ4_0_0009,PJ4_0_0012,PJ4_0_0015,PJ4_0_0018,PJ4_0_0021,PJ4_0_0024,PJ4_0_0027];
for i=1:10
    b4(1,i)=-0.0003+0.0003*i;
end
gudinguiji_PJ=[b4;a4];

figure(1)
plot(R_avg_zuiyou_PJ(1,:),R_avg_zuiyou_PJ(2,:),'-rd', ...%最优方案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','r')
hold on
plot(R_avg_wulubang_PJ(1,:),R_avg_wulubang_PJ(2,:),'-bs', ...%无鲁棒案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','b')
hold on
plot(R_avg_gudinPJ_PJ(1,:),R_avg_gudinPJ_PJ(2,:),'-go', ...%固定PJ案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','g')
hold on
plot(gudinguiji_PJ(1,:),gudinguiji_PJ(2,:),'-m^', ...%固定PJ案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','m')
grid on
xlabel('P_{j}^{avg}(w)')
ylabel('Max-min average secrecy rate (bps/Hz)')
legend('Proposed scheme','Non-robust','Without-JPC','Fixed trajectory')