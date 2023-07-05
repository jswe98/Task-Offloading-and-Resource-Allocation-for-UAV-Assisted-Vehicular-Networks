%最优方案，平均干扰功率为0.0007下的速率随PG变化
PG_0_0=xlsread('R_avg_zuiyou_PG_0.xlsx');
PG_0_3=xlsread('R_avg_zuiyou_PG_0.3.xlsx');
PG_0_6=xlsread('R_avg_zuiyou_PG_0.6.xlsx');
PG_0_9=xlsread('R_avg_zuiyou_PG_0.9.xlsx');
PG_1_2=xlsread('R_avg_zuiyou_PG_1.2.xlsx');
PG_1_5=xlsread('R_avg_zuiyou_PG_1.5.xlsx');
PG_1_8=xlsread('R_avg_zuiyou_PG_1.8.xlsx');
PG_2_1=xlsread('R_avg_zuiyou_PG_2.1.xlsx');
PG_2_4=xlsread('R_avg_zuiyou_PG_2.4.xlsx');
PG_2_7=xlsread('R_avg_zuiyou_PG_2.7.xlsx');

a=[PG_0_0,PG_0_3,PG_0_6,PG_0_9,PG_1_2,PG_1_5,PG_1_8,PG_2_1,PG_2_4,PG_2_7];
for i=1:10
    b(1,i)=-0.3+0.3*i;
end
R_avg_zuiyou_PG=[b;a];
%最优方案，平均干扰功率为0.0012下的速率随PG变化
JPG_0_0=xlsread('R_avg_zuiyou_PG_0_0.0012.xlsx');
JPG_0_3=xlsread('R_avg_zuiyou_PG_0.3_0.0012.xlsx');
JPG_0_6=xlsread('R_avg_zuiyou_PG_0.6_0.0012.xlsx');
JPG_0_9=xlsread('R_avg_zuiyou_PG_0.9_0.0012.xlsx');
JPG_1_2=xlsread('R_avg_zuiyou_PG_1.2_0.0012.xlsx');
JPG_1_5=xlsread('R_avg_zuiyou_PG_1.5_0.0012.xlsx');
JPG_1_8=xlsread('R_avg_zuiyou_PG_1.8_0.0012.xlsx');
JPG_2_1=xlsread('R_avg_zuiyou_PG_2.1_0.0012.xlsx');
JPG_2_4=xlsread('R_avg_zuiyou_PG_2.4_0.0012.xlsx');
JPG_2_7=xlsread('R_avg_zuiyou_PG_2.7_0.0012.xlsx');

Ja=[JPG_0_0,JPG_0_3,JPG_0_6,JPG_0_9,JPG_1_2,JPG_1_5,JPG_1_8,JPG_2_1,JPG_2_4,JPG_2_7];
for i=1:10
    Jb(1,i)=-0.3+0.3*i;
end
JR_avg_zuiyou_PG=[Jb;Ja];

%无鲁棒，平均干扰功率为0.0007下的速率随PG变化
PG2_0_0=xlsread('R_avg_wulubang_PG_0.xlsx');
PG2_0_3=xlsread('R_avg_wulubang_PG_0.3.xlsx');
PG2_0_6=xlsread('R_avg_wulubang_PG_0.6.xlsx');
PG2_0_9=xlsread('R_avg_wulubang_PG_0.9.xlsx');
PG2_1_2=xlsread('R_avg_wulubang_PG_1.2.xlsx');
PG2_1_5=xlsread('R_avg_wulubang_PG_1.5.xlsx');
PG2_1_8=xlsread('R_avg_wulubang_PG_1.8.xlsx');
PG2_2_1=xlsread('R_avg_wulubang_PG_2.1.xlsx');
PG2_2_4=xlsread('R_avg_wulubang_PG_2.4.xlsx');
PG2_2_7=xlsread('R_avg_wulubang_PG_2.7.xlsx');

a2=[PG2_0_0,PG2_0_3,PG2_0_6,PG2_0_9,PG2_1_2,PG2_1_5,PG2_1_8,PG2_2_1,PG2_2_4,PG2_2_7];
for i=1:10
    b2(1,i)=-0.3+0.3*i;
end
R_avg_wulubang_PG=[b2;a2];

%无鲁棒，平均干扰功率为0.0012下的速率随PG变化
JPG2_0_0=xlsread('R_avg_wulubang_PG_0_0.0012.xlsx');
JPG2_0_3=xlsread('R_avg_wulubang_PG_0.3_0.0012.xlsx');
JPG2_0_6=xlsread('R_avg_wulubang_PG_0.6_0.0012.xlsx');
JPG2_0_9=xlsread('R_avg_wulubang_PG_0.9_0.0012.xlsx');
JPG2_1_2=xlsread('R_avg_wulubang_PG_1.2_0.0012.xlsx');
JPG2_1_5=xlsread('R_avg_wulubang_PG_1.5_0.0012.xlsx');
JPG2_1_8=xlsread('R_avg_wulubang_PG_1.8_0.0012.xlsx');
JPG2_2_1=xlsread('R_avg_wulubang_PG_2.1_0.0012.xlsx');
JPG2_2_4=xlsread('R_avg_wulubang_PG_2.4_0.0012.xlsx');
JPG2_2_7=xlsread('R_avg_wulubang_PG_2.7_0.0012.xlsx');

Ja2=[JPG2_0_0,JPG2_0_3,JPG2_0_6,JPG2_0_9,JPG2_1_2,JPG2_1_5,JPG2_1_8,JPG2_2_1,JPG2_2_4,JPG2_2_7];
for i=1:10
    Jb2(1,i)=-0.3+0.3*i;
end
JR_avg_wulubang_PG=[Jb2;Ja2];

%固定干扰机功率，平均干扰机功率为0.0007下的速率随PG变化
PG3_0_0=xlsread('gudinPJ_R_avg_PG_0.xlsx');
PG3_0_3=xlsread('gudinPJ_R_avg_PG_0.3.xlsx');
PG3_0_6=xlsread('gudinPJ_R_avg_PG_0.6.xlsx');
PG3_0_9=xlsread('gudinPJ_R_avg_PG_0.9.xlsx');
PG3_1_2=xlsread('gudinPJ_R_avg_PG_1.2.xlsx');
PG3_1_5=xlsread('gudinPJ_R_avg_PG_1.5.xlsx');
PG3_1_8=xlsread('gudinPJ_R_avg_PG_1.8.xlsx');
PG3_2_1=xlsread('gudinPJ_R_avg_PG_2.1.xlsx');
PG3_2_4=xlsread('gudinPJ_R_avg_PG_2.4.xlsx');
PG3_2_7=xlsread('gudinPJ_R_avg_PG_2.7.xlsx');

a3=[PG3_0_0,PG3_0_3,PG3_0_6,PG3_0_9,PG3_1_2,PG3_1_5,PG3_1_8,PG3_2_1,PG3_2_4,PG3_2_7];
for i=1:10
    b3(1,i)=-0.3+0.3*i;
end
R_avg_gudinPJ_PG=[b3;a3];
%固定干扰机功率，平均干扰机功率为0.0012下的速率随PG变化
JPG3_0_0=xlsread('gudinPJ_R_avg_PG_0_0.0012.xlsx');
JPG3_0_3=xlsread('gudinPJ_R_avg_PG_0.3_0.0012.xlsx');
JPG3_0_6=xlsread('gudinPJ_R_avg_PG_0.6_0.0012.xlsx');
JPG3_0_9=xlsread('gudinPJ_R_avg_PG_0.9_0.0012.xlsx');
JPG3_1_2=xlsread('gudinPJ_R_avg_PG_1.2_0.0012.xlsx');
JPG3_1_5=xlsread('gudinPJ_R_avg_PG_1.5_0.0012.xlsx');
JPG3_1_8=xlsread('gudinPJ_R_avg_PG_1.8_0.0012.xlsx');
JPG3_2_1=xlsread('gudinPJ_R_avg_PG_2.1_0.0012.xlsx');
JPG3_2_4=xlsread('gudinPJ_R_avg_PG_2.4_0.0012.xlsx');
JPG3_2_7=xlsread('gudinPJ_R_avg_PG_2.7_0.0012.xlsx');

Ja3=[JPG3_0_0,JPG3_0_3,JPG3_0_6,JPG3_0_9,JPG3_1_2,JPG3_1_5,JPG3_1_8,JPG3_2_1,JPG3_2_4,JPG3_2_7];
for i=1:10
    Jb3(1,i)=-0.3+0.3*i;
end
JR_avg_gudinPJ_PG=[Jb3;Ja3];
%固定轨迹，平均干扰机功率为0.0007下的速率随PG变化
PG4_0_0=xlsread('gudinguiji_PG_0.xlsx');
PG4_0_3=xlsread('gudinguiji_PG_0.3.xlsx');
PG4_0_6=xlsread('gudinguiji_PG_0.6.xlsx');
PG4_0_9=xlsread('gudinguiji_PG_0.9.xlsx');
PG4_1_2=xlsread('gudinguiji_PG_1.2.xlsx');
PG4_1_5=xlsread('gudinguiji_PG_1.5.xlsx');
PG4_1_8=xlsread('gudinguiji_PG_1.8.xlsx');
PG4_2_1=xlsread('gudinguiji_PG_2.1.xlsx');
PG4_2_4=xlsread('gudinguiji_PG_2.4.xlsx');
PG4_2_7=xlsread('gudinguiji_PG_2.7.xlsx');

a4=[PG4_0_0,PG4_0_3,PG4_0_6,PG4_0_9,PG4_1_2,PG4_1_5,PG4_1_8,PG4_2_1,PG4_2_4,PG4_2_7];
for i=1:10
    b4(1,i)=-0.3+0.3*i;
end
R_avg_gudinguiji_PG=[b4;a4];
%固定轨迹，平均干扰机功率为0.0012下的速率随PG变化
JPG4_0_0=xlsread('gudinguiji_PG_0_0.0012.xlsx');
JPG4_0_3=xlsread('gudinguiji_PG_0.3_0.0012.xlsx');
JPG4_0_6=xlsread('gudinguiji_PG_0.6_0.0012.xlsx');
JPG4_0_9=xlsread('gudinguiji_PG_0.9_0.0012.xlsx');
JPG4_1_2=xlsread('gudinguiji_PG_1.2_0.0012.xlsx');
JPG4_1_5=xlsread('gudinguiji_PG_1.5_0.0012.xlsx');
JPG4_1_8=xlsread('gudinguiji_PG_1.8_0.0012.xlsx');
JPG4_2_1=xlsread('gudinguiji_PG_2.1_0.0012.xlsx');
JPG4_2_4=xlsread('gudinguiji_PG_2.4_0.0012.xlsx');
JPG4_2_7=xlsread('gudinguiji_PG_2.7_0.0012.xlsx');

Ja4=[JPG4_0_0,JPG4_0_3,JPG4_0_6,JPG4_0_9,JPG4_1_2,JPG4_1_5,JPG4_1_8,JPG4_2_1,JPG4_2_4,JPG4_2_7];
for i=1:10
    Jb4(1,i)=-0.3+0.3*i;
end
JR_avg_gudinguiji_PG=[Jb4;Ja4];

figure(1)
%在干扰机平均功率为0.0007下的各种方案的速率随PG变化
plot(R_avg_zuiyou_PG(1,:),R_avg_zuiyou_PG(2,:),'-rp', ...%最优方案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','r', ...
'MarkerFaceColor','r')
hold on
plot(R_avg_wulubang_PG(1,:),R_avg_wulubang_PG(2,:),'-bs', ...%无鲁棒案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','b', ...
'MarkerFaceColor','b')
hold on
plot(R_avg_gudinPJ_PG(1,:),R_avg_gudinPJ_PG(2,:),'-go', ...%固定PJ案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','g', ...
'MarkerFaceColor','g')
hold on
plot(R_avg_gudinguiji_PG(1,:),R_avg_gudinguiji_PG(2,:),'-m^', ...%固定轨迹案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','m', ...
'MarkerFaceColor','m')
%在干扰机平均功率为0.0007下的各种方案的速率随PG变化
plot(JR_avg_zuiyou_PG(1,:),JR_avg_zuiyou_PG(2,:),':rp', ...%最优方案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','r', ...
'MarkerFaceColor','r')
hold on
plot(JR_avg_wulubang_PG(1,:),JR_avg_wulubang_PG(2,:),':bs', ...%无鲁棒案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','b', ...
'MarkerFaceColor','b')
hold on
plot(JR_avg_gudinPJ_PG(1,:),JR_avg_gudinPJ_PG(2,:),':go', ...%固定PJ案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','g', ...
'MarkerFaceColor','g')
hold on
plot(JR_avg_gudinguiji_PG(1,:),JR_avg_gudinguiji_PG(2,:),':m^', ...%固定轨迹案下，速率随着PG变化
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','m', ...
'MarkerFaceColor','m')
grid on
xlabel('P_{g}^{avg}(w)')
ylabel('Max-min average  secrecy rate (bps/Hz)')
legend('Proposed scheme,','Non-robust','Without-JPC','Fixed trajectory')
