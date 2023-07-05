iterations_zuiyou=xlsread('iterations_zuiyou.xlsx');
iterations_wulubang=xlsread('iterations_wulubang.xlsx');
iterations_gudinguiji=xlsread('iterations_gudinguiji.xlsx');
iterations_gudingPJ=xlsread('iterations_gudingPJ.xlsx');

figure(1)
plot(iterations_zuiyou(1,:),iterations_zuiyou(2,:),'-rd', ...%最优方案下，
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','r')
hold on
plot(iterations_wulubang(1,:),iterations_wulubang(2,:),'-bs', ...%无鲁棒案下，
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','b')
hold on
plot(iterations_gudinguiji(1,:),iterations_gudinguiji(2,:),'-go', ...%固定PJ案下，
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','g')
hold on
plot(iterations_gudingPJ(1,:),iterations_gudingPJ(2,:),'-m^', ...%固定轨迹案下，
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','m')
grid on
xlabel('The number of iterations')
ylabel('Max-min average secrecy rate (bps/Hz)')
legend('Proposed scheme','Non-robust','Fixed trajectory','Without-JPC')
