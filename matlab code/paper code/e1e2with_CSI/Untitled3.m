clc
close all
clear
diffE_same_V=diffE_same_V( );
diffE_diff_V=diffE_diff_V( );
figure
plot((diffE_same_V(1,:)),'-+r','linewidth',2);
xlabel('iteration');
ylabel('EE');
hold on
plot((diffE_diff_V(1,:)),'-+b','linewidth',2);
xlabel('iteration');
ylabel('EE');
legend('same_V','diff_V')
set(gca,'XTickLabel',{'0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5'})
% set(gca,'XTickLabel',{'0.1','0.2','0.3','0.4','0.5'})