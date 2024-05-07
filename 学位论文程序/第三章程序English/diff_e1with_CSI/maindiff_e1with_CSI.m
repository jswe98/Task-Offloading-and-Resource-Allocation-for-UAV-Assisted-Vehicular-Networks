clc
close all
clear
zihao=18;
diffE_same_V=diffE_same_V( );
diffE_diff_V=diffE_diff_V( );
figure
plot((diffE_same_V(1,:)),'-+r','linewidth',2);
xlabel('iteration');
ylabel('EE');
hold on
grid on
plot((diffE_diff_V(1,:)),'-+b','linewidth',2);
% xlabel('��1');
% ylabel('ϵͳЧ��');
% legend('\fontname{����}��ͬ�ٶ�','\fontname{����}��ͬ�ٶ�')
legend('\fontname{Times New Roman}Same Speed','\fontname{Times New Roman}Different Speed')
set(gca,'FontName','Times New Roman')
xlabel('$\varepsilon_1$ ', 'Interpreter', 'latex'); 
set(gca,'ylim',[12.45,12.53]); 
% ylabel('\fontname{����}ϵͳЧ�� ');
ylabel('\fontname{Times New Roman}Utility');
ax = gca;
ax.YLabel.FontSize = zihao;
ax.XLabel.FontSize = zihao;
set(gca,'XTickLabel',{'0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5'})
% set(gca,'XTickLabel',{'0.1','0.2','0.3','0.4','0.5'})