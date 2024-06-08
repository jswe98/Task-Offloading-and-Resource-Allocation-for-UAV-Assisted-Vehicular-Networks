clc
close all
clear
for FF=1:5
V1=-10*ones(1,5)+10*FF*ones(1,5);%发送车CM速度
V1=15*ones(1,5)+5*FF*ones(1,5);%发送车CM速度
VV(FF,:)=V1;
end
imperfect=imperfect(V1);
perfect=perfect(V1);
figure

plot((perfect(1,:)),'-+r','linewidth',2);
xlabel('iteration');
ylabel('EE');
grid on
hold on
plot((imperfect(1,:)),'-+r','linewidth',2);
% xlabel('车辆速度 (m/s)','FontName','Times New Roman');
% xlabel('车辆速度 (m/s)', 'Interpreter', 'latex');  % 使用 LaTeX 语法，仅对 (m/s) 部分应用 Times New Roman 字体
% % xlabel('车辆速度 (m/s)');
% ylabel('系统效用');
%legend('U')
%legend('U','Imperfect CSI')
set(gca,'ylim',[12.5,13]); 
set(gca,'XTickLabel',{'0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5'})
set(gca,'XTickLabel',{'20','25','30','35','40','45','50','55','60'})
% set(gca,'YTickLabel',{'12.5','12.6','12.7','12.8','12.9','13'})
set(gca,'FontName','Times New Roman')
xlabel('\fontname{宋体}车辆速度 \fontname{Times New Roman}(m/s)');
ylabel('\fontname{宋体}系统效用 ');
xlabel('\nu(m/s)','FontName','Times New Roman');
ylabel('Utility','FontName','Times New Roman');
%ax = gca;  % 获取当前坐标轴对象
%ax.XAxis.FontName = 'Times New Roman';  % 设置 x 轴标签字体为 Times New Roman