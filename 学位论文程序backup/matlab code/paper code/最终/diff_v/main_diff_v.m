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
xlabel('V (m/s)','FontName','Times New Roman');
ylabel('Utility value','FontName','Times New Roman');
%legend('U')
%legend('U','Imperfect CSI')
set(gca,'ylim',[12.2,13]); 
set(gca,'XTickLabel',{'0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45','0.5'})
set(gca,'XTickLabel',{'20','25','30','35','40','45','50','55','60'})