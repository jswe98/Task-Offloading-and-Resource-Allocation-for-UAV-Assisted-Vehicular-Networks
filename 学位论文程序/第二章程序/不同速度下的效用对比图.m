 clc
close all
clear
%-------------------ϵͳ����-------------------%
 U_sum=[62.4410478106308;56.1304258741786;46.8790851952090;41.5126706181199;40.6134283788861];
 U_0=[2.38045445836699;3.39613165162366;5.96451476045873;8.12154902603442;8.56520533176818];
%[0.00624410478106308][0.00561304258741786][0.00468790851952090][0.00415126706181199][0.00406134283788861]
%[0.000238045445836699][0.000339613165162366][0.000596451476045873][0.000812154902603442][0.000856520533176818]
 
 speed=linspace(0,40,5);

figure 
plot(speed,U_sum,'-og','linewidth',2);
hold on
plot(speed,U_0,'-or','linewidth',2);

grid on
      
legend('U_sum','U_0')'%I_\infty-approximate'
xlim([0 40]);
ylim([0 70]);

xlabel('V (m/s)','FontSize',12);
ylabel('Utility value','FontSize',12);
   set(gca,'FontSize',10);