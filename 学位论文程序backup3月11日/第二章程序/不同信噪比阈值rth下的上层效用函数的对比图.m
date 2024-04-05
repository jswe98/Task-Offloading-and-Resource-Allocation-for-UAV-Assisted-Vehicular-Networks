clc
close all
clear
%-------------------ϵͳ����-------------------%
 U1=[8.17104095500104;8.53595067369873;8.67844721885958;8.68441889469840;8.60120176904096];
 U2=[8.28199988766811;8.64916704881360;8.79031588670405;8.79291948258303;8.70513158086664];
 U3=[8.39619145707464;8.76557969747731;8.90524207345868;8.90428605553440;8.81170740997560];
 U4=[8.51376029012299;8.88532726951379;9.02335420266773;9.01863479715680;8.92103223206717];
 U5=[8.63485981220045;9.00855660036039;9.14478803562334;9.13608827834433;9.03321441611430];
%[0.00624410478106308][0.00561304258741786][0.00468790851952090][0.00415126706181199][0.00406134283788861]
%[0.000238045445836699][0.000339613165162366][0.000596451476045873][0.000812154902603442][0.000856520533176818]
 
 speed=linspace(0.1,0.3,5);

figure 
plot(speed,U1,'-oc','linewidth',2);
hold on
plot(speed,U2,'-*g','linewidth',2);
plot(speed,U3,'-dk','linewidth',2);
plot(speed,U4,'-sr','linewidth',2);
plot(speed,U5,'-+b','linewidth',2);
grid on
      

      
legend('I_{th}=0.001','I_{th}=0.003','I_{th}=0.005','I_{th}=0.007','I_{th}=0.009')'%I_\infty-approximate'
ylim([8 9.2]);
xlim([0.1 0.3])

xlabel('\gamma_{th}','FontSize',12);
ylabel('U_0','FontSize',12);
   set(gca,'FontSize',10);