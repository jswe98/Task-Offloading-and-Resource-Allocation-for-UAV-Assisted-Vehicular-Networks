 % SPC=[15.5130471523253,15.4523882753053,15.2280919958316,15.5532755616198,15.4413593815687,15.5064568045820,15.7004673998357,15.9486259205149,15.2438250888021];
 % MPC=[41.4498312571311,41.5470185714044,41.8126706181199,42.0205352925397,41.4498312571311,41.5470185714044,41.8126706181199,42.0205352925397,41.4498312571311];
 
 clc
close all
clear
%-------------------ϵͳ����-------------------%
 SPC=[15.5130471523253;15.4523882753053;15.2280919958316;15.5532755616198;15.4413593815687;15.5064568045820;15.7004673998357;15.9486259205149;15.2438250888021];
 MPC=[8.47464131510829;13.6780249244839;9.41941613942833;12.8627178239961;11.8989845788862;15.2864490053887;8.21797450863422;15.5000079158587;12.3883983957634];
 %FEM=[15.5130471523253;15.4523882753053;15.2280919958316;15.5532755616198;15.4413593815687;15.5064568045820;15.7004673998357;15.9486259205149;15.2438250888021];


I=linspace(1,9,9);
x_total=[SPC,MPC];
bar(I,x_total)
ylim([5 20]);
xlim([0 10])  

legend('Robust game','D.C. programming')
xlabel('VUEs','FontSize',12);
ylabel('Transmission rate (Mb/s)','FontSize',12);
   set(gca,'FontSize',12);