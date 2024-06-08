clc
close all
clear
no_p=no_p();
no_f=no_f();
EE=EE();
without_rubst=without_J();
figure
plot((EE(:,1)),'-+r','linewidth',1);
xlabel('iteration');
ylabel('EE');
hold on
plot((no_f(:,1)),'-+b','linewidth',1);
xlabel('iteration');
ylabel('EE');
plot((without_rubst(:,1)),'-+c','linewidth',1);
xlabel('iteration');
ylabel('EE');
plot((no_p(:,1)),'-+g','linewidth',1);
legend('EE£¨p&f£©','EEno_f','without_speed','EEno_p')
xlabel('iteration');
% set(gca,'ylim',[0,12]); 
% set(gca,'xlim',[0,28]); 
ylabel('EE');