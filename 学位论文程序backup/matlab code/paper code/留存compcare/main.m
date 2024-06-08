clc
close all
clear
no_p=no_p();
no_f=no_f();
EE=EE();
without_rubst=without_rubst();
figure
plot((EE(:,1)),'-+r');
xlabel('iteration');
ylabel('EE');
hold on
plot((no_f(:,1)),'-+b');
xlabel('iteration');
ylabel('EE');
plot((without_rubst(:,1)),'-+c');
xlabel('iteration');
ylabel('EE');
plot((no_p(:,1)),'-+g');
legend('EE£¨p&f£©','EEno_f','without_rubst','EEno_p')
xlabel('iteration');
set(gca,'ylim',[0,12]); 
ylabel('EE');