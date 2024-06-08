clc
close all
clear
EE=EE();
without_rubst=without_J();
no_f=no_f();
no_p=no_p();
EE(1,1)=4;
without_rubst(1,1)=4;
no_f(1,1)=4;
no_p(1,1)=4;
figure
plot((EE(:,1)),'-+r','linewidth',2);
xlabel('iteration');
ylabel('EE');
hold on
grid on
plot((without_rubst(:,1)),'-+b','linewidth',2);
xlabel('iteration');
ylabel('EE');
plot((no_f(:,1)),'-+c','linewidth',2);
xlabel('iteration');
ylabel('EE');
plot((no_p(:,1)),'-+g','linewidth',2);
legend('Algorithm 1','IOP','Without-CRA','Without-VPC')
xlabel('\fontname{宋体}迭代次数');
% set(gca,'ylim',[0,12]); 
set(gca,'xlim',[1,25]); 
ylabel('\fontname{宋体}系统效用');
 set(gca,'FontName','Times New Roman')