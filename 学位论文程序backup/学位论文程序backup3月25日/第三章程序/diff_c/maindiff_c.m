clc
close all
clear
% nub_a=1.4;
% nub_b=0.1;
nub_a=1.4;
nub_b=0.1;
Z=[1,2,3,4,5];
ddd=0.42;
% Z=[1,2,3];
aqa=1000*(nub_a+nub_b*Z);
for CC=1:5
    CE=nub_a+nub_b*CC;
    EER(CC)=EE(CC,CE,ddd);
    no_pR(CC)=no_p(CC,CE,ddd);
    no_fR(CC)=no_f(CC,CE,ddd);
    without_JR(CC)=without_J(CC,CE,ddd);
end    

data = [EER;without_JR;no_fR;no_pR];
% data = [EER;no_fR;no_pR];
b = bar(data');
grid on
% set(b(1),'FaceColor',[224,84,28]/255)     
% set(b(2),'FaceColor',[162,26,84]/255)    
% set(b(3),'FaceColor',[240,180,28]/255)    
% set(b(4),'FaceColor',[50,24,60]/255)   
set(b(1),'FaceColor',[224,84,28]/255)     
set(b(2),'FaceColor',[162,26,84]/255)    
set(b(3),'FaceColor',[240,180,28]/255)    
set(b(4),'FaceColor',[50,24,60]/255)   

ch = get(b,'children');
% set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3)})
set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4),aqa(5)})
% set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4)})
legend('JEE','no_f','no_p');
legend('JEE','without_JR','no_f','no_p');
legend('Algorithm 1','IOP','Without-CRA','Without-VPC');
xlabel('任务量（Megzcycles）')
% xlabel('\epsilon1','fontsize',13);
ylabel('系统效用');
% 'distributed'  ;distributed
% 
% plot(EER,'-+r');
% hold on
% plot(without_JR,'-og');
% plot(no_fR,'-*b');
% plot(no_pR,'-sk');
% legend('JEE','without_JR','no_f','no_p');