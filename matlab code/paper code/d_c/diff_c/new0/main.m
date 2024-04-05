clc
close all
clear
nub_a=1.5;
nub_b=0.1;
Z=[1,2,3,4,5];
% Z=[1,2,3];
aqa=1000*(nub_a+nub_b*Z);
for CC=1:4
    CE=nub_a+nub_b*CC;
    EER(CC)=EE(CC,CE);
    no_pR(CC)=no_p(CC,CE);
    no_fR(CC)=no_f(CC,CE);
    without_JR(CC)=without_J(CC,CE);
end    
data = [EER;without_JR;no_fR;no_pR];
data = [EER;no_fR;no_pR];
b = bar(data');
ch = get(b,'children');
% set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3)})
set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4),aqa(5)})
% set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4)})
legend('JEE','without_JR','no_f','no_p');
xlabel('Task Load£¨Megzcycles£©')
ylabel('Utility');
% 'distributed'  ;distributed