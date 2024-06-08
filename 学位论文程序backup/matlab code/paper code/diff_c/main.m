clc
close all
clear
nub_a=1.1;
nub_b=0.4;
Z=[1,2,3];
aqa=1000*(nub_a+nub_b*Z);
no_p=no_p_diff_c(nub_a,nub_b);
no_f=no_f_diff_c(nub_a,nub_b);
EE=EEdiff_c(nub_a,nub_b);%这个其实是Distributed
JEE=JEEdiff_c(nub_a,nub_b);%这个是联合
% distributed=distributed(nub_a,nub_b);  JEEdiff_c
data = [JEE;EE;no_f;no_p];
b = bar(data');
ch = get(b,'children');
set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3)})
% set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4),aqa(5)})
legend('JEE','EE','no_p','no_f');
xlabel('Task Load（Megzcycles）')
ylabel('Utility');
% 'distributed'  ;distributed