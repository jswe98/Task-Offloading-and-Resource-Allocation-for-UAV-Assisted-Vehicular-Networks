clc
close all
clear
nub_a=1;
nub_b=0.4;
Z=[1,2,3];
aqa=1000*(nub_a+nub_b*Z);
no_p=no_p_diff_c(nub_a,nub_b);
no_f=no_f_diff_c(nub_a,nub_b);
EE=EEdiff_c(nub_a,nub_b);
data = [EE; no_f;no_p] ;
b = bar(data');
ch = get(b,'children');
set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3)})
legend('EE','no_p','no_f');
ylabel('Utility');