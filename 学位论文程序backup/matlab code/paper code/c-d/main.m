clc
close all
clear
CEE=2.5;
for FF=1:5 
dddd=0.35+0.05*FF;
aqa(FF)=1000*dddd;
EER(FF)=EE(FF,dddd,CEE);
no_pR(FF)=no_p(FF,dddd,CEE);
no_fR(FF)=no_f(FF,dddd,CEE);
% % % without_rubstR(FF)=without_rubst(FF,dddd);
without_JR(FF)=without_J(FF,dddd,CEE);
end
% data = [EER; no_pR;no_fR;without_rubstR;without_JR] ;
data = [EER; no_pR;no_fR;without_JR] ;
% data = [EER; no_pR;no_fR] ;
b = bar(data');
ch = get(b,'children');
set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4),aqa(5),})
legend('EE','no_p','no_f','without_JR]');
% legend('EE','no_p','no_f');
xlabel('Task Load��Megzcycles��')
ylabel('Utility');