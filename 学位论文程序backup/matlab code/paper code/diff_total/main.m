clc
close all
clear
for FF=1:4 
dddd=25+2*FF;
aqa(FF)=1000*dddd;
EER(FF)=EE(FF,dddd);
no_pR(FF)=no_p(FF,dddd);
no_fR(FF)=no_f(FF,dddd);
% % % without_rubstR(FF)=without_rubst(FF,dddd);
without_JR(FF)=without_J(FF,dddd);
end
% data = [EER; no_pR;no_fR;without_rubstR;without_JR] ;
data = [EER;without_JR; no_pR;no_fR] ;
% data = [EER; no_pR;no_fR] ;
b = bar(data');
ch = get(b,'children');
% set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4),aqa(5)})
set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4)})
legend('EE','without_JR]','no_p','no_f');
% legend('EE','no_p','no_f');
xlabel('Task Load£¨Megzcycles£©')
ylabel('Utility');