clc
close all
clear
zihao=18;
for FF=1:4 
dddd=25+2*FF;
aqa(FF)=0.1*dddd;
EER(FF)=EE(FF,dddd);
no_pR(FF)=no_p(FF,dddd);
no_fR(FF)=no_f(FF,dddd);
% % % without_rubstR(FF)=without_rubst(FF,dddd);
without_JR(FF)=without_J(FF,dddd);
end
% data = [EER; no_pR;no_fR;without_rubstR;without_JR] ;
data = [EER;without_JR;no_fR; no_pR] ;
% data = [EER; no_pR;no_fR] ;
b = bar(data');

grid on

% set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4),aqa(5)})
set(b(1),'FaceColor',[8,116,188]/255)     
set(b(2),'FaceColor',[224,84,28]/255)    
set(b(3),'FaceColor',[240,180,28]/255)    
set(b(4),'FaceColor',[128,44,140]/255) 

% set(b(1),'FaceColor',[224,84,28]/255)     
% set(b(2),'FaceColor',[162,26,84]/255)    
% set(b(3),'FaceColor',[240,180,28]/255)    
% set(b(4),'FaceColor',[50,24,60]/255)   
ch = get(b,'children');
% set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4),aqa(5)})
set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4)})
legend('EE','without_JR]','no_p','no_f');
legend('Algorithm 1','IOP','Without-CRA','Without-VPC');
% legend('EE','no_p','no_f');
% xlabel('\fontname{宋体}总计算能力的阈值 \fontname{Times New Roman}(GHz)')
xlabel('\fontname{Times New Roman}Threshold value of computing total(GHz)')
% ylabel('\fontname{宋体}系统效用');
ylabel('\fontname{Times New Roman}Utility');
ax = gca;
ax.YLabel.FontSize = zihao;
ax.XLabel.FontSize = zihao;
 set(gca,'FontName','Times New Roman')