clc
close all
clear
for FF=1:4 
%%dddd=0.35+0.05*FF; backup
dddd=0.35+0.05*FF;
aqa(FF)=1000*dddd;
EER(FF)=EE(FF,dddd);
no_pR(FF)=no_p(FF,dddd);
no_fR(FF)=no_f(FF,dddd);
% % % without_rubstR(FF)=without_rubst(FF,dddd);
without_JR(FF)=without_J(FF,dddd);
end
% data = [EER; no_pR;no_fR;without_rubstR;without_JR] ;
data = [EER;no_fR;without_JR;no_pR] ;
% data = [EER; no_pR;no_fR] ;
b = bar(data');
grid on
% set(b(1),'FaceColor',[224,84,28]/255)     
% set(b(2),'FaceColor',[162,26,84]/255)    
% set(b(3),'FaceColor',[240,180,28]/255)    
% set(b(4),'FaceColor',[50,24,60]/255)   
set(b(1),'FaceColor',[126,153,244]/255)     
set(b(2),'FaceColor',[204,124,113]/255)    
set(b(3),'FaceColor',[122,182,86]/255)  
set(b(4),'FaceColor',[240,180,28]/255)  

ch = get(b,'children');
% set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4),aqa(5)})
set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4)})
legend('Algorithm 3-1','Without-CRA','IOP','Without-VPC');
% legend('EE','no_p','no_f');
xlabel('\fontname{宋体}任务输入 \fontname{Times New Roman}(kb)')
ylabel('\fontname{宋体}系统效用');
 set(gca,'FontName','Times New Roman')