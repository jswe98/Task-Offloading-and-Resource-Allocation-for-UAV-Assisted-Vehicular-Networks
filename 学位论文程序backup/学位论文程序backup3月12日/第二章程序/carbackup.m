clc
close all
clear
digits(3);
EEE=1;
programmatic=4;
for EEm = 1:4
%   M=  EEm;

    switch EEm 
    case 1  
                              M=4;  %固定价格 
    case 2  
                               M=5;  %优化价格
    case 3  
                               M=6;  %优化价格
    case 4  
                               M=7;  %优化价格
    otherwise  
                              disp('Invalid value of programmatic'); 
    end

%  EE=EEE;
    E=0.1*EEE;
e=E;
%-------------------仿真参数------------------%
h=27; %无人机高度
% M=10; %车辆数目
N=10; %无人机数目
L0=0.9; %慢衰落系数
MUk=0.5;
%g=rand(1); %信道随机变量
g=0.955477890177557;
w=10; %带宽GHz
Delta=1e-9;  %背景噪声
% E=0.8;%第一个阈值
e2=9;%SINR阈值设为9
% wwwww=(-e2)/((1/log(2))*log(1-E));
% Pr=1-exp(-(e2*I)/(P*G));

% W=[1900,140,400,200 500,430,610,90,330,405;1900,500,580,430,480,550,210,260,104,361];%车辆位置
% q=[2000,110,300,240,510,430,210,130,500,405;2000,420,580,430,480,550,210,460,344,361];%无人机初始位置

W=[1900,140,400,200 500,430,210,90,370,405;1900,500,580,430,480,550,210,260,104,361];%车辆位置
q=[2000,140,400,200,430,670,210,90,370,405;2000,500,580,430,480,720,210,260,104,361];%无人机初始位置

G=ones(M,N);
 for n=1:N     
  for m=1:M
d(m,n)=sqrt(h^2+norm([q(1,n) q(2,n)]-[W(1,m) W(2,m)]))
G(m,n)=g*L0*(d(m,n)^(-1.4));  %车车通信信道
  end 
 end
  for m=1:M
G(m,m)=25*G(m,m)      ;
  end
X=G+MUk;%xi部分  
  
 %%%%%%%%%%%%%%%%%%求解功率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pmax=0.5;%FUE最大功率值
T=20;
Lamda=zeros(T,M);
Mu=zeros(T,M);


global T
T =20;
SINR=zeros(T,M);
P=zeros(T,M);
 C=zeros(T,M);
P0=0.378*exp(e^2);
% P(1,:)=[0 7.5 7.5 7.5 7.5 7.5 7.5 7.5 7.5 7.5];
P(1,:)=7.5*(linspace(1,1,M)-[1, zeros(1, M-1)]);
%  C(1,:)=[0 400 400 400 400 400 400 400 400 400];
 C(1,:)=400*(linspace(1,1,M)-[1, zeros(1, M-1)]);
Lamda(1,:)=0.1;
Mu(1,:)=40;
rth=0.5;

 for t=1:T-1
     
    for m=2:M
        I(t,m)=P0*G(1,m)+P(t,:)*G(:,m)+Delta-P(t,m)*G(m,m);
        SINR(t,m)=P(t,m)*G(m,m)/I(t,m); 
    end
  
 
 for m=2:M  
      %A(t+1,m)= -rth*I(t,m)/(log(1-e))-P(t,m)*G(m,m)
      Lamda(t+1,m)=Lamda(t,m)+0.001/sqrt(t)*(-rth*I(t,m)/(log(1-e))-P(t,m)*G(m,m))
     b(t,m)=I(t,m)/G(m,m)     ;
     a(t,m)=w/(C(t,m)*G(m,1)-Lamda(t,m)*G(m,m));
  
 end   

       P(t+1,:)=(a(t,:)-b(t,:))*exp(e^2); %

%-------------------上层博弈------------------%
  for m=2:M
 z(t+1,m)=w*G(m,m)/(I(t,m)*a(t,m));
 Mu(t+1,m)=Mu(t,m)+0.001/sqrt(t)*(C(t,m)*G(m,1)*(a(t,m)-b(t,m))-Mu(t,m)*w*log(w*G(m,m)/(I(t,m)*a(t,m))))  ;   

c(t,m)=sqrt(w^2*Mu(t,m)^2-4*w*Lamda(t,m)*(Mu(t,m)-1)^2*I(t,m))
%c(t,m)=w^2*Mu(t,m)^2-4*w*Lamda(t,m)*(Mu(t,m)-1)^2*I(t,m)
%d(t,m)=(Mu(t,m)-1)*I(t,m)
d(t,m)=(Mu(t,m)-1)
g(t,m)=c(t,m)/d(t,m)+Lamda(t,m);
h(t,m)=Mu(t,m)/(2*Mu(t,m)-1)*I(t,m);
f(t,m)=g(t,m)-h(t,m);


C(t+1,m)=G(m,m)/G(m,1)*f(t,m);
%C(t+1,m)=G(m,m)/G(m,1)*(c(t,m)/((Mu(t,m)-1)*I(t,m))+(Lamda(t,m)-w*Mu(t,m)/(2*Mu(t,m)-1)*I(t,m)));
  end 

 C(t+1,1)=0; 
 end
 for m=1:M
     if m==1
%        Pr(m)=(P0)*(1-exp(-(e2*0.2)/(P0*G(m,m)))); 
        Pr(m)=0.1*(P(T-1,2))*(1-exp(-(e2*I(T-1,2))/(P(T-1,2)*G(2,2))))
     else
 Pr(m)=0.1*(P(T-1,m))*(1-exp(-(e2*I(T-1,m))/(P(T-1,m)*G(m,m)))); 
     end
 end


  speed=linspace(0.1,0.1*EEE,EEE);
 sumsinr=zeros(1,M);
      for m=1:1
        I00(T-1,m)=sum(P(T-1,:)*G(:,m)+Delta-P(T-1,m)*G(m,m),2); 
     end
    I0=I00(T-1);

 for t=1:T-1
     for m=1:M
     if m==1
       sumsinr(1,m)=(1/log(2))*log(1+(P0*G(1,m))/(0.01/w*I0+Delta));     
    else
     sumsinr(1,m)=w*(1/log(2))*log(SINR(T-1,m))-log(2*w)/L0; 
     end
     end
 end   
 
 for i=1:M
    if i==1
        Gdc = G(i,i);
    else
    Gdc=G(i,i)/I(T-1,i);
    end
    f = @(x) (1/log(2))*log(1+x*Gdc/Delta);
    g = @(x) 2*w+x*Gdc*I(T-1,i);
dc_sum(i)=dc_programming(g,f);
end

binarySearchForMin_sum=zeros(1,M);
for i=1:M
    if i==1
        GbinarySearchForMin = G(i,i);
        fx = @(x) -((1/log(2))*log(1+x*GbinarySearchForMin/Delta)-2*w+x*GbinarySearchForMin)-log(w);
        binarySearchForMin_sum(i)=binarySearchForMin(fx);
    else
    GbinarySearchForMin=G(i,i)/I(T-1,i);
    fx = @(x) -((1/log(2))*log(1+x*GbinarySearchForMin/Delta)-2*w+x*GbinarySearchForMin*I(T-1,i));
binarySearchForMin_sum(i)=binarySearchForMin(fx);
    end
end
%   Prr(EE,:)=(Pr);
%  Pee(EE,:)=(P(T,:));
%  sumsinrEE(EE)= sum(sumsinr,1)
sumsinrEE=zeros(3,programmatic);

sumsinrEE(1,EEm)= sum(sumsinr,2)

end 

%  data= sumsinrEE
% xxxx = 1:M;
% % bar(xxxx,data','stacked')
% bar(xxxx,data')
%  dc= [true true true true true true true true true true];
% %   dc= [false true false false false false false false false false];
%  Prrr=Prr(:,dc)
%  
 
 
 
%{
 
figure 
% plot(speed,Prrr(:,1),'-*r','linewidth',4);   %宏用户   D.C.
hold on
grid on
box on
% plot(speed,Prrr(:,2),'-*b','linewidth',1.5);   %宏用户   D.C.
% plot(speed,Prrr(:,3),'-*b','linewidth',1.5);   %宏用户   D.C.
% plot(speed,Prrr(:,4),'-*b','linewidth',1.5);   %宏用户   D.C.
% plot(speed,Prrr(:,5),'-*b','linewidth',1.5);   %宏用户   D.C.
% plot(speed,Prrr(:,6),'-*y','linewidth',1.5);   %宏用户   D.C.
% plot(speed,Prrr(:,7),'-*g','linewidth',1.5);   %宏用户   D.C.
% plot(speed,Prrr(:,8),'-*k','linewidth',1.5);   %宏用户   D.C.
% plot(speed,Prrr(:,9),'-*b','linewidth',1.5);   %宏用户   D.C.
% plot(speed,Prrr(:,10),'-*b','linewidth',1.5);   %宏用户   D.C.

plot(speed,Prrr(:,2),'-*m','linewidth',1);   %宏用户   D.C.
plot(speed,Prrr(:,3),'-*c','linewidth',1);   %宏用户   D.C.
% plot(speed,Prrr(:,4),'-*b','linewidth',1);   %宏用户   D.C.
plot(speed,Prrr(:,4),'-*', 'Color',[0.5, 0, 0.5]);
plot(speed,Prrr(:,5),'-*r','linewidth',1);   %宏用户   D.C.
plot(speed,Prrr(:,6),'-*y','linewidth',1);   %宏用户   D.C.
plot(speed,Prrr(:,7),'-*g','linewidth',1);   %宏用户   D.C.
% plot(speed,Prrr(:,8),'-*b','linewidth',1);   %宏用户   D.C.
plot(speed,Prrr(:,8),'-*', 'Color',[0, 0.5, 0]);   %宏用户   D.C.
% plot(speed,Prrr(:,9),'-*b','linewidth',1);   %宏用户   D.C.
plot(speed,Prrr(:,9),'-*', 'Color',[1, 0.5, 0]);   %宏用户   D.C.
plot(speed,Prrr(:,10),'-*k','linewidth',1);   %宏用户   D.C.
plot(speed,Prrr(:,1),'-*b','linewidth',1.5);   %宏用户   D.C.


% plot(speed,Prrr(:,2),'-or','linewidth',1);   %宏用户   D.C.        1
% plot(speed,Prrr(:,3),'-*b','linewidth',1);   %宏用户   D.C.        2
% % plot(speed,Prrr(:,4),'-*b','linewidth',1);   %宏用户   D.C.     
% plot(speed,Prrr(:,4),'-k','linewidth',1);     %    3
% plot(speed,Prrr(:,5),'-d', 'Color',[1, 0.5, 0],'linewidth',1);   %宏用户   D.C. 橘色是这个4
% plot(speed,Prrr(:,6),'-og','linewidth',1);   %宏用户   D.C.  5   [0.5, 0, 0.5] 表示线条的颜色为紫色。  'Color',[0.5, 0, 0.5]','linewidth',1
% plot(speed,Prrr(:,7),'-*','Color',[0.5, 0, 0.5]','linewidth',1);   %宏用户   D.C. 6
% % plot(speed,Prrr(:,8),'-*b','linewidth',1);   %宏用户   D.C. 
% plot(speed,Prrr(:,8),'-+y','linewidth',1);   %宏用户   D.C. 7
% % plot(speed,Prrr(:,9),'-*b','linewidth',1);   %宏用户   D.C. 
% plot(speed,Prrr(:,9),'-d', 'Color',[0, 0.5, 0],'linewidth',1);   %宏用户   D.C.  8  '-+b','linewidth',1
% plot(speed,Prrr(:,10),'-sc','linewidth',1);   %宏用户   D.C.    9
% plot(speed,Prrr(:,1),'-*b','linewidth',1.5);   %宏用户   D.C.   cue
set(gca,'xlim',[0,0.1*EEE]);
set(gca,'ylim',[0,0.1*EEE]);
legend('V1','V2','V3','V4','V5','V6','V7','V8','V9','CUE')
% plot([0,0.1*EEE],[0,10*EEE],'k','linewidth',3);
xlabel('中断概率阈值', 'FontSize', 14);
ylabel('真实的中断概率', 'FontName', 'SimSun', 'FontSize', 14);
kk=1;
x0=0.15;y0=0.151;
b=y0-kk*x0;
x=-5:20;
y=kk*x+b;
plot(x,y,'k','linewidth',3)


% plot(C(:,2),'-or','linewidth',1.5);
% hold on
% plot(C(:,3),'-*b','linewidth',1.5);
% plot(C(:,4),'-sk','linewidth',1.5);
% plot(C(:,5),'-d','Color',[1, 0.5, 0],'linewidth',1.5);
% plot(C(:,6),'-og','linewidth',1.5);
% plot(C(:,7),'-*', 'Color',[0.5, 0, 0.5]','linewidth',1.5);
% plot(C(:,8),'-+b','linewidth',1.5);
% plot(C(:,9),'-d', 'Color',[0, 0.5, 0],'linewidth',1.5);
% plot(C(:,10),'-sc','linewidth',1.5);

%}



% figure 
% TT=linspace(1,10,10);
% x_total=[sumsinr',binarySearchForMin_sum',dc_sum'];
% 
% data = x_total';
% b = bar(TT,data');
% ylim([5 18]);
% xlim([0 11])  
% grid on
% 
% % set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4),aqa(5)})
% set(b(1),'FaceColor',[240,1,233]/255)     
% set(b(2),'FaceColor',[0,250,8]/255)    
% set(b(3),'FaceColor',[255,255,0]/255)    
% % set(b(4),'FaceColor',[128,44,140]/255) 
% 
% ch = get(b,'children');
% set(gca,'XTickLabel',{'CUE','V1','V2','V3','V4','V5','V6','V7','V8','V9'})
% legend('鲁棒博弈论','二分法','D.C.规划')
% xlabel('通信用户','FontSize',12);
% ylabel('传输速率 (Mb/s)','FontSize',12);


%{
for i=1:M
if i==1
        Gdc = G(i,i);
    else
    Glyapunov=G(i,i)/I(T-1,i);
    lyapunov = @(x) -((1/log(2))*log(1+x*Glyapunov/Delta)-2*w+x*Glyapunov*I(T-1,i));
    
%     glyapunov = @(x) 2*w+x*Gdc*I(T-1,i);
%     flyapunov =lyapunov-glyapunov;
lyapunov_sum(i)=lyapunovOptimization(lyapunov);
end
end
%}
