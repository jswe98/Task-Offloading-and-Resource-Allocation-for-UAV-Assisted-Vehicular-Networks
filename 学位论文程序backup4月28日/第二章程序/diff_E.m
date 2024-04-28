clc
close all
clear
digits(3);
EEE=6;
global T
for EE=1:EEE
E=0.1*EE;
e=E;
%-------------------仿真参数------------------%
h=27; %无人机高度
M=10; %车辆数目
N=10; %无人机数目
L0=0.9; %慢衰落系数
MUk=0.5; %g=rand(1); %信道随机变量
g=0.955477890177557;
w=10; %带宽GHz
Delta=1e-9;  %背景噪声
e2=9;%SINR阈值设为9
W=[1900,140,400,200 500,430,610,90,330,405;1900,500,580,430,480,550,210,260,104,361];%车辆位置
q=[2000,110,300,240,510,430,210,130,500,405;2000,420,580,430,480,550,210,460,344,361];%无人机初始位置
G=ones(M,N);
d=ones(M,N);
    for n=1:N     
        for m=1:M
            d(m,n)=sqrt(h^2+norm([q(1,n) q(2,n)]-[W(1,m) W(2,m)]));
            G(m,n)=g*L0*(d(m,n)^(-1.4));  %车车通信信道
        end 
    end
    for m=1:M
        G(m,m)=25*G(m,m)      ;
    end
X=G+MUk;%xi部分  
%%%%%%%%%%%%%%%%%%求解功率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pmax=0.5;%FUE最大功率值
T =20;
Lamda=zeros(T,10);
Mu=zeros(T,10);
SINR=zeros(T,10);
P=zeros(T,10);
I=zeros(T,10);
b=zeros(T,10);
a=b;f=a;c=f;z=c;g=z;h=g;
P0=0.378*exp(e^2);
P(1,:)=[0 7.5 7.5 7.5 7.5 7.5 7.5 7.5 7.5 7.5];
C(1,:)=[0 400 400 400 400 400 400 400 400 400];
Lamda(1,:)=0.1;
Mu(1,:)=40;
rth=0.5;
    for t=1:T-1    
        for m=2:10
            I(t,m)=P0*G(1,m)+P(t,:)*G(m,:)'+Delta-P(t,m)*G(m,m);
            SINR(t,m)=P(t,m)*G(m,m)/I(t,m); 
        end
        for m=2:10    
            Lamda(t+1,m)=Lamda(t,m)+0.001/sqrt(t)*(-rth*I(t,m)/(log(1-e))-P(t,m)*G(m,m))
            b(t,m)=I(t,m)/G(m,m)     ;
            a(t,m)=w/(C(t,m)*G(m,1)-Lamda(t,m)*G(m,m)); 
        end   
            P(t+1,:)=(a(t,:)-b(t,:))*exp(e^2)-4*exp(e^2); %
%-------------------上层博弈------------------%
        for m=2:10
            z(t+1,m)=w*G(m,m)/(I(t,m)*a(t,m));
            Mu(t+1,m)=Mu(t,m)+0.001/sqrt(t)*(C(t,m)*G(m,1)*(a(t,m)-b(t,m))-Mu(t,m)*w*log(w*G(m,m)/(I(t,m)*a(t,m))))  ;   
            c(t,m)=sqrt(w^2*Mu(t,m)^2-4*w*Lamda(t,m)*(Mu(t,m)-1)^2*I(t,m))
            d(t,m)=(Mu(t,m)-1)
            g(t,m)=c(t,m)/d(t,m)+Lamda(t,m);
            h(t,m)=Mu(t,m)/(2*Mu(t,m)-1)*I(t,m);
            f(t,m)=g(t,m)-h(t,m)
            C(t+1,m)=G(m,m)/G(m,1)*f(t,m);
        end 
            C(t+1,1)=0; 
     end
     for m=1:10
         if m==1
            Pr(m)=0.1*(P(T-1,2))*(1-exp(-(e2*I(T-1,2))/(P(T-1,2)*G(2,2))))
         else
            Pr(m)=0.1*(P(T-1,m))*(1-exp(-(e2*I(T-1,m))/(P(T-1,m)*G(m,m)))); 
         end
     end
%  Prr=zeros(EEE,10);
%  Pee=zeros(EEE,10);
 Prr(EE,:)=(Pr);
 Pee(EE,:)=(P(T,:))
 end
speed=linspace(0.1,0.1*EEE,EEE);
sumsinr=zeros(1,M);
dc= [true true true true true true true true true true];
Prrr=Prr(:,dc)
figure 
hold on
grid on
box on
plot(speed,Prrr(:,2),'-*m','linewidth',1);       %车用户  
plot(speed,Prrr(:,3),'-*c','linewidth',1);       %车用户   
plot(speed,Prrr(:,4),'-*', 'Color',[0.5,0,0.5]); %车用户
plot(speed,Prrr(:,5),'-*r','linewidth',1);       %车用户  
plot(speed,Prrr(:,6),'-*y','linewidth',1);       %车用户  
plot(speed,Prrr(:,7),'-*g','linewidth',1);       %车用户  
plot(speed,Prrr(:,8),'-*', 'Color',[0,0.5,0]);   %车用户  .
plot(speed,Prrr(:,9),'-*', 'Color',[1,0.5,0]);   %车用户  
plot(speed,Prrr(:,10),'-*k','linewidth',1);      %车用户  
plot(speed,Prrr(:,1),'-*b','linewidth',1.5);     %宏用户 
set(gca,'xlim',[0,0.1*EEE]);
set(gca,'ylim',[0,0.1*EEE]);
legend('V1','V2','V3','V4','V5','V6','V7','V8','V9','CUE')
xlabel('\fontname{宋体}中断概率阈值', 'FontSize', 14);
ylabel('\fontname{宋体}真实的中断概率', 'FontName', 'SimSun', 'FontSize', 14);
kk=1;
x0=0.15;y0=0.151;
b=y0-kk*x0;
x=-5:20;
y=kk*x+b;
plot(x,y,'k','linewidth',3)
set(gca,'FontName','Times New Roman')

%{
     for m=1:1
        I00(T-1,m)=sum(P(T-1,:)*G(m,:)'+Delta-P(T-1,m)*G(m,m),2); 
     end
    I0=I00(T-1);
 for t=1:T-1
     for m=1:10
     if m==1
       sumsinr(1,m)=w*(1/log(2))*log(1+(P0*G(1,m))/(0.01/w*I0+Delta));     
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
