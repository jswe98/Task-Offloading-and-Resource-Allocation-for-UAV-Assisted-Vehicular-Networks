clc
close all
clear

%-------------------仿真参数------------------%
h=30; %无人机高度
M=10; %车辆数目
N=10; %无人机数目
L0=0.9; %慢衰落系数
%g=rand(1); %信道随机变量
g=0.955477890177557;
w=10; %带宽GHz
Delta=1e-9;  %背景噪声


W=[1900,140,400,200 500,430,210,90,370,405;1900,500,580,430,480,550,210,260,104,361];%车辆位置
q=[2000,140,400,200,430,670,210,90,370,405;2000,500,580,430,480,720,210,260,104,361];%无人机初始位置

G=ones(M,N);
 for n=1:N     
  for m=1:M
d(m,n)=sqrt(h^2+norm([q(1,n) q(2,n)]-[W(1,m) W(2,m)]));
G(m,n)=g*L0*(d(m,n)^(-1.4));  %车车通信信道
  end 
 end
  for m=1:M
G(m,m)=25*G(m,m)      
  end
 %%%%%%%%%%%%%%%%%%求解功率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pmax=0.5;%FUE最大功率值
T=20;
Lamda=zeros(T,10);
Mu=zeros(T,10);

e=0.2;
global T
T =20;
SINR=zeros(T,10);
P=zeros(T,10);
P0=0.001
P(1,:)=[0 7.5 7.5 7.5 7.5 7.5 7.5 7.5 7.5 7.5];
 C(1,:)=[0 400 400 400 400 400 400 400 400 400];
Lamda(1,:)=0.1;
Mu(1,:)=40;
rth=0.5

 for t=1:T-1
     
    for m=2:10
        I(t,m)=P0*G(1,m)+P(t,:)*G(m,:)'+Delta-P(t,m)*G(m,m);
        SINR(t,m)=P(t,m)*G(m,m)/I(t,m); 
    end
  
 
 for m=2:10    
      %A(t+1,m)= -rth*I(t,m)/(log(1-e))-P(t,m)*G(m,m)
      Lamda(t+1,m)=Lamda(t,m)+0.001/sqrt(t)*(-rth*I(t,m)/(log(1-e))-P(t,m)*G(m,m));
     b(t,m)=I(t,m)/G(m,m)     
     a(t,m)=w/(C(t,m)*G(m,1)-Lamda(t,m)*G(m,m))
  
 end   

     P(t+1,:)=(a(t,:)-b(t,:));
    
%-------------------上层博弈------------------%
  for m=2:10
 z(t+1,m)=w*G(m,m)/(I(t,m)*a(t,m))
 Mu(t+1,m)=Mu(t,m)+0.001/sqrt(t)*(C(t,m)*G(m,1)*(a(t,m)-b(t,m))-Mu(t,m)*w*log(w*G(m,m)/(I(t,m)*a(t,m))));     

c(t,m)=sqrt(w^2*Mu(t,m)^2-4*w*Lamda(t,m)*(Mu(t,m)-1)^2*I(t,m))
%c(t,m)=w^2*Mu(t,m)^2-4*w*Lamda(t,m)*(Mu(t,m)-1)^2*I(t,m)
%d(t,m)=(Mu(t,m)-1)*I(t,m)
d(t,m)=(Mu(t,m)-1);
g(t,m)=c(t,m)/d(t,m)+Lamda(t,m)
h(t,m)=Mu(t,m)/(2*Mu(t,m)-1)*I(t,m)
f(t,m)=g(t,m)-h(t,m)


C(t+1,m)=G(m,m)/G(m,1)*f(t,m)
%C(t+1,m)=G(m,m)/G(m,1)*(c(t,m)/((Mu(t,m)-1)*I(t,m))+(Lamda(t,m)-w*Mu(t,m)/(2*Mu(t,m)-1)*I(t,m)));
  end 
 C(t+1,1)=0; 
 end
 sumsinr=zeros(1,M);
 
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

figure 
TT=linspace(1,10,10);
x_total=[sumsinr',binarySearchForMin_sum',dc_sum'];

data = x_total';
% xxxx = 1:10;
% bar(xxxx,data','stacked')

b = bar(TT,data');
ylim([5 18]);
xlim([0 11])  
grid on

% set(gca,'XTickLabel',{aqa(1),aqa(2),aqa(3),aqa(4),aqa(5)})
% set(b(1),'FaceColor',[240,1,233]/255)     
% set(b(2),'FaceColor',[0,250,8]/255)    
% set(b(3),'FaceColor',[255,255,0]/255)    

set(b(1),'FaceColor',[126,153,244]/255)     
set(b(2),'FaceColor',[204,124,113]/255)    
set(b(3),'FaceColor',[122,182,86]/255)  
% set(b(4),'FaceColor',[128,44,140]/255) 

ch = get(b,'children');
set(gca,'XTickLabel',{'CUE','V1','V2','V3','V4','V5','V6','V7','V8','V9'})
legend('\fontname{宋体}鲁棒博弈论','\fontname{宋体}二分法','\fontname{Times New Roman}D.C.\fontname{宋体}规划')
xlabel('\fontname{宋体}通信用户','FontSize',14);
ylabel('\fontname{宋体}传输速率 \fontname{Times New Roman}(Mb/s)', 'FontName', 'SimSun','FontSize',14);
 set(gca,'FontName','Times New Roman')

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
