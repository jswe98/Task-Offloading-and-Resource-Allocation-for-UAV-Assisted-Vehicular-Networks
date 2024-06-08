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
 PPP=P*1e-3;
 sumsinr=zeros(1,M);
 for t=1:T-1
     for m=2:10
     sumsinr(1,m)=w*(1/log(2))*log(SINR(T-1,m)); 
     end
 end    
figure 

plot(PPP(:,2),'-om','linewidth',1.5);
hold on
plot(PPP(:,3),'-*c','linewidth',1.5);
plot(PPP(:,4),'-s', 'Color',[0.5, 0, 0.5],'linewidth',1.5);  %,'-*', 'Color',[0.5, 0, 0.5]
plot(PPP(:,5),'-r','linewidth',1.5);
plot(PPP(:,6),'-oy','linewidth',1.5);
plot(PPP(:,7),'-*g','linewidth',1.5);
plot(PPP(:,8),'-+', 'Color',[0, 0.5, 0]','linewidth',1.5);
plot(PPP(:,9),'-d', 'Color',[1, 0.5, 0],'linewidth',1.5);
plot(PPP(:,10),'-*k','linewidth',1);

% plot(P(:,2)*1e-3,'-or','linewidth',1.5);
% hold on
% plot(P(:,3)*1e-3,'-*b','linewidth',1.5);
% plot(P(:,4)*1e-3,'-sk','linewidth',1.5);
% plot(P(:,5)*1e-3,'-dc','linewidth',1.5);
% plot(P(:,5)*1e-3,'-d','Color',[1, 0.5, 0],'linewidth',1.5);
% plot(P(:,6)*1e-3,'-og','linewidth',1.5);
% % plot(P(:,7)*1e-3,'-*r','linewidth',1.5);
% plot(P(:,7)*1e-3,'-*', 'Color',[0.5, 0, 0.5]','linewidth',1.5);
% plot(P(:,8)*1e-3,'-+b','linewidth',1.5);
% % plot(P(:,9)*1e-3,'-dk','linewidth',1.5);
% plot(P(:,9)*1e-3,'-d', 'Color',[0, 0.5, 0],'linewidth',1.5);
% plot(P(:,10)*1e-3,'-sc','linewidth',1.5);
% %plot(P(:,11),'-dc');

legend('P1','P2','P3','P4','P5','P6','P7','P8','P9')
xlabel('迭代次数', 'FontName', 'SimSun','FontSize',14);
ylabel('功率(瓦)', 'FontName', 'SimSun','FontSize',14);  
xlim([0 20])


figure 
plot(C(:,2),'-om','linewidth',1.5);
hold on
plot(C(:,3),'-*c','linewidth',1.5);
plot(C(:,4),'-s', 'Color',[0.5, 0, 0.5],'linewidth',1.5);  %,'-*', 'Color',[0.5, 0, 0.5]
plot(C(:,5),'-r','linewidth',1.5);
plot(C(:,6),'-oy','linewidth',1.5);
plot(C(:,7),'-*g','linewidth',1.5);
plot(C(:,8),'-+', 'Color',[0, 0.5, 0]','linewidth',1.5);
plot(C(:,9),'-d', 'Color',[1, 0.5, 0],'linewidth',1.5);
plot(C(:,10),'-*k','linewidth',1);
%plot(P(:,11),'-dc');


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

legend('C1','C2','C3','C4','C5','C6','C7','C8','C9')
xlabel('\fontname{宋体}迭代次数', 'FontName', 'SimSun','FontSize',14);
ylabel('\fontname{宋体}价格', 'FontName', 'SimSun','FontSize',14);   
xlim([0 20])
 set(gca,'FontName','Times New Roman')
% plot(speed,Prrr(:,2),'-*m','linewidth',1);   %宏用户   D.C.
% plot(speed,Prrr(:,3),'-*c','linewidth',1);   %宏用户   D.C.
% plot(speed,Prrr(:,4),'-*', 'Color',[0.5, 0, 0.5]);
% plot(speed,Prrr(:,5),'-*r','linewidth',1);   %宏用户   D.C.
% plot(speed,Prrr(:,6),'-*y','linewidth',1);   %宏用户   D.C.
% plot(speed,Prrr(:,7),'-*g','linewidth',1);   %宏用户   D.C.
% plot(speed,Prrr(:,8),'-*', 'Color',[0, 0.5, 0]);   %宏用户   D.C.
% plot(speed,Prrr(:,9),'-*', 'Color',[1, 0.5, 0]);   %宏用户   D.C.
% plot(speed,Prrr(:,10),'-*k','linewidth',1);   %宏用户   D.C.
% plot(speed,Prrr(:,1),'-*b','linewidth',1.5);   %宏用户   D.C.

