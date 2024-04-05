clc
close all
clear
%-------------------系统参数-------------------%研究了在某个时间段5个车簇在通信，每个车簇有一个发射CM车与簇头CH接收车进行通信
%-------------------系统参数-------------------%发送车CM速度32，30，34，34，30.相应的CH接收车速度30，34，35，36，33
  L=zeros(5,5);
for F=1:6
V1=28*ones(1,5)+2*F*ones(1,5);%发送车CM速度
% V1=[0,0,0,0,0];
V2=[30,30,30,30,30];%接收车CH速度
v=zeros(5,5);%车与车之间的相对速度矩阵
h=zeros(5,5);%车与车之间的先前时刻的值
g=zeros(5,5);%误差，服从参数为1指数分布
G=zeros(5,5);%车与车之间信道增益矩阵
I=zeros(5,5);%车与车之间干扰信道增益矩阵（对角线全为0元素，因为对角线是有效信道增益）
for i=1:5
    for j=1:5
        v(i,j)=abs(V1(i)-V2(j));
        h(i,j)=1+0.3*rand(1);
    end
end
% [6,1,4,0,3;0,7,2,6,9;8,1,6,2,1;10,3,8,4,1;0,7,2,6,9];相对速度矩阵
%h=[0.845418353066450,0.995224514302217,0.687960311535798,0.741632431164121,0.706496751675251;0.716599213586070,0.675355816604331,0.609143849450247,0.779761672779394,0.697456077357881;0.947490617895953,0.811444306977048,0.965654067256766,0.989572070823563,0.834170383464932;0.647590153069472,0.970613035572538,0.837424308992463,0.953446121969907,0.769790372435062;0.842902954353616,0.628305432990921,0.969908972473279,0.856831766456680,0.641799867159030];
h=[1.21026765186132,1.00186773622345,1.11230368038275,1.27044868578780,1.09550344842318;1.17912487696444,1.08933855474170,1.03750431438572,1.11650670059712,1.24530632319216;1.29435271127300,1.25859695058615,1.02514626462880,1.10131359710114,1.07083874991164;1.09534161637965,1.29533452164231,1.16447526737693,1.22477537780633,1.25255552041723;1.05006691209611,1.27092928278688,1.03153723834667,1.22352792893222,1.21881158135137];
Delta=1e-6;%背景噪声
T=0.001;%所有车与车移动信道CSI采样周期（V2V）
c=3*1e+8;%光速
fc=5.9*1e+9;%多普勒频移中心频率
j1=2*pi*fc*T/c*v;%0阶贝塞尔函数参数（V2V）
epsi=besselj(0,j1);%求贝塞尔函数值（V2V）
a=(epsi.^2).*h;%车车通信移动信道快衰落的模型中的CSI反馈部分（先前时刻的信道状态信息）
l=0.0001;%车车通信慢衰落系数
G1=l*a;%车车通信衰落（先前时刻的信道增益（也是测得的平均信道增益））
G2=l*(1-epsi.^2);%车车通信信道衰落误差部分的平均信道增益
G=G1+G2;%车车通信总的平均信道增益矩阵（该矩阵用于求取平均信噪比，长期速率和）
  for i=1:5
      G(i,i)=100*G(i,i);%有效信道增益的倍数扩大（簇内通信）
  end 
MUk=0.5;
alfak=l*(1-epsi.^2);
betak=l*(1-epsi.^2);
X=G1+MUk*alfak+betak;%%xi部分
SIGMA=0.5;
E=0.1;
H=X+SIGMA*sqrt(-2*log(E))*alfak; 
% ----------构造贝恩斯坦近似干扰约束的系数矩阵----------% 
I=H;%贝恩斯坦近似约束矩阵
    for i=1:5
      I(i,i)=0;%构造干扰约束中的系数矩阵
    end
    
global N
N =25;
arfa=zeros(N,5);
beta=zeros(N,5);
SINR=zeros(N,5);
Pmax=log(0.5)*ones(N,5);
B=1;
C1=linspace(1,1,5);
GF=linspace(1.5,1.5,5);%功放系数
PC=0.05;%电路中固定消耗的功率
S=zeros(N,5);
s=zeros(N,5);
Gama=0.1;%原来是
Ith=2e-6;%干扰阈值Gama要给相应的存储空间
Temp=zeros(40,5);%每次迭代收敛的最优功率暂存区 40应为Z
Z=1;

while 1
P=zeros(N,5);
P(1,:)=[-8,-8,-8,-8,-8]; 
Lamda=zeros(N,5);
Lamda(1,:)=linspace(20,20,5);
for k=1:N
    for i=1:5
    I1(i)=exp(P(k,:))*G(:,i)+Delta-G(i,i)*exp(P(k,i));
    SINR(k,i)=G(i,i)*exp(P(k,i))/I1(i);     
    arfa(k,i)=SINR(k,i)/(1+SINR(k,i));
    Total_T(k)=B*dot(log(C1+SINR(k,:))/log(2),C1);
    end
    for i=1:5
        for j=1:5
            if j~=i
               S(k,j)=B/log(2)*arfa(k,j)*G(i,j)*SINR(k,j)/((exp(P(k,j))*G(j,j)))+Lamda(k,j)*I(i,j);%%到底是G(i,j)还是G(j,i),应该是G(i,j)。当j不变是第j列，代表的是各个CM向某个簇头传输的信息。
            else
               S(k,j)=0;
            end
         end
          s(k,i)=sum(S(k,:))+Gama(Z)*GF(i);
    end
%  Total_T(k+1)=B*dot(log(C1+SINR(k,:))/log(2),C1); 
P(k+1,:)=min(log(B.*arfa(k,:)/log(2))-log(s(k,:)),Pmax(k,:));
     if max(abs(exp(P(k+1,:))-exp(P(k,:))))<=1e-5
        P(k+1,:)=P(k,:);
     end
     for i=1:5
         Lamda(k+1,i)=max(0,Lamda(k,i)+1000000/k*(exp(P(k,:))*I(:,i)-Ith));
         w(k,i)=exp(P(k,:))*I(:,i)-Ith;
     end 
     if  max(abs(w(k,:)))<=3e-7    %该数值可改，表示乘子在一定的程度下就认为已经稳定了
         Lamda(k+1,:)=Lamda(k,:);
     end
end
D(Z)=Total_T(N);
Temp(Z,:)=P(N,:);
error =abs(D(Z)-Gama(Z)*(GF*exp(Temp(Z,:)')+PC));
Gama(Z+1)=D(Z)/(GF*exp(Temp(Z,:)')+PC);
      if  error<=1e-3
          L(F,:)=exp(Temp(Z,:));
          Q(F)=D(Z);
          A(F)=Gama(Z+1);
break;
%       else
%        Lamda=zeros(N,5);
%     Lamda(1,:)=linspace(4,6,5);
%         P=zeros(N,5);
%         P =Temp(Z,:); 
      end 
          Z=Z+1;
end
end
   
%{
plot(10*log(exp(P(:,1)))/log(10),'-+r');
hold on
plot(10*log(exp(P(:,2)))/log(10),'-og');
plot(10*log(exp(P(:,3)))/log(10),'-*b');
plot(10*log(exp(P(:,4)))/log(10),'-sk');
plot(10*log(exp(P(:,5)))/log(10),'-dc');
legend('P1','P2','P3','P4','P5')
xlabel('iteration');
ylabel('Power conversion value(dB)');
%}
  Q1=Q;
  A1=A;
plot(exp(P(:,1)),'-+r');
hold on
plot(exp(P(:,2)),'-og');
plot(exp(P(:,3)),'-*b');
plot(exp(P(:,4)),'-sk');
plot(exp(P(:,5)),'-dc');

legend('P1','P2','P3','P4','P5')
xlabel('iteration');
ylabel('Power conversion value(W)');

figure
plot(Lamda(:,1),'-+r');
hold on
plot(Lamda(:,2),'-og');
plot(Lamda(:,3),'-*b');
plot(Lamda(:,4),'-sk');
plot(Lamda(:,5),'-dc');
legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5')
xlabel('iteration');
ylabel('Multiplier convergence');

figure
plot(Gama,'-sb','linewidth',1.5) 
hold on
plot(D,'-^r','linewidth',1.5)

hold on

[AX,H1,H2] =plotyy(Z,Gama,Z,D,@plot);% 获取坐标轴、图像句柄
set(AX(1),'XColor','k','YColor','b','LineWidth',1.5);
set(AX(2),'XColor','k','YColor','r','LineWidth',1.5)

set(get(AX(1),'ylabel'),'color','blue','string', 'Average EE(bits/Joule/Hertz)','fontsize',12);
set(get(AX(2),'ylabel'),'color','red','string', 'Average SE(bits/sec/Hertz)','fontsize',12);

xlabel('iteration','fontsize',12);
legend('EE','SE')


    FF=[0,2,4,6,8,10];
    figure
 [AX,H1,H2] = plotyy(FF,A1,FF,Q1,'plot');

set(AX(1),'XColor','k','YColor','b','LineWidth',1.5);
set(AX(2),'XColor','k','YColor','r','LineWidth',1.5);

HH1=get(AX(1),'Ylabel');
set(HH1,'String','Average EE(bits/Hz/Joule)','fontsize',12);
set(HH1,'color','b');

HH2=get(AX(2),'Ylabel');
set(HH2,'String','Average SE(bits/s/Hz)','fontsize',12);
set(HH2,'color','r');

set(H1,'marker','s')
set(H1,'LineStyle','-','LineWidth',1.5);
set(H1,'color','b');

set(H2,'LineStyle','-','LineWidth',1.5);
set(H2,'color','r');
set(H2,'marker','^')
% hold on 
% EE=[1.53634734393827,1.53579347534572,1.53541361465141,1.53510426222248,1.53482959098909];
% SE=[0.231773787999855,0.231780191199803,0.231784581226184,0.231788155496835,0.231791328385694];

legend([H1,H2],{'EE';'SE'});
xlabel('\Delta\upsilon (m/s)','fontsize',12);
legend('EE,C=C^*','SE,C=C^*')


% set(H1,'Linestyle','--');
% 
% set(H2,'Linestyle',':');

% set(gcf,'color','white')
% set(gca,'xcolor','r');
% set(gca,'ycolor','r','linewidth',1.5) %设置边框宽度  
