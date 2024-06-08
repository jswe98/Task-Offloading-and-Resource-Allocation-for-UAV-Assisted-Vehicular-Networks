clc
close all
clear

%-------------------仿真参数------------------%
R=500;%基站半径
r=30;%D2D通信半径

z=0.01;%功率价格
p00=1e-4;%未被复用信道时的宏用户发射功率
Delta=1e-6;%噪声
Klambda=0.1;%迭代步长
Knu=0.1;
Keta=0.1;
Ith=1e-3;%干扰阈值
rth=1e-1;%信噪比阈值
T1=0.0001;%V2V移动信道CSI采样周期T1<T2
T2=0.0005;%V2I信道链路CSI采样周期
c=3*1e+8;%光速
fc=6*1e+9;%多普勒频移中心频率
tau=100;%数据包长度
varsigma=200;%数据包流量强度
Dmax=0.01;%时延阈值
Omax=0.9;%递包率阈值
e1=0.1;%中断概率
e2=0.1;
e3=0.1;

%-------------------信道增益------------------%发送车速度32，30，34，30.相应的接收车32，30，34，30
V=[32,30,34,30];%发送车速度
v1=zeros(4,4);%不同簇内车与车之间的相对速度矩阵
v2=zeros(1,4);%车与基站之间的相对速度矩阵
h1=zeros(4,4);%车与车之间的先前时刻的值
h2=zeros(1,4);%车与基站之间的先前时刻的值
g1=zeros(4,4);%误差，服从参数为1指数分布
g2=zeros(1,4);%误差，服从参数为1指数分布

for i=1:4
    v2(i)=V(i);
    h2(i)=rand(1);
    g2(i)=min(exprnd(1),1);
    for j=1:4
        v1(i,j)=abs(V(i)-V(j));
        h1(i,j)=0.6+0.4*rand(1);
        g1(i,j)=min(exprnd(1),1);
    end
end


v1=[0,2,2,2;2,0,4,0;2,4,0,4;2,0,4,0;];
g1=[1,0.630579842689046,0.0552694482734952,0.934748595037084;0.0289143824639772,0.797788371141981,0.375948324473220,0.430720222787889;1,0.0204423033898084,0.519237069647440,1;0.939898370678918,1,0.981569615142164,0.832315647916305;];
h1=[0.877761563873822,0.603903459954064,0.711756786101186,0.962577306604001;0.934868254150822,0.622773154174099,0.832988120694515,0.887773101038399;0.832632833285146,0.623061744485949,0.713929490744247,0.984864412404504;0.973159158327232,0.709286683199985,0.758843537097381,0.652445882817202;];
h2=[0.672270237457429,0.0248552338448721,0.726914550912652,0.193039815977027;];
g2=[0.840461852451740,0.398335390200032,0.983906875129379,1;];
v2=[32,30,34,30;];

d0=[50];%宏用户到基站的距离
H0=g2(1,1)*d0.^(-1.0);%g0,0
d1=[10 13 9 12];%簇成员到簇头的距离
H1=g2(1,1:4).*d1.^(-1.0);%gi,i


j1=2*pi*fc*T1/c*v1;%0阶贝塞尔函数参数（V2V）
j2=2*pi*fc*T2/c*v2;%0阶贝塞尔函数参数（V2I）
epsi1=besselj(0,j1);%求贝塞尔函数值（V2V）
epsi2=besselj(0,j2);%求贝塞尔函数值（V2I）
a1=(epsi1.^2).*(h1.^2)+(1-epsi1.^2).*g1;%（V2V）通信移动信道快衰落的模型
a2=(epsi2.^2).*(h2.^2)+(1-epsi2.^2).*g2;%（V2I）通信移动信道快衰落的模型
l1=0.003;%车车通信慢衰落系数
l2=0.01;%V2I慢衰落系数
G1=l1.*a1;%gi,j
G2=l2.*a2;%gi,0



  A = zeros(5,5);
  G= zeros(5,5);
 % vpa(A,5);%设置A中每个元素为d小数位精度 
  A(1,1)=H0;
  A(2,2)=H1(1);
  A(3,3)=H1(2);
  A(4,4)=H1(3);
  A(5,5)=H1(4);
  A(1,2:5)=G2;
  A(2,3:5)=G1(1,2:4);
  A(3,4:5)=G1(2,3:4);
  A(4,5)=G1(3,4);
  A(2,1)=0;
  A(3,1:2)=0;
  A(4,1:3)=0;
  A(5,1:4)=0;
  G=triu(A,0)+tril(A',-1);
  
   
  %-------------------下层博弈------------------%
  SIGMA=0.5;
  
  GG1=l1.*(epsi1.^2).*(h1.^2);%求先前信道增益^g i,j
  GG2=l2.*(epsi2.^2).*(h2.^2);
  B = zeros(5,5);
  GG= zeros(5,5);
 % vpa(A,5);%设置A中每个元素为d小数位精度 
  B(1,1)=H0;
  B(2,2)=H1(1);
  B(3,3)=H1(2);
  B(4,4)=H1(3);
  B(5,5)=H1(4);
  B(1,2:5)=GG2;
  B(2,3:5)=GG1(1,2:4);
  B(3,4:5)=GG1(2,3:4);
  B(4,5)=GG1(3,4);
  B(2,1)=0;
  B(3,1:2)=0;
  B(4,1:3)=0;
  B(5,1:4)=0;
  GG=triu(B,0)+tril(B',-1);
  
  GGG1=l1.*sqrt(1-epsi1.^2);%alfai,betai
  GGG2=l2.*sqrt(1-epsi2.^2);
  C = zeros(5,5);
  alfai=zeros(5,5);
  betai=zeros(5,5);
 % vpa(A,5);%设置A中每个元素为d小数位精度 
  C(1,1)=H0;
  C(2,2)=H1(1);
  C(3,3)=H1(2);
  C(4,4)=H1(3);
  C(5,5)=H1(4);
  C(1,2:5)=GGG2;
  C(2,3:5)=GGG1(1,2:4);
  C(3,4:5)=GGG1(2,3:4);
  C(4,5)=GGG1(3,4);
  C(2,1)=0;
  C(3,1:2)=0;
  C(4,1:3)=0;
  C(5,1:4)=0;
  alfai=triu(C,0)+tril(C',-1); 
  betai=alfai;

  
  
  MUi=0.5;
  W=0.001;
  Klamda=0.01
  X=GG+MUi*alfai+betai ;%%d2d-v用户到蜂窝内基站的先前的信道增益


E=0.01;
C=1/sqrt(-2*log(E));
global N
N =20;
arfa=zeros(N,5);
SINR=zeros(N,5);
P=zeros(N,5);
P(1,:)=[0.002 0.002 0.002 0.002 0.002];
Pmax=2;
Lamda2=zeros(N,5);
Lamda3=zeros(N,5);
Lamda4=zeros(N,5);
Lamda5=zeros(N,5);
NU=zeros(N,5);
Fai=zeros(N,5);
Lamda2(1,:)=1;
Lamda3(1,:)=1;
Lamda4(1,:)=1;
Lamda5(1,:)=1;

F=zeros(N,5);
F(1,:)=[0 70 70 70 70];

NU(1,:)=0.1;
Fai(1,:)=0.1;


MU2=zeros(N,1);
MU3=zeros(N,1);
MU4=zeros(N,1);
MU5=zeros(N,1);

 S=zeros(N,5);
 s=zeros(N,1);
 C1=linspace(1,1,5);
 for t=1:N
    for i=1:5
        I1(i)=G(i,:)*P(t,:)'+Delta-G(i,i)*P(t,i);
        SINR(t,i)=G(i,i)*P(t,i)/I1(i); 
    arfa(t,i)=SINR(t,i)/(1+SINR(t,i));
    end
  
 
 for i=1:5     
     MU2(t)=C*sum(Lamda2(t,:));
     MU3(t)=C*sum(Lamda3(t,:));
     MU4(t)=C*sum(Lamda4(t,:));
     MU5(t)=C*sum(Lamda5(t,:));
    % s(i,t+1)=max(s(i,t)+k1/sqrt(t)*(A*(G(i,:)*p-G(i,i)*p(i)+h)-G(i,i)*p(i)),0);
    % u(i,t+1)=max(u(i,t)+k2/sqrt(t)*(B*(G(i,:)*p-G(i,i)*p(i)+h)-G(i,i)*p(i)),0);
        
     NU(t+1,i)=NU(t,i)+1/t*((1-varsigma*Dmax)*I1(i)/(tau*Dmax*log(1-e2))-G(i,i)*P(t,i));%v的更新
     Fai(t+1,i)=Fai(t,i)+1/t*(rth*I1(i)/(log(Omax)*log(1-e3))-G(i,i)*P(t,i));
     
      B(t,i)=MU2(t)*X(2,i)+MU3(t)*X(3,i)+MU4(t)*X(4,i)+MU5(t)*X(5,i)
     +Lamda2(t,i)*sqrt(5)*SIGMA*alfai(2,i)+Lamda3(t,i)*sqrt(5)*SIGMA*alfai(3,i)
     +Lamda4(t,i)*sqrt(5)*SIGMA*alfai(4,i)+Lamda5(t,i)*sqrt(5)*SIGMA*alfai(5,i)
  A(t,i)=NU(t,i)*G(i,i)+Fai(t,i)*G(i,i)
 
 f(t,i)=-NU(t,i)*G(i,i)-Fai(t,i)*G(i,i)+MU2(t)*X(2,i)+MU3(t)*X(3,i)+MU4(t)*X(4,i)+MU5(t)*X(5,i)
     +Lamda2(t,i)*sqrt(5)*SIGMA*alfai(2,i)+Lamda3(t,i)*sqrt(5)*SIGMA*alfai(3,i)
    +Lamda4(t,i)*sqrt(5)*SIGMA*alfai(4,i)+Lamda5(t,i)*sqrt(5)*SIGMA*alfai(5,i);
     
     p(t,i)=I1(i)/G(i,i)
     
     a(t,i)=W/(F(t,i)*G(1,i)+f(t,i))
     

 end   

     P(t+1,:)=a(t,:)-p(t,:);
    
     Lamda2(t+1,:)=Lamda2(t,:)+1/t*(sqrt(5)*SIGMA*alfai(2,i).*P(t,:)+(P(t,:)*X'-Ith)/C);
     Lamda3(t+1,:)=Lamda3(t,:)+1/t*(sqrt(5)*SIGMA*alfai(3,i).*P(t,:)+(P(t,:)*X'-Ith)/C);
     Lamda4(t+1,:)=Lamda4(t,:)+1/t*(sqrt(5)*SIGMA*alfai(4,i).*P(t,:)+(P(t,:)*X'-Ith)/C);
     Lamda5(t+1,:)=Lamda5(t,:)+1/t*(sqrt(5)*SIGMA*alfai(5,i).*P(t,:)+(P(t,:)*X'-Ith)/C);
     %Fai(i+1,:)=Fai(i,:)+1/i*(0.9*(Ith+Delta)/log(0.9)/log(0.9)*ones(1,5)-exp(P(i,:)).*GG2)



%-------------------上层博弈------------------%
  for i=2:5
      K=rth/(rth*G(1,:)*G(:,1)/G(i,i)+G(1,1));
           b(t,i)= f(t,i)*G(i,i)*G(i,1);
     e(t,i)=K*f(t,i)*G(1,1)*G(i,1)*G(i,1); 
     d(t,i)=z*K*G(i,i)*G(i,1);
     u(t,i)=K*G(1,1)*G(i,1)*G(i,1)/G(i,i)-G(i,1);
    % j(t,i)=K*Delta*Delta*G(1,1)*G(i,1);
     g(t,i)=K*Delta*G(1,1)*G(i,1)-I1(i)*u(t,i);
     h(t,i)=f(t,i)/G(i,1)
     j(t,i)=sqrt(W*(b(t,i)-e(t,i)+d(t,i))/(g(t,i)))/G(i,1)
 end   
P(t+1,1)=K*(W./(F(t,:).*G(1,:)+f(t,:))-I1(i)/G(i,i))*G(:,1)+K*Delta;
  
%F(t+1,:)=sqrt(W*(b(t,:)-e(t,:)+d(t,:))/g(t,:)/G(:,1))- 0.1*h(t,:);
F(t+1,:)=j(t,:)-h(t,:);
 end

figure 
plot(P(:,1),'-og');
hold on
plot(P(:,2),'-+r');
plot(P(:,3),'-*b');
plot(P(:,4),'-sk');
plot(P(:,5),'-dc');

legend('CU','P1','P2','P3','P4')
xlabel('iteration');
ylabel('Power conversion value(W)');  
xlim([0 20])



figure 
plot(F(:,2),'-+r');
hold on
plot(F(:,3),'-*b');
plot(F(:,4),'-sk');
plot(F(:,5),'-dc');

legend('1','2','3','4')
xlabel('iteration');
ylabel('Interference price');  
xlim([0 20])