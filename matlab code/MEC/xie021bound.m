 clc
close all
clear
%-------------------系统参数-------------------%发送车速度32，30，34，34，30.相应的D2D对接收车32，30，34，34，30
V=[32,30,34,34,30];%发送车速度
v1=zeros(5,5);%车与车之间的相对速度矩阵
v2=zeros(1,5);%车与基站之间的相对速度矩阵
h1=zeros(5,5);%车与车之间的先前时刻的值
h2=zeros(1,5);%车与基站之间的先前时刻的值
g1=zeros(5,5);%误差，服从参数为1指数分布
g2=zeros(1,5);%误差，服从参数为1指数分布
 G= zeros(5,5);%车与车之间信道增益矩阵
for i=1:5
    v2(i)=rand(1)*V(i);
    h2(i)=0.6+0.4*rand(1);
    g2(i)=min(exprnd(1),1);
    for j=1:5
        v1(i,j)=abs(V(i)-V(j));
        h1(i,j)=0.6+0.4*rand(1);
        g1(i,j)=min(exprnd(1),1);
    end
end
g1=[1,1,1,1,0.193944467428701;0.247071917326137,0.00667900504869646,0.429831762570529,1,1;0.499989141185540,0.912457204764465,0.0828545933129404,1,1;0.316515884647952,1,0.147265296729244,0.888397941568164,0.287355068625972;0.692602227188954,0.0236457740563112,0.914085998106791,1,0.0363978204503782];
 h1=[0.707376063316427,0.792403814371603,0.619440919296757,0.703378523100252,0.842277225673208;0.943751141721211,0.681518740418938,0.637449811287940,0.686070610446946,0.735875122029881;0.833785936127243,0.885857146835165,0.943476408775190,0.900335555797250,0.918726315155221;0.753335972407940,0.954913872011155,0.655286537026870,0.768698593019774,0.983656536227459;0.753832879687560,0.828102347548962,0.797142072406048,0.997996104728119,0.866129884296064];
 h2=[0.843523343836463,0.920894806411973,0.980303193098356,0.844275180592760,0.693406322250299];
 g2=[0.734149023801544,0.345151061941124,0.929864392850716,0.351265692707954,1];
v2=[31.0582955144346,11.9631134127137,9.90779637426617,15.4771414293975,25.5544585676781];
% v1=[0,2,2,2,2;2,0,4,4,0;2,4,0,0,4;2,4,0,0,4;2,0,4,4,0];%发射汽车到接收汽车之间的相对速度设置为一样的，可以构建一个车与车之间的相对速度矩阵.
% v2=[25,22,20,20,23];%发射车到基站之间的相对速度，为车的速度乘以一个cosa
T1=0.0001;%车与车移动信道CSI采样周期T1<T2,T2为车到基站移动信道的采样CSI（V2V）
T2=0.0005;%基站采集与之通信的发射车干扰信道链路的CSI周期(V2I)
c=3*1e+8;%光速
fc=5.9*1e+9;%多普勒频移中心频率
j1=2*pi*fc*T1/c*v1;%0阶贝塞尔函数参数（V2V）
j2=2*pi*fc*T2/c*v2;%0阶贝塞尔函数参数（V2I）
epsi1=besselj(0,j1);%求贝塞尔函数值（V2V）
epsi2=besselj(0,j2);%求贝塞尔函数值（V2I）
a1=(epsi1.^2).*h1+(1-epsi1.^2).*g1;%车车通信移动信道快衰落的模型
a2=(epsi2.^2).*h2;%（V2I）通信移动信道快衰落的模型中的CSI反馈部分（先前时刻的信道状态信息）
  l1=0.001;%车车通信慢衰落系数
  l2=0.005;%V2I慢衰落系数
  G=l1*a1;%车车通信衰落
  G2=l2*a2;%
  Delta=1e-6;%背景噪声
  for i=1:5
      G(i,i)=20*G(i,i);
  end
  d=[290 250 320 350 400];
  H0=l1*d.^(-3);%蜂窝用户到蜂窝内小汽车的距离
  p0=-1;%设置的蜂窝用户的发射功率 
  MUk=0.5;
  alfak=l2*sqrt(1-epsi2.^2);
  betak=l2*sqrt(1-epsi2.^2);
 X=G2+MUk*alfak+betak ;%%d2d-v用户到蜂窝内基站的先前的信道增益
  SIGMA=0.5;
   D=zeros(5,1);
D1=zeros(5,5);
  
  for Z=1:5
   E=Z/10; %0.1 [12.0483178226206]; 0.5，`[12.0978389067218]  0.01,[12.0011555595393] 0.05, [12.0324329381173]
   H=X+SIGMA*sqrt(-2*log(E))*alfak; 
global N
N =15;
arfa=zeros(N,5);
beta=zeros(N ,5);
SINR=zeros(N,5);
P=zeros(N,5);
P(1,:)=[-9 -9 -9 -9 -9];
Pmax=log(0.5)*ones(N,5);
L=zeros(N,1);
L(1)=25;
 W=0.001;
 S=zeros(N,5);
 s=zeros(N,1);
 Ith=1e-4;%干扰阈值
 Total_T=zeros(N,1);
C=linspace(1,1,5);
for k=1:N
  for i=1:5
     %  w=exp(P(k,:)')H;观测量
        I1(i)=G(i,:)*exp(P(k,:)')+exp(p0)*H0(i)+Delta-G(i,i)*exp(P(k,i));
        SINR(k,i)=G(i,i)*exp(P(k,i))/I1(i); 
           arfa(k,i)=SINR(k,i)/(1+SINR(k,i));
  end
   for i=1:5
        for j=1:5
            if j~=i
                S(k,j)=W*arfa(k,j)*G(j,i)*SINR(k,j)/((exp(P(k,j))*G(j,j)));
                 else
                S(k,j)=0;
            end
        end
        s(k)=sum(S(k,:));
        f(k,i)=L(k)*H(i)+(1/log(2))*s(k);
  
  end
   Total_T(k)=dot(log(C+SINR(k,:))/log(2),C);
P(k+1,:)=log(W.*arfa(1,:)/log(2))-log(f(k,:));
 L(k+1)=L(k)+10000/k*(exp(P(k,:))*H'-Ith);
w(k)=exp(P(k,:))*H'-Ith;
if  abs(w(k))<=1e-5
     L(k+1)= L(k);
end

end
D1(Z,:)=exp(P(N,:));
D(Z)=Total_T(N);
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
plot(L,'-+');
legend('L')
xlabel('iteration');
ylabel('Multiplier convergence');

 figure 
      plot(Total_T,'-*r');
      legend('sum rate')
xlabel('iteration');
ylabel('convergence');
