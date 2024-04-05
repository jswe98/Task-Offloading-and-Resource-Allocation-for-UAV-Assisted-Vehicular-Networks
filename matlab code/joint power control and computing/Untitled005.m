clc
close all
clear
V1=[30,32,28,34,30];%发送车的速度
V2=[0,0,0,0,0];%接收基站的速度
v=zeros(5,5);%车与基站之间的相对速度矩阵
h=zeros(5,5);%车与基站之间的先前时刻的值
for i=1:5
    for j=1:5
        v(i,j)=abs(V1(i)-V2(j));
        h(i,j)=1+0.3*rand(1);
    end
end
%%h=[1.21026765186132,1.00186773622345,1.11230368038275,1.27044868578780,1.09550344842318;1.17912487696444,1.08933855474170,1.03750431438572,1.11650670059712,1.24530632319216;1.29435271127300,1.25859695058615,1.02514626462880,1.10131359710114,1.07083874991164;1.09534161637965,1.29533452164231,1.16447526737693,1.22477537780633,1.25255552041723;1.05006691209611,1.27092928278688,1.03153723834667,1.22352792893222,1.21881158135137];
h=[1.02789667806120,1.13904677432867,1.00279975360833,1.27450777352074,1.19282252173993;1.00042571746072,1.00911558174848,1.06254106699160,1.13648984350063,1.03817980999754;1.00259430636343,1.21812388034456,1.10623493981698,1.23413378382881,1.13099699249059;1.13096643471168,1.01476395410695,1.01488957045582,1.02733005262392,1.17821110943324;1.07232521655071,1.25241073059164,1.25716382922718,1.28908366028541,1.14666993584818];
Delta=1e-6;%背景噪声
f5 = 1; f4 = f5; f3 = f4; f2 = f3; f1 = f2;
t_imax=0.3+4*log(5); %执行时候最大可以容忍的时间
Tc=1;%transmission latency between RSU to Cloud is set to fix
bate1=0.2;%网络边缘的系数
bate2=0.8;%云端的系数
f_ba=2;%边缘的固定算力
d_up=1;%上传的数据量
c_exe=1;%执行的数据量
tau=100;%数据包长度
varsigma=200;%数据包流量强度
Dmax=0.1;%时延阈值
e2=0.1;%第二个阈值
B=1;%带宽
T=0.0005;%基站采集与之通信的发射车干扰信道链路的CSI周期(V2I)
c=3*1e+8;%光速
fc=5.9*1e+9;%多普勒频移中心频率
j1=2*pi*fc*T/c*v;%0阶贝塞尔函数参数（V2I）
epsi=besselj(0,j1);%求贝塞尔函数值（V2I）
a=(epsi.^2).*h;%车与基站通信移动信道快衰落的模型中的CSI反馈部分（先前时刻的信道状态信息）
l=0.0001;%车与基站通信慢衰落系数
G1=l*a;%V2I通信衰落（先前时刻的信道增益（也是测得的平均信道增益））
G2=l*(1-epsi.^2);%v2i通信信道衰落误差部分的平均信道增益
G=G1+G2;%车与基站通信总的平均信道增益矩阵（该矩阵用于求取平均信噪比，长期速率和）
  for i=1:5
      G(i,i)=100*G(i,i);%有效信道增益的倍数扩大（簇内通信）
  end 
MUk=0.5;
alfak=l*(1-epsi.^2);
betak=l*(1-epsi.^2);
X=G1+MUk*alfak+betak;%xi部分
SIGMA=0.5;
E=0.1;%第一个阈值
H=X+SIGMA*sqrt(-2*log(E))*alfak; 
I=H;%贝恩斯坦近似约束矩阵
for i=1:5
    I(i,i)=0;
end%构造矩阵对角线为0
global N
N=90;%迭代次数
arfa=zeros(N,5);
beta=zeros(N,5);
SINR=zeros(N,5);
Pmax=(5*log(0.5)-1)*ones(N,5);
C1=linspace(1,1,5);
GF=linspace(1.5,1.5,5);%%%%%%%%%%
S=zeros(N,5);
s=zeros(N,5);
Gama=0.1;%原来是
Ith=2e-6;%干扰阈值Gama要给相应的存储空间
Temp=zeros(N,5);%每次迭代收敛的最优功率暂存区 40应为Z
Z=1;%迭代次数每次加一
I1=zeros(N,5);
Total_T=zeros(N,5);
% for sub1=1:N
% while 1
P=zeros(N,5);
Lamda=zeros(N,5);
Miu=zeros(N,5);
fi=zeros(N,5);
u_iexe=zeros(N,5);
P(1,:)=[-8,-8,-8,-8,-8]; 
Lamda(1,:)=linspace(30,30,5);
Miu(1,:)=linspace(50,50,5);
for k=1:N
    fi(k,:)=[f1,f2,f3,f4,f5];%这是每个云给边缘分配的算力
    u_iexe(k,:)=(t_imax*C1-((Tc+bate1*c_exe/f_ba)*C1+(bate2*c_exe)./fi(k,:)))/t_imax*(B/d_up/log(2));%
    for i=1:5
    I1(i)=exp(P(k,:))*G(:,i)+Delta-G(i,i)*exp(P(k,i));
    SINR(k,i)=G(i,i)*exp(P(k,i))/I1(i);      
    arfa(k,i)=SINR(k,i)/(1+SINR(k,i));
    Total_T(k)=dot(log(C1+SINR(k,:)),u_iexe(k,:)); %这是第一个子问题的和的形式
    end
    for i=1:5
        for j=1:5
            if j~=i
               S(k,j)=u_iexe(k,j)*arfa(k,j)*G(i,j)*SINR(k,j)/((exp(P(k,j))*G(j,j)))+Lamda(k,j)*I(i,j);%%到底是G(i,j)还是G(j,i),应该是G(i,j)。当j不变是第j列，代表的是各个CM向某个簇头传输的信息。
            else
               S(k,j)=-Miu(k,j)*G(j,j);
            end
        end
               s(k,i)=sum(S(k,:));
    end
%%P(k+1,:)=min(log(ZZ.*arfa(k,:)/log(2))-log(s(k,:)),Pmax(k,:));%%这里还需要改,这是功率p的迭代
P(k+1,:)=log(u_iexe(k,:).*arfa(k,:)/log(2))-log(s(k,:));
%     if max(abs(exp(P(k+1,:))-exp(P(k,:))))<=1e-5
%        P(k+1,:)=P(k,:);
%     end  
    for i=1:5
         Lamda(k+1,i)=max(0,Lamda(k,i)+100000/k*(exp(P(k,:))*I(:,i)-Ith));
    end 
%     if  max(abs(exp(P(k,:))*I(:,i)-Ith))<=3e-7    %该数值可改，表示乘子在一定的程度下就认为已经稳定了
%         Lamda(k+1,:)=Lamda(k,:);
%     end
    for i=1:5
        Miu(k+1,i)=max(0,Miu(k,i)+100000/k*(-exp(P(k,i))*G(i,i)+(exp(log(2)*(1+tau*varsigma)/(tau*(Dmax-c_exe./fi(k,i))))-1)*(exp(P(k,:))*I(:,i))/log(1/(1-e2))));%%这里需要换乘求导后的结果
    end  
%     if  max(abs((exp(log(2)*(1+tau*varsigma)/(tau*(Dmax-c_exe./fi(k,i))))-1)*(Ith+Delta)/log(1/(1-e2))-exp(P(k,:))*G(:,i)))<=3e-7    %该数值可改，表示乘子在一定的程度下就认为已经稳定了
%         Miu(k+1,:)=Miu(k,:);
%     end
end    
Z=Z+1;%表示一共迭代了多少次就收敛了
% end   %sub1 end
% aaaaaaaabaaaaa=G(i,i)*exp(P(k,i))