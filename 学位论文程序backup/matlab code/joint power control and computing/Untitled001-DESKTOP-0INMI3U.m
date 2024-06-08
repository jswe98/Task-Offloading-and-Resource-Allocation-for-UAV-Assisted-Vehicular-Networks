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
Delta=1e-6;%背景噪声
tau=100;%数据包长度
varsigma=200;%数据包流量强度
Dmax=0.01;%时延阈值
e2=0.1;
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
MUk=0.5;
alfak=l*(1-epsi.^2);
betak=l*(1-epsi.^2);
X=G1+MUk*alfak+betak;%xi部分
SIGMA=0.5;
E=0.1;
H=X+SIGMA*sqrt(-2*log(E))*alfak; 
I=H;%贝恩斯坦近似约束矩阵
for i=1:5
    I(i,i)=0;
end%构造矩阵对角线为0
N=30;%迭代次数
arfa=zeros(N,5);
beta=zeros(N,5);
SINR=zeros(N,5);
Pmax=log(0.5)*ones(N,5);
C1=linspace(1,1,5);
S=zeros(N,5);
s=zeros(N,5);
Gama=0.1;%原来是
Ith=2e-6;%干扰阈值Gama要给相应的存储空间
Temp=zeros(N,5);%每次迭代收敛的最优功率暂存区 40应为Z
Z=1;%迭代次数每次加一
I1=zeros(N,5);
Total_T=zeros(N,5);
for sub1=1:N
P=zeros(N,5);
Lamda=zeros(N,5);
Miu=zeros(N,5);
P(1,:)=[-8,-8,-8,-8,-8]; 
Lamda(1,:)=linspace(20,20,5);
Miu(1,:)=linspace(20,20,5);
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
               S(k,j)=B/log(2)*arfa(k,j)*G(i,j)*SINR(k,j)/((exp(P(k,j))*G(j,j)))+Lamda(k,j)*I(i,j)+Miu(k,j)*G(i,j);%%到底是G(i,j)还是G(j,i),应该是G(i,j)。当j不变是第j列，代表的是各个CM向某个簇头传输的信息。
            else
               S(k,j)=0;
            end
        end
               s(k,i)=sum(S(k,:));
    end
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
    for i=1:5
        Miu(k+1,i)=2*Miu(k,i);
        n(k,i)=2*Miu(k,i)*5;
    end  
    if  max(abs(n(k,:)))<=3e-7    %该数值可改，表示乘子在一定的程度下就认为已经稳定了
        Miu(k+1,:)=Miu(k,:);
    end
end    

Z=Z+1;%表示一共迭代了多少次就收敛了
end   %sub1 end

