clc
close all
clear
V1=[20,32,22,38,30];%发送车的速度
V2=[0,0,0,0,0];%接收基站的速度
v=zeros(5,5);%车与基站之间的相对速度矩阵
h=zeros(5,5);%车与基站之间的先前时刻的值
D=zeros(3,1);
for i=1:5
    for j=1:5
        v(i,j)=abs(V1(i)-V2(j));
        h(i,j)=1+0.3*rand(1);
    end
end
h=[1.19329543905811,1.11358281479808,1.24347413748474,1.15984767663984,1.10521813107307;1.28170046859997,1.26278284344790,1.16504690286953,1.18674252580037,1.17611341135943;1.06232268781991,1.09037389908385,1.14127700455528,1.06914644806347,1.25329263780862;1.05842928687011,1.06777653429172,1.05121241414436,1.06829928934497,1.13070960523117;1.09333068599512,1.27701389263097,1.12906221739888,1.05544489603724,1.27146429060397];
for CC=1:3
    CE=1.5+0.3*CC;
    tt=[1,2,3];
    aqa=1.5+0.3*tt;
Delta=1e-6;%背景噪声
f5 = 5; f4 = f5; f3 = f4; f2 = f3; f1 = f2;
t_imax=0.25; %执行时候最大可以容忍的时间
Tc=0.001;%transmission latency between RSU to Cloud is set to fix
f_ba=5;%边缘的固定算力
d_up=0.5;%上传的数据量
c_exe=CE;%执行的数据量
tau=100;%数据包长度
f_total=30;%云分给边缘的计算资源
varsigma=200;%数据包流量强度
Dmax=0.2;%时延阈值
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
diagG=-diag(G1)/9;%SINR阈值设为9
G11=zeros(5,5);
 for i=1:5
     G11(i,i)=diagG(i,1);
 end
 for i=1:5
      G1(i,i)=0;
 end
 G1=G11+G1;
X=G1+MUk*alfak+betak;%xi部分
SIGMA=0.5;
E=0.1;%第一个阈值
H=X+SIGMA*sqrt(-2*log(E))*alfak; 
I=H;%贝恩斯坦近似约束矩阵
global N
N=50;%迭代次数
arfa=zeros(N,5);
SINR=zeros(N,5);
Pmax=-9*ones(N,5);
C1=linspace(1,1,5);
S=zeros(N,5);
s=zeros(N,5);
Ith=2e-6;%干扰阈值Gama要给相应的存储空间
Temp=zeros(N,5);%每次迭代收敛的最优功率暂存区 40应为Z
Z=1;%迭代次数每次加一
I1=zeros(N,5);
Total_T=zeros(N,5);
P=zeros(N,5);
Lamda=zeros(N,5);
Miu=zeros(N,5);
Ksai=zeros(N,5);
Fai=zeros(N,5);
fi=zeros(N,5);
sumfai=zeros(N,5);
sumfi=zeros(N,5);
u_iexe=zeros(N,5);
P(1,:)=linspace(-6.5,-6.5,5);
Lamda(1,:)=linspace(15,15,5);
Miu(1,:)=linspace(20,20,5);
for bcd=1:100
for k=1:N
    fi(k,:)=[f1,f2,f3,f4,f5];%这是每个云给边缘分配的算力
    u_iexe(k,:)=(t_imax*C1-(Tc+(c_exe./(f_ba*C1+fi(k,:)))))/t_imax*(1/d_up/log(2));%这是执行时间的效用多乘了一个上传前缀
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
               S(k,j)=-Miu(k,j)*tau*(log(1-e2)-G(j,j));
            end
        end
               s(k,i)=sum(S(k,:));
    end
%     P(k+1,:)=min(log(u_iexe(k,:).*arfa(k,:)/log(2))-log(s(k,:)),Pmax(k,:));
    P(k+1,:)=log(u_iexe(k,:).*arfa(k,:))-log(s(k,:));
    for i=1:5
         Lamda(k+1,i)=max(0,Lamda(k,i)+1000000/k*(exp(P(k,:))*G1(:,i)+Delta));
    end 
    for i=1:5
        A=(varsigma/tau+1/tau/(Dmax-c_exe/fi(k,i)));%z这是固定f后的常数项
        test1=power(2,A)*(Ith+Delta);
        Miu(k+1,i)=max(0,Miu(k,i)+100/k*(test1+exp(P(k,i))*(log(1-e2)-G(j,j))));%%这里需要换乘求导后的结果
    end  
end
for u=1:N
    fi(1,:)=[f1,f2,f3,f4,f5];
    Ksai(1,:)=linspace(3000,3000,5);
    Fai(1,:)=linspace(10,10,5);
    for i=1:5
        RRR=(exp(P(k,:))*G(i,i))/(I1(i)+Delta);
        Ri=(log(C1+RRR))/log(2);
        D1=C1./(tau*Ri-varsigma);
        BR=Ri.*c_exe/(t_imax*d_up); 
        BRR=BR+c_exe.*Ksai(u,:);
        sumfai(u,i)=sum(Fai(u,:));
        fi(u+1,:)=max(0,sqrt(BRR/sumfai(u,i))-f_ba);
    end    
   f1= fi(u+1,1);
   f2= fi(u+1,2);
   f3= fi(u+1,3);
   f4= fi(u+1,4);
   f5= fi(u+1,5);
   sumfi(u,1)=sum(fi(u,:));       
    for i=1:5       
        Ksai(u+1,i)=Ksai(u,i)+1000/u*(c_exe/(fi(u,i)+f_ba)+D1(1,i)-Dmax);
    end    
    for i=1:5       
        sumfi(u,i)=sum(fi(u,:));
        Fai(u+1,i)=max(0,Fai(u,i)+1/u*(sumfi(u,i)-f_total));
    end  
end    
end
 fi(1,:)=linspace(6,6,5);
D(CC)=Total_T(N);
end
data = [D(1,1), D(2,1), D(3,1)];
b = bar(data);
ch = get(b,'children');
set(ch,'FaceVertexCData',[4;2;3;1;5;6])
cc1=1000*aqa(1);
cc2=1000*aqa(2);
cc3=1000*aqa(3);
set(gca,'XTickLabel',{cc1,cc2,cc3})
ylabel('EE');
% R=zeros(N,5);
%     for i=1:5
%         R=(log(1+(exp(P(:,i))*G(j,j))/(Ith+Delta)))/log(2);  
%         U=1-(Tc+c_exe./(fi(:,i)+f_ba))/t_imax;
%         Total_T=R.*U/d_up;
%     end

% plot(exp(P(:,1)),'-+r');
% hold on
% plot(exp(P(:,2)),'-og');
% plot(exp(P(:,3)),'-*b');
% plot(exp(P(:,4)),'-sk');
% plot(exp(P(:,5)),'-dc');
% legend('P1','P2','P3','P4','P5')
% xlabel('iteration');
%  set(gca,'ylim',[0.00150,0.00185]); 
% ylabel('Power conversion value(W)');
% 
% figure
% plot((fi(:,1)),'-+r');
% hold on
% plot((fi(:,2)),'-og');
% plot((fi(:,3)),'-*b');
% plot((fi(:,4)),'-sk');
% plot((fi(:,5)),'-dc');
% legend('F1','F2','F3','F4','F5')
% xlabel('iteration');
% set(gca,'ylim',[5.987,6.011]); 
% % set(gca,'xlim',[20,30]);
% ylabel('computing');
% 
% figure
% plot((Total_T(:,1)),'-+r');
% legend('F1')
% xlabel('iteration');
% ylabel('EE');
