function [Q] = rate(V1,p)
%-------------------系统参数-------------------%发送车速度32，30，34，34，30.相应的D2D对接收基站0，0，0，0，0
% V1=[32,30,34,34,30];%发送车CM速度
V2=[0,0,0,0,0];%接收基站BS速度
v=zeros(5,5);%车与基站之间的相对速度矩阵
h=zeros(5,5);%车与基站之间的先前时刻的值
% g=zeros(5,5);%误差，服从参数为1指数分布
% G=zeros(5,5);%车与车之间信道增益矩阵
I=zeros(5,5);%车与基站之间干扰信道增益矩阵（对角线全为0元素，因为对角线是有效信道增益）
for i=1:5
    for j=1:5
        v(i,j)=abs(V1(i)-V2(j));
        h(i,j)=1+0.3*rand(1);
    end
end
    
Delta=1e-6;%背景噪声
T=0.0005;%基站采集与之通信的发射车干扰信道链路的CSI周期(V2I)
c=3*1e+8;%光速
fc=5.9*1e+9;%多普勒频移中心频率
j1=2*pi*fc*T/c*v;%0阶贝塞尔函数参数（V2I）
epsi=besselj(0,j1);%求贝塞尔函数值（V2I）
a=(epsi.^2).*h;%车与基站通信移动信道快衰落的模型中的CSI反馈部分（先前时刻的信道状态信息）
 l=0.0001;%车与基站通信慢衰落系数
  %W=0.001;
G1=l*a;%V2I通信衰落（先前时刻的信道增益（也是测得的平均信道增益））
G2=l*(1-epsi.^2);%v2i通信信道衰落误差部分的平均信道增益
G=G1+G2;%车与基站通信总的平均信道增益矩阵（该矩阵用于求取平均信噪比，长期速率和）
  for i=1:5
      G(i,i)=100*G(i,i);%有效信道增益的倍数扩大（簇内通信）
  end
%------计算信噪比初始值------%  
%  p=[-6 -6 -6 -6 -6];
%  M=zeros(1,5);
%  Q=zeros(1,5);%逼近的beta参数初始赋值
sinr=zeros(1,5);
Q=zeros(1,5);%速率的初始赋值
ganrao=zeros(25,5);
 s=zeros(25,5);
 global N
 N=25;
 for k=1:N
 for i=1:5
            for j=1:5
            if j~=i
                 ganrao(k,j)=(exp(p(k,j))*G(j,i));
                  else
                 ganrao(k,j)=0;
            end
            end
        s(k,i)=sum(ganrao(k,:));
 end
 end
 

for i=1:5
  I(i)=G(i,:)*exp(p')+Delta-G(i,i)*exp(p(i));
 sinr(i)=G(i,i)*exp(p(i))/I(i);
 Q(i)=log(1+sinr(i));
%   M(i)= sinr(i)/(1+ sinr(i));
%   Q(i)=log(1+sinr(i))-M(i)*log(sinr(i));
 %MH=G(i,:)*exp(p');
 %Q=exp(p0)*H0(i);
end 
% for i=1:5
%         for j=1:5
%             if j~=i
%                 ganrao(k,j)=(exp(P(k,j))*G(j,i)
%                  else
%                 S(k,j)=0;
%             end
%         end
%         s(k)=sum(ganrao(k,:));
%                 f(k,i)=L(k)*H(i)+(1/log(2))*s(k);
%   
%   end

