function [GVR,GVU,distanceVR,distanceVU,CARposition,L] = CARinfo(H,T,M,det,roadlength,roadwidth,UAVcvx)
Vposition = zeros(M,3,T);% 定义车车的初始位置
ll=zeros(1,M);
lright=1;%车道索引车向右 right
lleft=-1;%车道索引车向左
Vvelocity = 10;  % 车车速度为40

Vvelocity=zeros(1,M);  % 每辆车的速度都不一样
for i=1:M
    Vvelocity(i)=randi([10, 20]);
end

while 1
for m = 1:M
random_number = randi([0, 1]); % 生成0到1之间的随机整数（0代表-1，1代表1）
if random_number == 0
    ll(m) = lright;   % 向右跑的车
    random_number = randi([0, 1]);
    if random_number == 0
        Vposition_x=0;   
    else
        Vposition_x=roadlength/2;
    end   
    random_numbery = (roadwidth/2) * rand;  % 生成4到8之间的随机数
else
    ll(m) = lleft; %向左跑的车
    random_number = randi([0, 1]);
    if random_number == 0
        Vposition_x=roadlength/2;
    else
        Vposition_x=roadlength/4;
    end 
    random_numbery = roadwidth/2+(roadwidth/2) * rand; % 生成0到4之间的随机数
end
Vposition_x=randi([200,roadlength-200]);  %随机车辆的起点，之前的未使用
Vposition_x=randi([200,roadlength+200]);  %随机车辆的起点，之前的未使用
%      Vposition(m,:,1)=[Vposition_x; roadwidth/2+roadwidth/2*ll(m); 0];
     Vposition(m,:,1)=[Vposition_x; random_numbery+roadwidth/2*ll(m); 0];%第二个是车辆的上下位置
     for t=2:T
%      Vposition(m,:,t)=[Vposition(m,1,t-1)+ll(m)*det*Vvelocity; roadwidth/2+roadwidth/2*ll(m); 0];
     Vposition(m,:,t)=[Vposition(m,1,t-1)+ll(m)*det*Vvelocity(m); Vposition(m,2); 0];
     end
end
      if  sum(ll == 1)==0.5*M;
          break; 
          elseif sum(ll == 1)==0.5*M+0.5;
          break; 
          elseif sum(ll == 1)==0.5*M-0.5;
          break;     
      end 
end


Varin=load ('八十秒三千米.mat');%这两行都要一起用，其中CARposition是默认的不能改动
Vposition = Varin.CARposition;  %调用了之前的数据S

RSU=zeros(M,3,T);
distanceVR=zeros(M,T);
for t=1:T
    for m=1:M
    RSU(m,:,t)=[0,0,0];
    distanceVR(m,t) = norm(RSU(m,:,t)-Vposition(m,:,t)); 
    end
end
V1=20*ones(1,M);%发送车的速度
V2=zeros(1,1);%接收基站的速度
v=zeros(M,1);%车与基站之间的相对速度矩阵
h=zeros(M,1);%车与基站之间的先前时刻的值
for m=1:M 
        v(m)=abs(V1(m)-V2);
        h(m)=1+0.3*rand(1);    
end
TC=0.0005;%基站采集与之通信的发射车干扰信道链路的CSI周期(V2I)
c=3*1e+8;%光速
fc=5.9*1e+9;%多普勒频移中心频率
j1=2*pi*fc*TC/c*v;%0阶贝塞尔函数参数（V2I）
epsi=besselj(0,j1);%求贝塞尔函数值（V2I）er
l=5;%车与基站通信慢衰落系数
%先定义Gvr的初值，这是在第一个时隙的快衰落加慢衰落
PLvr=zeros(M,1);
GVR=zeros(M,T);
for m=1:M
    PLvr(m,1) =l/(distanceVR(m,1)^2);%路径损耗的倒数
    GVR(m,1)=l*(PLvr(m,1))^2; 
end
for t=2:T
    for m=1:M
    PLvr(m,1) =l/(distanceVR(m,1)^2);%路径损耗的倒数
    GVR(m,t)=((epsi(m)*GVR(m,t-1))^2+(1-epsi(m))^2)*l/(distanceVR(m,t)^2);%前一时刻的相关系数加上新的慢衰落
    end
end

UAVpositions=zeros(M,3,T);
GVU=zeros(M,T);
% 创建要添加的一列
newColumn = linspace(H,H,T)';
UAVposition= [UAVcvx, newColumn];
distanceVU=zeros(M,T);
for t=1:T
    for m=1:M
      UAVpositions(m,:,t)=UAVposition(t,:);
%       distanceVU(m,t)=sqrt((Vposition(m,1,t)-UAVposition(t,1)).^2+(Vposition(m,2,t)-UAVposition(t,2)).^2+(Vposition(m,3,t)-UAVposition(t,3)).^2);
      distanceVU(m,t) = norm(UAVpositions(1,:,t)-Vposition(m,:,t)); 
    end
end
for t=1:T
    for m=1:M
    GVU(m,t) =l/(distanceVU(m,t)^2);%路径损耗的倒数
    end
end

CARposition=Vposition;
L=ll;


