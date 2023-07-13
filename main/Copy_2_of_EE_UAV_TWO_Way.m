clc
close all
clear
% format long  
format short 
T=60;%共有70个时隙
H=40;%无人机高度
M=6;%车的数量
% UAVposition = zeros(T,3);
UAVposition(1,:) = [0;0;H];% 定义无人机初始位置
UAVpositions=zeros(M,3,T);
det=1;%solt length
TTT=T*det;%time span
lright=1;%车道索引车向右 right
lleft=-1;%车道索引车向左
%UAV coordinate setup
VUAV_max=30;%无人机速度约束
D=VUAV_max*det;%maximum distance of UAV
% 定义车车与无人机的距离
Vposition = zeros(M,3,T);% 定义车车的初始位置
Vvelocity = 15;  % 车车速度为40
ll=zeros(1,M);
for m = 1:M
random_number = randi([0, 1]); % 生成0到1之间的随机整数（0代表-1，1代表1）
if random_number == 0
    ll(m) = lright;   % 右车
    random_number = randi([0, 1]);
    if random_number == 0
        Vposition_x=0;
    else
        Vposition_x=200;
    end        
else
    ll(m) = lleft;
    random_number = randi([0, 1]);
    if random_number == 0
        Vposition_x=400;
    else
        Vposition_x=800;
    end 
end
     Vposition(m,:,1)=[Vposition_x; 2+2*ll(m); 0];
     for t=2:T
     Vposition(m,:,t)=[Vposition(m,1,t-1)+ll(m)*det*Vvelocity; 2+2*ll(m); 0];
     end
end   % 定义车车的初始位置
% 定义车车与路边单元的距离
RSU=zeros(M,3,T);
distanceVR=zeros(M,T);

for t=1:T
    for m=1:M
    RSU(m,:,t)=[0,0,0];
    distanceVR(m,t) = norm(RSU(m,:,t)-Vposition(m,:,t)); 
    end
end
distanceVU=zeros(M,1,T);

V1=20*ones(1,M);%发送车的速度
V2=zeros(1,1);%接收基站的速度
v=zeros(M,1);%车与基站之间的相对速度矩阵
h=zeros(M,1);%车与基站之间的先前时刻的值
for m=1:M 
        v(m)=abs(V1(m)-V2);
        h(m)=1+0.3*rand(1);    
end
 %h=[1.17065638186208;1.15200009731163;1.21035565946212;1.15704945577761];
B=1;%带宽
f = 2.4e9;   % 频率（Hz）
TC=0.0005;%基站采集与之通信的发射车干扰信道链路的CSI周期(V2I)
c=3*1e+8;%光速
fc=5.9*1e+9;%多普勒频移中心频率
j1=2*pi*fc*TC/c*v;%0阶贝塞尔函数参数（V2I）
epsi=besselj(0,j1);%求贝塞尔函数值（V2I）er
l=5;%车与基站通信慢衰落系数
%先定义Gvr的初值，这是在第一个时隙的快衰落加慢衰落
PLvr=zeros(M,1);
GVR=zeros(M,T);
GVU=zeros(M,T);
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

Delta=1e-6;%背景噪声是
C1=linspace(1,1,M);
%求解吞吐量第一个时隙的
P=zeros(T,M);
R=zeros(T,M);
for t = 2:T
%     UAVposition(t,:)=[UAVposition(1,1)+t*abs(UAVposition(T,1)-UAVposition(1,1))/T,UAVposition(1,2)+t*abs(UAVposition(T,2)-UAVposition(1,2))/T,H];
%     UAVposition(t,:)=UAVposition(t-1,:)+[500/(t-1),500/(t-1),0];
    UAVposition(t,:)=UAVposition(t-1,:)+[500/(t-1),500/(t-1),0];
end
for t=1:T
    for m=1:M
%     UAVpositions(m,:,t)=UAVposition(t,:);
      distanceVU(m,t)=sqrt((Vposition(m,1,t)-UAVposition(t,1)).^2+(Vposition(m,2,t)-UAVposition(t,2)).^2+(Vposition(m,3,t)-UAVposition(t,3)).^2);
%       distanceVU(m,:,t) = norm(UAVpositions(1,:,t)-Vposition(m,:,t)); 
    end
end

cvx_begin
cvx_solver SDPT3
variable UAVposition(T,3)   
expression R(T,M)
sumr=sum(R(:))
%先定义Gvu的初值
for t=1:T
    for m=1:M
    GVU(m,t) =l/(distanceVU(m,t)^2);%路径损耗的倒数
    end
end
%使用信噪比进行判断
X=zeros(M,T);
for t=1:T
    for m=1:M
    if  GVU(m,t)<GVR(m,t) %车车向基站通信为1
        X(m,t)=1;
    else
        X(m,t)=0;
    end;
    end
end
for t=1:T
    for m=1:M
      P(t,m)=0.003;%功率的初值
      R(t,m)=X(m,t)'*log(1+P(t,m)*GVR(m,t)'/Delta)+(1-X(m,t)')*log(1+P(t,m)*GVU(m,t)'/Delta);
    end
end
maximize sumr
subject to

for t=2:T-1
    norm([UAVposition(t+1,1) UAVposition(t+1,2)]-[UAVposition(t,1) UAVposition(t,2)])<= D;
end

cvx_end
sumr=sum(R(:))
qqqq=UAVposition(T,:)
dense_matrix = full(UAVposition);
UAV=UAVposition'
figure;
plot(UAV(1,:), UAV(2,:));
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('无人机轨迹');

% distanceVR=distanceVR';
% distanceVU=distanceVU';
% GVR=GVR';
% GVU=GVU';
% X=X';
