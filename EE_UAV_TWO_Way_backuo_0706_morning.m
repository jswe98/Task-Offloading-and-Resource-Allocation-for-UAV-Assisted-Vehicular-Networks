clc
close all
clear
format long  
% format short 
N=70;%共有70个时隙
M=4;%车的数量
det=1;%solt length
T=N*det;%time span
lright=1;%车道索引车向右 right
lleft=-1;%车道索引车向左
%UAV coordinate setup
UAVposition = zeros(N,3);
UAVposition(1,:) = [0;0;80];% 定义无人机初始位置
VUAV_max=30;%无人机速度约束
D=VUAV_max*det;%maximum distance of UAV
VUAU_x=30*rand();
VUAU_y=30*rand();
% 模拟无人机位置变化C1=linspace(1,1,N);
for t = 2:N
    UAVvelocity = [VUAU_x; VUAU_y];  % 假设无人机在每个时刻沿着固定方向移动，这里假设速度为[1;1;1]
    UAVvelocity = UAVvelocity / norm(UAVvelocity) * VUAV_max;    % 缩放方向向量以控制飞行距离小于30
    UAVposition(t,:) = [UAVposition(t-1,1)+ det.*UAVvelocity(1);UAVposition(t-1,2)+ det.*UAVvelocity(2) ;100];   % 更新位置
%     UAVposition(t,:) = [UAVposition(t-1,1:2) + det.*velocity(1:2), 100];
end
%{
figure;
plot3(UAVposition(1,:), UAVposition(2,:), UAVposition(3,:));
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('无人机轨迹');
%}
Vposition = zeros(4,3,N);% 定义车车的初始位置
% Vposition(:,1) = [0; 0; 0];% 定义车车的初始位置
Vvelocity = 15;  % 车车速度为40
ll=zeros(1,4);
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
     for t=2:N
     Vposition(m,:,t)=[Vposition(m,1,t-1)+ll(m)*det*Vvelocity; 2+2*ll(m); 0];
     end
end   % 定义车车的初始位置
% 定义车车与路边单元的距离
RSU=zeros(M,3,N);
distanceVR=zeros(M,1,N);
 % 定义车车与无人机的距离
distanceVU=zeros(M,1,N);
UAVpositions=zeros(M,3,N);
for t=1:N
    for m=1:M
    RSU(m,:,t)=[0,0,0];
    distanceVR(m,:,t) = norm(RSU(m,:,t)-Vposition(m,:,t)); 
    UAVpositions(m,:,t)=UAVposition(t,:);
    distanceVU(m,:,t) = norm(UAVpositions(m,:,t)-Vposition(m,:,t)); 
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
 %h=[1.17065638186208;1.15200009731163;1.21035565946212;1.15704945577761];
B=1;%带宽
f = 2.4e9;   % 频率（Hz）
T=0.0005;%基站采集与之通信的发射车干扰信道链路的CSI周期(V2I)
c=3*1e+8;%光速
fc=5.9*1e+9;%多普勒频移中心频率
j1=2*pi*fc*T/c*v;%0阶贝塞尔函数参数（V2I）
epsi=besselj(0,j1);%求贝塞尔函数值（V2I）er
l=5;%车与基站通信慢衰落系数
%先定义Gvr的初值，这是在第一个时隙的快衰落加慢衰落
PLvr=zeros(M,1);
GVR=zeros(M,1,N);
GVU=zeros(M,1,N);
for m=1:M
    PLvr(m,1,1) =l/(distanceVR(m,1,1)^2);%路径损耗的倒数
    GVR(m,1,1)=l*(PLvr(m,1,1))^2; 
end
for t=2:N
    for m=1:M
    PLvr(m,1,1) =l/(distanceVR(m,1,1)^2);%路径损耗的倒数
    GVR(m,1,t)=((epsi(m)*GVR(m,1,t-1))^2+(1-epsi(m))^2)*l/(distanceVR(m,1,t)^2);
    end
end
%先定义Gvu的初值
for t=1:N
    for m=1:M
    GVU(m,1,t) =l/(distanceVU(m,1,t)^2);%路径损耗的倒数
    end
end
%使用信噪比进行判断
X=zeros(M,N);
for t=1:N
    for m=1:M
    if  GVU(m,1,t)>GVR(m,1,t)
        X(m,t)=1;
    else
        X(m,t)=0;
    end;
    end
end
%{
% 定义车道的宽度和车的初始位置
lane_width = 3;
initial_positions = [1.5 0.5]; % 车的初始位置分别为 (1.5, 0) 和 (0.5, 0)
velocity = 40;  % 速度为40
% 模拟车的位置变化
for t = 2:N
    % 更新每辆车的位置
    Vposition(:, t) = Vposition(:, t-1) + velocity * dt;
end
for t=1:N
QUAV(t)=[VUAV_max*t,VUAV_max*t,100];
end
%vehicle coordinate setup
for i=1:4
    x_user(i)=6;
    y_user(i)=2;
end
x_user(1)=6;
x_user(2)=8;
x_user(3)=8;
y_user(1)=2;
y_user(2)=1.5;
y_user(3)=2;\
%}