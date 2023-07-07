function [PLVR,PLVU] = vehiclePL(T,M,lright,lleft,Vvelocity,V1,det,VUAV_max,H,q_cvx)%2b是散射分量平均功率，omega是LOS分量
D=VUAV_max*det;%maximum distance of UAV
UAVposition = zeros(T,3);
UAVposition(1,:) = [0;0;H];% 定义无人机初始位置
for t = 2:T
%     UAVvelocity = [VUAU_x(t); VUAU_y(t)];  % 假设无人机在每个时刻沿着固定方向移动，这里假设速度为[1;1;1]
q_cvx(:,t)=q_cvx(:,t-1)+[500/(N-1),500/(N-1),H]';
 %     UAVvelocity = UAVvelocity / norm(UAVvelocity) * D;    % 缩放方向向量以控制飞行距离小于30
%     UAVposition(t,:) = [UAVposition(t-1,1)+ det.*UAVvelocity(1);UAVposition(t-1,2)+ det.*UAVvelocity(2) ;100];   % 更新位置
%     UAVposition(T,:)=UAVposition(1,:);
    %     UAVposition(t,:) = [UAVposition(t-1,1:2) + det.*velocity(1:2), 100];
end

Vposition = zeros(M,3,T);% 定义车车的初始位置
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
RSU=zeros(M,3,T);
distanceVR=zeros(M,T);
 % 定义车车与无人机的距离
distanceVU=zeros(M,T);
UAVpositions=zeros(M,3,T);

for t=1:T
    for m=1:M
    RSU(m,:,t)=[0,0,0];
    distanceVR(m,t) = norm(RSU(m,:,t)-Vposition(m,:,t)); 
    UAVpositions(m,:,t)=UAVposition(t,:);
    distanceVU(m,t) = norm(UAVpositions(m,:,t)-Vposition(m,:,t)); 
    end
end
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
%先定义Gvu的初值
for t=1:T
    for m=1:M
    GVU(m,t) =l/(distanceVU(m,t)^2);%路径损耗的倒数
    end
end
figure;
plot(UAVposition(1,:), UAVposition(2,:));
%plot3(UAVposition(1,:), UAVposition(2,:), UAVposition(3,:));
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('无人机轨迹');
PLVR=GVR';%ShadowedRician distribution
PLVU=GVU';
