clc
close all
clear
T=60;%共有70个时隙
M=4;%车的数量
det=1;%solt length
TTT=T*det;%time span
lright=1;%车道索引车向右 right
lleft=-1;%车道索引车向左
Vvelocity = 15;  % 车车速度为40
V1=20*ones(1,M);%发送车的速度
H=80;%无人机高度
VUAV_max=30;%无人机速度约束
q_cvx=zeros(3,T);
UAVposition = zeros(T,3);
UAVposition(1,:) = [0;0;H];% 定义无人机初始位置
cvx_begin
cvx_solver sedumi
variable UAVposition(T,3)   
% for t=1:T
% VUAU_x(t)=30*rand();
% VUAU_y(t)=30*rand();
% end
[PLVR,PLVU]=vehiclePL(T,M,lright,lleft,Vvelocity,V1,det,VUAV_max,UAVposition);
%使用信噪比进行判断
X=zeros(M,T);
for t=1:T
    for m=1:M
    if  PLVU(t,m)<PLVR(t,m) %车车向基站通信为1
        X(t,m)=1;
    else
        X(t,m)=0;
    end;
    end
end
Delta=1e-6;%背景噪声是
C1=linspace(1,1,M);
%求解吞吐量第一个时隙的
P=zeros(T,M);
R=zeros(T,M);
for t=1:T
    for m=1:M
      P(t,m)=0.003;%功率的初值
      R(t,m)=X(t,m)*log(1+P(t,m)*PLVR(t,m)/Delta)+(1-X(t,m))*log(1+P(t,m)*PLVU(t,m)/Delta);
    end
end
%求解吞吐量和R SUM
for t=2:T
    RSUM(t)=sum(R(t,:));
end
%  expression Rate_N(1,N)
expression RSUM(T)
maximize RSUM(T)
subject to
for t=2:T-1
    norm([q_cvx(1,n) q_cvx(2,n)]-[q_cvx(1,n+1) q_cvx(2,n+1)])<=Vmax*delta;
end
cvx end