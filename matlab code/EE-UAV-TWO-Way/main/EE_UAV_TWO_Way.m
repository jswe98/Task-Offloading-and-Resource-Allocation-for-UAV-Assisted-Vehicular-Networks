clc
close all
clear
format long  
% format short 
T=60;%����70��ʱ϶
M=4;%��������
det=1;%solt length
TTT=T*det;%time span
lright=1;%�������������� right
lleft=-1;%��������������
%UAV coordinate setup
H=80;
UAVposition = zeros(T,3);
UAVposition(1,:) = [0;0;H];% �������˻���ʼλ��
VUAV_max=30;%���˻��ٶ�Լ��
D=VUAV_max*det;%maximum distance of UAV
VUAU_x=30*rand();
VUAU_y=30*rand();
for t=1:T
VUAU_x(t)=30*rand();
VUAU_y(t)=30*rand();
end
% ģ�����˻�λ�ñ仯C1=linspace(1,1,N);
for t = 2:T
    UAVvelocity = [VUAU_x(t); VUAU_y(t)];  % �������˻���ÿ��ʱ�����Ź̶������ƶ�����������ٶ�Ϊ[1;1;1]
    UAVvelocity = UAVvelocity / norm(UAVvelocity) * D;    % ���ŷ��������Կ��Ʒ��о���С��30
    UAVposition(t,:) = [UAVposition(t-1,1)+ det.*UAVvelocity(1);UAVposition(t-1,2)+ det.*UAVvelocity(2) ;H];   % ����λ��
%     UAVposition(T,:)=UAVposition(1,:);
    %     UAVposition(t,:) = [UAVposition(t-1,1:2) + det.*velocity(1:2), 100];
end

figure;
plot3(UAVposition(1,:), UAVposition(2,:), UAVposition(3,:));
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('���˻��켣');

Vposition = zeros(M,3,T);% ���峵���ĳ�ʼλ��
% Vposition(:,1) = [0; 0; 0];% ���峵���ĳ�ʼλ��
Vvelocity = 15;  % �����ٶ�Ϊ40
ll=zeros(1,M);
for m = 1:M
random_number = randi([0, 1]); % ����0��1֮������������0����-1��1����1��
if random_number == 0
    ll(m) = lright;   % �ҳ�
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
end   % ���峵���ĳ�ʼλ��
% ���峵����·�ߵ�Ԫ�ľ���
RSU=zeros(M,3,T);
distanceVR=zeros(M,T);
 % ���峵�������˻��ľ���
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
V1=20*ones(1,M);%���ͳ����ٶ�
V2=zeros(1,1);%���ջ�վ���ٶ�
v=zeros(M,1);%�����վ֮�������ٶȾ���
h=zeros(M,1);%�����վ֮�����ǰʱ�̵�ֵ
for m=1:M 
        v(m)=abs(V1(m)-V2);
        h(m)=1+0.3*rand(1);    
end
 %h=[1.17065638186208;1.15200009731163;1.21035565946212;1.15704945577761];
B=1;%����
f = 2.4e9;   % Ƶ�ʣ�Hz��
TC=0.0005;%��վ�ɼ���֮ͨ�ŵķ��䳵�����ŵ���·��CSI����(V2I)
c=3*1e+8;%����
fc=5.9*1e+9;%������Ƶ������Ƶ��
j1=2*pi*fc*TC/c*v;%0�ױ���������������V2I��
epsi=besselj(0,j1);%����������ֵ��V2I��er
l=5;%�����վͨ����˥��ϵ��
%�ȶ���Gvr�ĳ�ֵ�������ڵ�һ��ʱ϶�Ŀ�˥�����˥��
PLvr=zeros(M,1);
GVR=zeros(M,T);
GVU=zeros(M,T);
for m=1:M
    PLvr(m,1) =l/(distanceVR(m,1)^2);%·����ĵĵ���
    GVR(m,1)=l*(PLvr(m,1))^2; 
end
for t=2:T
    for m=1:M
    PLvr(m,1) =l/(distanceVR(m,1)^2);%·����ĵĵ���
    GVR(m,t)=((epsi(m)*GVR(m,t-1))^2+(1-epsi(m))^2)*l/(distanceVR(m,t)^2);%ǰһʱ�̵����ϵ�������µ���˥��
    end
end
%�ȶ���Gvu�ĳ�ֵ
for t=1:T
    for m=1:M
    GVU(m,t) =l/(distanceVU(m,t)^2);%·����ĵĵ���
    end
end
%ʹ������Ƚ����ж�
X=zeros(M,T);
for t=1:T
    for m=1:M
    if  GVU(m,t)<GVR(m,t) %�������վͨ��Ϊ1
        X(m,t)=1;
    else
        X(m,t)=0;
    end;
    end
end
Delta=1e-6;%����������
C1=linspace(1,1,M);
%�����������һ��ʱ϶��
P=zeros(T,M);
R=zeros(T,M);
for t=1:T
    for m=1:M
      P(t,m)=0.003;%���ʵĳ�ֵ
      R(t,m)=X(m,t)'*log(1+P(t,m)*GVR(m,t)'/Delta)+(1-X(m,t)')*log(1+P(t,m)*GVU(m,t)'/Delta);
    end
end
%�����������R SUM
RSUM=sum(R(:));



distanceVR=distanceVR';
distanceVU=distanceVU';
GVR=GVR';
GVU=GVU';
X=X';

%{
% ���峵���Ŀ��Ⱥͳ��ĳ�ʼλ��
lane_width = 3;
initial_positions = [1.5 0.5]; % ���ĳ�ʼλ�÷ֱ�Ϊ (1.5, 0) �� (0.5, 0)
velocity = 40;  % �ٶ�Ϊ40
% ģ�⳵��λ�ñ仯
for t = 2:N
    % ����ÿ������λ��
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