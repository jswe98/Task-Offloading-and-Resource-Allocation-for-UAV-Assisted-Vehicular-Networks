function [PLVR,PLVU] = vehiclePL(T,M,lright,lleft,Vvelocity,V1,det,VUAV_max,H,q_cvx)%2b��ɢ�����ƽ�����ʣ�omega��LOS����
D=VUAV_max*det;%maximum distance of UAV
UAVposition = zeros(T,3);
UAVposition(1,:) = [0;0;H];% �������˻���ʼλ��
for t = 2:T
%     UAVvelocity = [VUAU_x(t); VUAU_y(t)];  % �������˻���ÿ��ʱ�����Ź̶������ƶ�����������ٶ�Ϊ[1;1;1]
q_cvx(:,t)=q_cvx(:,t-1)+[500/(N-1),500/(N-1),H]';
 %     UAVvelocity = UAVvelocity / norm(UAVvelocity) * D;    % ���ŷ��������Կ��Ʒ��о���С��30
%     UAVposition(t,:) = [UAVposition(t-1,1)+ det.*UAVvelocity(1);UAVposition(t-1,2)+ det.*UAVvelocity(2) ;100];   % ����λ��
%     UAVposition(T,:)=UAVposition(1,:);
    %     UAVposition(t,:) = [UAVposition(t-1,1:2) + det.*velocity(1:2), 100];
end

Vposition = zeros(M,3,T);% ���峵���ĳ�ʼλ��
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
V2=zeros(1,1);%���ջ�վ���ٶ�
v=zeros(M,1);%�����վ֮�������ٶȾ���
h=zeros(M,1);%�����վ֮�����ǰʱ�̵�ֵ
for m=1:M 
        v(m)=abs(V1(m)-V2);
        h(m)=1+0.3*rand(1);    
end
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
figure;
plot(UAVposition(1,:), UAVposition(2,:));
%plot3(UAVposition(1,:), UAVposition(2,:), UAVposition(3,:));
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('���˻��켣');
PLVR=GVR';%ShadowedRician distribution
PLVU=GVU';
