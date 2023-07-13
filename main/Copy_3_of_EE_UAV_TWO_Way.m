clc
close all
clear
% format long  
format short 
T=60;%����70��ʱ϶
H=80;%���˻��߶�
Delta=1e-6;%����������
l=10;%�����վͨ����˥��ϵ��
M=4;%��������
UAVposition = zeros(T,3);
distanceVU=zeros(M,1,T);
UAVposition(1,:) = [0;0;H];% �������˻���ʼλ��
UAVposition(T,:)=[800,4,H];
UAVpositions=zeros(M,3,T);
for t = 2:T
    UAVposition(t,:)=UAVposition(t-1,:)+[800/(t-1),800/(t-1),0];
end
det=1;%solt length
TTT=T*det;%time span
lright=1;%�������������� right
lleft=-1;%��������������
VUAV_max=30;%���˻��ٶ�Լ��
D=VUAV_max*det;%maximum distance of UAV
Vposition = zeros(M,3,T);% ���峵���ĳ�ʼλ��
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
GVU=zeros(M,T);
P=zeros(T,M);
R=zeros(T,M);


cvx_solver SDPT3
cvx_begin
% cvx_solver SeDuMi
variable UAVpositionCVX(T,3)   
% variable distanceVU(M,T)  
expression R(T,M)
for t=1:T
    for m=1:M
%     UAVpositions(m,:,t)=square_pos(UAVpositionCVX(t,:));
%       distanceVU(m,t)=square_pos(((Vposition(m,1,t)-UAVpositionCVX(t,1)).^2+(Vposition(m,2,t)-UAVpositionCVX(t,2)).^2+(Vposition(m,3,t)-UAVpositionCVX(t,3)).^2));
      distanceVU(m,:,t) = square_pos(UAVpositionCVX(t,1).^2); 
     GVU(m,t) =l/(distanceVU(m,t)^2);%·����ĵĵ���
    end
end

for t=1:T
    for m=1:M
      P(t,m)=0.003;%���ʵĳ�ֵ
      R(t,m)=log(1+P(t,m)*GVU(m,t)'/Delta);
    end
end
for t=1:T
RSUM(t)=sum(R(t,:));
end
maximize sum(RSUM(t))
subject to

for t=1:T-1
    norm([UAVpositionCVX(t+1,1) UAVpositionCVX(t+1,2)]-[UAVpositionCVX(t,1) UAVpositionCVX(t,2)])<= D;
end
cvx_end
qqqq=UAVpositionCVX(T,:)
UAV=zeros(3,T);
dense_matrix = full(UAVpositionCVX);
% figure;
% plot3(UAV(1,:), UAV(2,:), UAV(3,:));
% grid on;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('���˻��켣');

