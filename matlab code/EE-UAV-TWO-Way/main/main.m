clc
close all
clear
T=60;%����70��ʱ϶
M=4;%��������
det=1;%solt length
TTT=T*det;%time span
lright=1;%�������������� right
lleft=-1;%��������������
Vvelocity = 15;  % �����ٶ�Ϊ40
V1=20*ones(1,M);%���ͳ����ٶ�
H=80;%���˻��߶�
VUAV_max=30;%���˻��ٶ�Լ��
q_cvx=zeros(3,T);
UAVposition = zeros(T,3);
UAVposition(1,:) = [0;0;H];% �������˻���ʼλ��
cvx_begin
cvx_solver sedumi
variable UAVposition(T,3)   
% for t=1:T
% VUAU_x(t)=30*rand();
% VUAU_y(t)=30*rand();
% end
[PLVR,PLVU]=vehiclePL(T,M,lright,lleft,Vvelocity,V1,det,VUAV_max,UAVposition);
%ʹ������Ƚ����ж�
X=zeros(M,T);
for t=1:T
    for m=1:M
    if  PLVU(t,m)<PLVR(t,m) %�������վͨ��Ϊ1
        X(t,m)=1;
    else
        X(t,m)=0;
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
      R(t,m)=X(t,m)*log(1+P(t,m)*PLVR(t,m)/Delta)+(1-X(t,m))*log(1+P(t,m)*PLVU(t,m)/Delta);
    end
end
%�����������R SUM
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