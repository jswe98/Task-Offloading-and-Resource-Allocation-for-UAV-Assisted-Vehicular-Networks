function [CARposition,l,GVR,distanceVR,Vposition] = CARinfo(T,M,det,roadlength,roadwidth)
Vposition = zeros(M,3,T);% ���峵���ĳ�ʼλ��
ll=zeros(1,M);
lright=1;%�������������� right
lleft=-1;%��������������
Vvelocity = 15;  % �����ٶ�Ϊ40
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
        Vposition_x=roadlength;
    end 
end
     Vposition(m,:,1)=[Vposition_x; roadwidth/2+roadwidth/2*ll(m); 0];
     for t=2:T
     Vposition(m,:,t)=[Vposition(m,1,t-1)+ll(m)*det*Vvelocity; roadwidth/2+roadwidth/2*ll(m); 0];
     end
end   % ���峵���ĳ�ʼλ��
RSU=zeros(M,3,T);
distanceVR=zeros(M,T);
for t=1:T
    for m=1:M
    RSU(m,:,t)=[0,0,0];
    distanceVR(m,t) = norm(RSU(m,:,t)-Vposition(m,:,t)); 
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
TC=0.0005;%��վ�ɼ���֮ͨ�ŵķ��䳵�����ŵ���·��CSI����(V2I)
c=3*1e+8;%����
fc=5.9*1e+9;%������Ƶ������Ƶ��
j1=2*pi*fc*TC/c*v;%0�ױ���������������V2I��
epsi=besselj(0,j1);%����������ֵ��V2I��er
l=5;%�����վͨ����˥��ϵ��
%�ȶ���Gvr�ĳ�ֵ�������ڵ�һ��ʱ϶�Ŀ�˥�����˥��
PLvr=zeros(M,1);
GVR=zeros(M,T);
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

CARposition=Vposition;
l=ll;


