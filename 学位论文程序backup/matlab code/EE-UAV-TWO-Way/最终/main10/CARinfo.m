function [GVR,GVU,distanceVR,distanceVU,CARposition,L] = CARinfo(H,T,M,det,roadlength,roadwidth,UAVcvx)
Vposition = zeros(M,3,T);% ���峵���ĳ�ʼλ��
ll=zeros(1,M);
lright=1;%�������������� right
lleft=-1;%��������������
Vvelocity = 10;  % �����ٶ�Ϊ40

Vvelocity=zeros(1,M);  % ÿ�������ٶȶ���һ��
for i=1:M
    Vvelocity(i)=randi([10, 20]);
end

while 1
for m = 1:M
random_number = randi([0, 1]); % ����0��1֮������������0����-1��1����1��
if random_number == 0
    ll(m) = lright;   % �����ܵĳ�
    random_number = randi([0, 1]);
    if random_number == 0
        Vposition_x=0;   
    else
        Vposition_x=roadlength/2;
    end   
    random_numbery = (roadwidth/2) * rand;  % ����4��8֮��������
else
    ll(m) = lleft; %�����ܵĳ�
    random_number = randi([0, 1]);
    if random_number == 0
        Vposition_x=roadlength/2;
    else
        Vposition_x=roadlength/4;
    end 
    random_numbery = roadwidth/2+(roadwidth/2) * rand; % ����0��4֮��������
end
Vposition_x=randi([200,roadlength-200]);  %�����������㣬֮ǰ��δʹ��
Vposition_x=randi([200,roadlength+200]);  %�����������㣬֮ǰ��δʹ��
%      Vposition(m,:,1)=[Vposition_x; roadwidth/2+roadwidth/2*ll(m); 0];
     Vposition(m,:,1)=[Vposition_x; random_numbery+roadwidth/2*ll(m); 0];%�ڶ����ǳ���������λ��
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


Varin=load ('��ʮ����ǧ��.mat');%�����ж�Ҫһ���ã�����CARposition��Ĭ�ϵĲ��ܸĶ�
Vposition = Varin.CARposition;  %������֮ǰ������S

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

UAVpositions=zeros(M,3,T);
GVU=zeros(M,T);
% ����Ҫ��ӵ�һ��
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
    GVU(m,t) =l/(distanceVU(m,t)^2);%·����ĵĵ���
    end
end

CARposition=Vposition;
L=ll;


