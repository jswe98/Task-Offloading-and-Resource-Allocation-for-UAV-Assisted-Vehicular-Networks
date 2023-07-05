%确定gama(y),窃听者估计偏差e1,e2,平均功率P.avg%峰�?功率P.max,中断概率lou,与放缩后的中断概率lou_,无人机飞行高度H，合法用户的位置qu,
%窃听者位置的估计值，q1,q2;无人机的�?��飞行速度V_max;划分的时隙d;无人机起点与终点坐标QI,QF
%确定无人机初始轨迹Qb_0,功率的�?数T_0
%确定松弛变量Z_0；F_0；U_0的初始�?
P_avg=0.001;
P_max=4*P_avg;
y=1e+8;
e1=5;
e2=35;
lou=0.1;
lou_=1-(1-lou)^0.5;
H=100;
qu=[0;0];
q1=[-200;150] ;
q2=[200;-150];
V_max=10;
d=1;
QI=[0;240];
QF=[0;-240];
Qb_0=[0;-240];
T_0=[];
F_0=[];
U_0=[];
cvx=0;
OPT=[];
% 初始化无人机轨迹的可行点Qb_0，以6m/s飞行
  %初始化无人机可行点的横坐�?
  Qb_x0=zeros(1,80);
  %初始化无人机可行点的纵坐�?
  Qb_y0 =[];
for i = 1:80
    a = 240-6*i;
     Qb_y0=[Qb_y0,a];
end
Qb_0=[Qb_x0;Qb_y0];
%Z的初始可行点Z_0,与不同窃听�?下的可行点Z_1,Z_2
Z_1=[];
Z_2=[];
%确定初始可行点Z_0
for i=Qb_0 %对矩阵的列向量进行遍�?
    a=norm([e1^2;0;0;e1^2;sqrt(2)*e1*(q1-i)],2);
    Z_1=[Z_1,a];
end
for i=Qb_0
    a=norm([e2^2;0;0;e2^2;sqrt(2)*e2*(q2-i)],2);
    Z_2=[Z_2,a];
end
%初始化功率�?数（T_0)，松弛变量Fai（F)的初始点F_0,以及不同窃听者下的初始点F_1,F_2
F_1=[];
F_2=[];
% 初始化功率�?数（T_0)的可行点为平均功率的�?��
T_0=1/(0.5*P_avg)*ones(1,80);
%确定初始点F_0
for i= 1:80
    a=y/(T_0(i)*(2*(e1^2)-Z_1(i)*sqrt(-2*log(lou_))+H^2+(norm(Qb_0(:,i)-q1,2))^2));
    F_1=[F_1,a];
end
for i= 1:80
   a=y/(T_0(i)*(2*(e2^2)-Z_2(i)*sqrt(-2*log(lou_))+H^2+(norm(Qb_0(:,i)-q2,2))^2));
    F_2=[F_2,a];
end
for i= 1:80
    if(F_1(i)>F_2(i))
        F_0=[F_0,F_1(i)];
    else F_0=[F_0,F_2(i)];
    end
end
%初始可行点U_0
for i=1:80
    a=y/(T_0(i)*(H^2+(norm(Qb_0(:,i)-qu,2))^2));
    U_0=[U_0,a];
end
%调用CVX，求解凸化后的问�?

for n=1:100
    cvx_begin
variable Qb(2,80)
variable T(1,80)
variable F(1,80)
variable U(1,80)
variable Z1(1,80)
variable Z2(1,80)
expression obj(1,80)
for i=1:80
    a=log(1+U(i))/log(2)-(log(1+F_0(i))/log(2)+(F(i)-F_0(i))*1/(log(2)*(1+F_0(i))));
    obj(i)=1/80*a;
end
expression st_1_1(1,80)
for i=1:80
    a=2*e1^2-sqrt(-2*log(lou_))*Z1(i)-Qb(:,i)'*Qb(:,i)+2*Qb_0(:,i)'*Qb(:,i)+q1'*q1-2*Qb(:,i)'*q1+H^2-y*prod_inv([F(i) T(i)]);
    st_1_1(i)=a;
end
expression st_1_2(1,80)
for i=1:80
    a=2*e2^2-sqrt(-2*log(lou_))*Z2(i)-Qb(:,i)'*Qb(:,i)+2*Qb_0(:,i)'*Qb(:,i)+q2'*q2-2*Qb(:,i)'*q2+H^2-y*prod_inv([F(i),T(i)]);
    st_1_2(i)=a;
end
expression st_2_1(1,80)
for i=1:80
    st_2_1(i)=norm([e1^2;0;0;e1^2;2^0.5*e1*(q1-Qb(:,i))])-Z1(i);
end
expression st_2_2(1,80)
for i=1:80
    st_2_2(i)=norm([e2^2;0;0;e2^2;2^0.5*e2*(q2-Qb(:,i))])-Z2(i);
end
expression st_3(1,80)
for i=1:80
    st_3(i)=H^2+Qb(:,i)'*Qb(:,i)+qu'*qu-2*Qb(:,i)'*qu-y/(U_0(i)*T_0(i))+((U(i)-U_0(i))*y)/(T_0(i)*U_0(i)^2)+((T(i)-T_0(i))*y)/(U_0(i)*T_0(i)^2);
end
expression st_4(1,79)
for i=1:79
    st_4(i)=norm(Qb(:,i+1)-Qb(:,i))-V_max*d;
end
expression st_5_1(1,80)
for i=1:80
    st_5_1(i)=1/80*inv_pos(T(i));
end
maximize(sum(obj))  
subject to 
st_1_1>=zeros(1,80);
st_1_2>=zeros(1,80);
st_2_1<=zeros(1,80);
st_2_2<=zeros(1,80);
st_3<=zeros(1,80);
st_4<=zeros(1,79);
F>zeros(1,80);
U>=zeros(1,80);
norm(Qb(:,1)-QI)-V_max*d<=0;
norm(Qb(:,80)-QF)-V_max*d<=0;
sum(st_5_1)<=P_avg;
T>=(1/P_max)*ones(1,80);
cvx_end
Qb_0=Qb;
Z_1=Z1;
Z_2=Z2;
F_0=F;
T_0=T;
U_0=U;
OPT(:,n)=[n;cvx_optval]
if ((cvx_optval-cvx)<=0.001)
    break;
else
    cvx=cvx_optval;
end
end
Qb0_x0=zeros(1,80);
 Qb0_y0 =[];
for i = 1:80
    a = 240-6*i;
     Qb0_y0=[Qb0_y0,a];
end
Qb0=[Qb0_x0;Qb0_y0];
figure(1)
plot(Qb(1,:),Qb(2,:),'-kv', ...
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','g', ...
'MarkerFaceColor','b')
hold on
plot(Qb0(1,:),Qb0(2,:),'-kv', ...
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','g', ...
'MarkerFaceColor','y')
hold on
plot(q1(1,1),q1(2,1),'^', ...
'MarkerSize',7, ...
'MarkerEdgeColor','k', ...
'MarkerFaceColor','r')
hold on
plot(q2(1,1),q2(2,1),'^', ...
'MarkerSize',7, ...
'MarkerEdgeColor','k', ...
'MarkerFaceColor','r')
hold on
plot(qu(1,1),qu(2,1),'h', ...
'MarkerSize',7, ...
'MarkerEdgeColor','k', ...
'MarkerFaceColor','g')
title('trajectory')
xlabel('m')
ylabel('m')
text(10,0,'user')
text(-170,150,'EVE1 fangcha=5(m)')
text(150,-150,'EVE2 fangcha=30(m)')
P=[];
for i=1:80
    a=1/T(i);
    P(:,i)=[i;a];
end
figure(2)
plot(P(1,:),P(2,:),'.b-')
xlabel('Time(s)')
ylabel('Power(w)')
figure(3)
plot(OPT(1,:),OPT(2,:),'.b-')
xlabel('Number of iterations')
ylabel('Average secrecy rate (bps/Hz)')
