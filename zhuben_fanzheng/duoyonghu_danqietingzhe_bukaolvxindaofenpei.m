lou=1e-7;
sigma=1e-14;
fai=3;
y=1e+8;
VH_max=10;
VV_max=5;
t=1;
k=0.557;
PG_avg=2;
PG_max=4*PG_avg;
PJ_avg=0.01;
PJ_max=4*PJ_avg;
HG=[0;0];
HJ_I=[-300;-100];
HJ_F=[300;-100];
VJ_I=70;
VJ_F=70;
HU=[-200,200;200,200];
HE=[0;-200];
QE=30;
H_NFZ1=[-150;-150];
Q_NFZ1=10;
H_NFZ2=[150;-150];
Q_NFZ2=10;
cvx3=-10;
cvx2=-10;
BCD=-10;
%初始化无人机轨迹
for i=1:120
    a1=[-300+5*i;-100];
    HJ_0(:,i)=a1;
end
VJ_0=80*ones(1,120);
%初始化无人机功率
PJ_0=0.006*ones(1,120);
%初始化基站功率
PG_0=0.5*ones(2,120);

M_E0=106400*ones(1,120);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% for BCD_N=1:100
    %优化无人机轨迹
for m=1:2
    for n=1:120
         c_m(m,n)=exp(-k)*PG_0(m,n)*y/((norm(HG-HU(:,m)))^(fai));
    end
end
for m=1:2
    for n=1:120
        e_m(m,n)=PG_0(m,n)*y/((norm(HG-HE)-QE)^(fai));
    end
end

for sub3=1:100
cvx_begin
variable eta
variable HJ(2,120)
variable VJ(1,120)
variable D(2,120)
variable M_Um(2,120)
variable M_E(1,120)
variable lamda(1,120)
variable F_E(1,120)
expression a1(2,120)

for m=1:2
    for n=1:120
    a1(m,n)=log(1+c_m(m,n)-c_m(m,n)*PJ_0(n)*y*inv_pos(y*PJ_0(n)+M_Um(m,n)))/log(2)-D(n);
    end
end
maximize(eta)  
subject to
sum(a1(1,:))*(1/120)-eta>=0;
sum(a1(2,:))*(1/120)-eta>=0;
for m=1:2
    for n=1:120
        VJ_0(n)^(2)+(VJ(n)-VJ_0(n))*2*VJ_0(n)+(norm(HJ_0(:,n)-HU(:,m)))^(2)+2*(HJ(:,n)-HJ_0(:,n))'*(HJ_0(:,n)-HU(:,m))-M_Um(m,n)>=0;
    end
end
for m=1:2
    for n=1:120
    log(1+(e_m(m,n)/(1+PJ_0(n)*y/M_E0(n))))/log(2)+(M_E(n)-M_E0(n))*y*PJ_0(n)*e_m(m,n)/(log(2)*(PJ_0(n)*y+M_E0(n)+e_m(m,n)*M_E0(n))*(PJ_0(n)*y+M_E0(n)))-D(m,n)<=0;
    end
end
for n=1:120
    [lamda(n)-1,0,HJ(1,n)-HE(1,1);0,lamda(n)-1,HJ(2,n)-HE(2,1);HJ(1,n)-HE(1,1),HJ(2,n)-HE(2,1),(-1)*lamda(n)*QE^(2)+F_E(n)]==semidefinite(3);
end
for n=1:120
    lamda(n)>=0;
end
for n=1:120
     -(square_pos(norm(HJ(:,n)-HE))+VJ(n)^(2)-M_E(n))-F_E(n)>=0;
end
for n=1:119
    norm(HJ(:,n+1)-HJ(:,n))-VH_max*t<=0;
end
for n=1:119
    norm(VJ(n+1)-VJ(n))-VV_max*t<=0;
end
norm(HJ(:,1)-HJ_I)-VH_max*t<=0;
norm(HJ(:,120)-HJ_F)-VH_max*t<=0;
VJ(1)^(2)+VJ_I^(2)-2*VJ(1)*VJ_I<=VV_max^(2)*t;
VJ(120)^(2)+VJ_F^(2)-2*VJ(120)*VJ_F<=VV_max^(2)*t;
20*ones(1,120)<=VJ<=120*ones(1,120);
cvx_end
HJ_0=HJ;
VJ_0=VJ;
M_E0=M_E;
sub3_cishu=sub3
if ((cvx_optval-cvx3)<=0.01)
    break;
else
    cvx3=cvx_optval;
      end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:2
y_m(m)=exp(-k)*y/(norm(HG-HU(:,m)))^(fai);
end
y_e=y/(norm(HG-HE)-QE)^(fai);
for m=1:2
    for n=1:120
        h_Jm(m,n)=y/(VJ_0(n)^(2)+(norm(HJ_0(:,n)-HU(:,m)))^(2));
    end
end
for n=1:120
    h_Je(n)=y/(VJ_0(n)^(2)+(norm(HJ_0(:,n)-HE)+QE)^(2));
end
for sub2=1:100
cvx_begin
variable eta
variable PG(2,120)
variable PJ(1,120)
variable B(2,120)
expression a2(2,120)

for m=1:2
    for n=1:120
        a2(m,n)=log(1+PJ(n)*h_Jm(m,n)+y_m(m)*PG(m,n))/log(2)-log(1+PJ_0(n)*h_Jm(m,n))/log(2)-(PJ(n)-PJ_0(n))*h_Jm(m,n)/(log(2)*(1+PJ_0(n)*h_Jm(m,n)))-B(m,n);
    end
end
maximize(eta)  
subject to
sum(a2(1,:))*(1/120)-eta>=0;
sum(a2(2,:))*(1/120)-eta>=0;
for m=1:2
    for n=1:120
        log(1+PJ_0(n)* h_Je(n)+PG_0(m,n)*y_e)/log(2)+y_e*(PG(m,n)-PG_0(m,n))/(log(2)*(1+PJ_0(n)* h_Je(n)+PG_0(m,n)*y_e))+h_Je(n)*(PJ(n)-PJ_0(n))/(log(2)*(1+PJ_0(n)* h_Je(n)+PG_0(m,n)*y_e))-log(1+PJ(n)*h_Je(n))/log(2)-B(m,n)<=0;
    end
end
PG<=PG_max*ones(2,120);
PG>=zeros(2,120);
sum(PG(:))*(1/120)-PG_avg<=0;
PJ<=PJ_max*ones(1,120);
PJ>=zeros(1,120);
sum(PJ)*(1/120)-PJ_avg<=0;
cvx_end
PJ_0=PJ;
PG_0=PG;
sub2_cishu=sub2
if ((cvx_optval-cvx2)<=0.01)
    break;
else
    cvx2=cvx_optval;
end
end
% BCD_cishu=BCD_N %(鏄剧ず鎬婚棶棰樿凯浠ｆ鏁?
% if(eta-BCD<=0.01)
%  break;
% else
%     BCD=eta;
% end
% end


HJ0_y0=-100*ones(1,120);
HJ0_x0=[];
for i = 1:120
    a = -300+5*i;
     HJ0_x0=[HJ0_x0,a];
end
HJ0=[HJ0_x0;HJ0_y0];
figure(1)
plot(HJ(1,:),HJ(2,:),'-kv', ...%绘制优化后的轨迹（蓝色下三角，黑线）
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','g', ...
'MarkerFaceColor','b')
hold on
plot(HJ0(1,:),HJ0(2,:),'-kv', ...%绘制初始轨迹（黄色下三角，黑线）
'LineWidth',1, ...
'MarkerSize',5, ...
'MarkerEdgeColor','g', ...
'MarkerFaceColor','y')
hold on
HU1=[HU(1,1);HU(2,1)];
plot(HU1(1,1),HU1(2,1),'^', ...%用户一（红色上三角）
'MarkerSize',7, ...
'MarkerEdgeColor','k', ...
'MarkerFaceColor','r')
hold on
HU2=[HU(1,2);HU(2,2)];
plot(HU2(1,1),HU2(2,1),'^', ...%用户二（红色上三角）
'MarkerSize',7, ...
'MarkerEdgeColor','k', ...
'MarkerFaceColor','r')
hold on
plot(HE(1,1),HE(2,1),'v', ...%窃听者一（黑色下三角）
'MarkerSize',7, ...
'MarkerEdgeColor','k', ...
'MarkerFaceColor','k')

plot(HG(1,1),HG(2,1),'h', ...%基站（绿色六角星）
'MarkerSize',7, ...
'MarkerEdgeColor','k', ...
'MarkerFaceColor','g')
title('trajectory')
xlabel('m')
ylabel('m')
legend('HJ0','HJ','HE','HU1','HU2','HG')
% text(10,0,'user1')
% text(10,0,'user2')
% text(-170,150,'EVE1 fangcha=5(m)')
% text(150,-150,'EVE2 fangcha=30(m)')

P_UAV=[];
for i=1:120
    a=PJ(i);
    P_UAV(:,i)=[i;a];
end
figure(2)
plot(P_UAV(1,:),P_UAV(2,:),'.b-')
xlabel('Time(s)')
ylabel('Power(w)')
legend('P_UAV')
PG1=[];
for i=1:120
    a=PG(1,i);
    PG1(:,i)=[i;a];
end
figure(3)
plot(PG1(1,:),PG1(2,:),'.b-')
hold on

PG2=[];
for i=1:120
    a=PG(2,i);
    PG2(:,i)=[i;a];
end
plot(PG2(1,:),PG2(2,:),'.r-')
xlabel('Time(s)')
ylabel('Power(w)')
legend('PG1','PG2')
% for i=1:BCD_N
%     CNM(i)=i;
% end
% OPT=[CNM;MMP];
% figure(4)
% plot(OPT(1,:),OPT(2,:),'.b-')
% xlabel('Number of iterations')
% ylabel('Average secrecy rate (bps/Hz)')
% 
VC=[];
for i=1:120
    a=VJ(i);
    VC(:,i)=[i;a];
end
figure(5)
plot(VC(1,:),VC(2,:),'.b-')
xlabel('Time(s)')
ylabel('Height(m)')