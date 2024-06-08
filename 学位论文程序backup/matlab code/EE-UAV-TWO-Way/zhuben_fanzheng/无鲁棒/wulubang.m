T=70;
MMP=[];
fai=3;
y=1e+8;
VH_max=10;%无人机速度约束
% VV_max=5;
t=1;
k1=0.557;
PG_avg=2.7;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PG_max=2*PG_avg;
PJ_avg=0.0012;
PJ_max=4*PJ_avg;
HG=[0;0];
HJ_I=[-200;-100];%无人机水平起点、终点位置
HJ_F=[200;-100];
% VJ_I=70;%无人机垂直起点、终点位置
% VJ_F=70;
HU=[-150,200;200,200];%合法用户位置
HE=[-100,50;-350,-380];%窃听者估计位置
% QE=[10,30];%估计误差
HE_T=[-90,90;-360,-382];
H_NFZ=[-100,50;-200,-250];%飞行禁区中心位置
Q_NFZ=[45,30];%飞行禁区半径
BCD=-20;
for i=1:T %初始化无人机水平轨迹
    a=[-200+(400/T)*i;-100];
    HJ_0(:,i)=a;
end
% HJ_0=xlsread('HJ.xlsx');
VJ_0=70*ones(1,T);%初始化无人机垂直轨迹
PJ_0=0.0007*ones(1,T);%初始化无人机功率
PG_0=0.6*ones(2,T);%初始化基站功率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_E0=160000*ones(2,T);
% M_E0=xlsread('M_E.xlsx');
for BCD_N=1:3
    %优化无人机轨迹
for m=1:2
    for n=1:T
         c_m(m,n)=exp(-k1)*PG_0(m,n)*y/((norm(HG-HU(:,m)))^(fai));
    end
end
for m=1:2
    for k=1:2
        for n=1:T
             e_m_k(m,k,n)=PG_0(m,n)*y/((norm(HG-HE(:,k)))^(fai));
        end
    end
end
 
cvx1=-20;
for sub1=1:100
cvx_begin
variable eta
variable HJ(2,T)
% variable VJ(1,90)
variable D(2,T)
variable M_Um(2,T)
variable M_E(2,T)
% variable lamda(2,90)
% variable F_E(2,90)
expression a1(2,T)

for m=1:2
    for n=1:T
    a1(m,n)=log(1+c_m(m,n)-c_m(m,n)*PJ_0(n)*y*inv_pos(y*PJ_0(n)+M_Um(m,n)))/log(2)-D(m,n);
    end
end
maximize(eta)  
subject to
sum(a1(1,:))*(1/T)-eta>=0;
sum(a1(2,:))*(1/T)-eta>=0;
for m=1:2
    for n=1:T
        VJ_0(n)^(2)+(norm(HJ_0(:,n)-HU(:,m)))^(2)+2*(HJ(:,n)-HJ_0(:,n))'*(HJ_0(:,n)-HU(:,m))-M_Um(m,n)>=0; %VJ_0(n)^(2)+(VJ(n)-VJ_0(n))*2*VJ_0(n)+
    end
end
for m=1:2
    for k=1:2
        for n=1:T
            log(1+(e_m_k(m,k,n)/(1+PJ_0(n)*y/M_E0(k,n))))/log(2)+(M_E(k,n)-M_E0(k,n))*y*PJ_0(n)*e_m_k(m,k,n)/(log(2)*(PJ_0(n)*y+M_E0(k,n)+e_m_k(m,k,n)*M_E0(k,n))*(PJ_0(n)*y+M_E0(k,n)))-D(m,n)<=0;
        end
    end
end

% for k=1:2
%     for n=1:T
%         [lamda(k,n)-1,0,HJ(1,n)-HE(1,k);0,lamda(k,n)-1,HJ(2,n)-HE(2,k);HJ(1,n)-HE(1,k),HJ(2,n)-HE(2,k),(-1)*lamda(k,n)*QE(k)^(2)+F_E(k,n)]==semidefinite(3);%%不确定
%     end
% end

% for k=1:2
%     for n=1:T
%         lamda(k,n)>=0;
%     end
% end

for k=1:2
    for n=1:T
       (VJ_0(n)^2)+square_pos(norm(HJ(:,n)-HE(:,k)))-M_E(k,n)<=0;% -(square_pos(norm(HJ(:,n)-HE(:,k)))+VJ(n)^(2)-M_E(k,n))-F_E(k,n)>=0;
    end
end

for n=1:(T-1)
    norm(HJ(:,n+1)-HJ(:,n))-VH_max*t<=0;
end
% for n=1:(T-1)
%     norm(VJ(n+1)-VJ(n))-VV_max*t<=0;
% end
norm(HJ(:,1)-HJ_I)-VH_max*t<=0;
HJ(:,T)==HJ_F;
% norm(VJ(1)-VJ_I)<=VV_max*t;
% norm(VJ(T)-VJ_F)<=VV_max*t;
% 60*ones(1,T)<=VJ<=120*ones(1,T);%飞行高度约束
for d=1:2
    for n=1:T
        (norm(HJ_0(:,n)-H_NFZ(:,d)))^2+2*(HJ_0(:,n)-H_NFZ(:,d))'*(HJ(:,n)-HJ_0(:,n))>=Q_NFZ(d)^2;%飞行禁区约束
    end
end
cvx_end
sub1_cishu=sub1
if (abs((cvx_optval-cvx1)/cvx_optval)<=0.01)
  %************************************************************************  
%      xlswrite('HJ.xlsx',HJ)
%       xlswrite('M_E.xlsx',M_E)
%**************************************************************************8  
break;
else
    cvx1=cvx_optval;
      HJ_0=HJ;
%       VJ_0=VJ;
      M_E0=M_E;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:2
y_m(m)=exp(-k1)*y/(norm(HG-HU(:,m)))^(fai);
end
for k=1:2
    y_e(k)=y/(norm(HG-HE(:,k)))^(fai);
end
for m=1:2
    for n=1:T
        h_Jm(m,n)=y/(VJ_0(n)^(2)+(norm(HJ_0(:,n)-HU(:,m)))^(2));
    end
end
for k=1:2
    for n=1:T
        h_Je(k,n)=y/(VJ_0(n)^(2)+(norm(HJ_0(:,n)-HE(:,k)))^(2));
    end
end
cvx2=-20;
for sub2=1:100
cvx_begin
variable eta
variable PG(2,T)
variable PJ(1,T)
variable B(2,T)
expression a2(2,T)

for m=1:2
    for n=1:T
        a2(m,n)=log(1+PJ(n)*h_Jm(m,n)+y_m(m)*PG(m,n))/log(2)-log(1+PJ_0(n)*h_Jm(m,n))/log(2)-(PJ(n)-PJ_0(n))*h_Jm(m,n)/(log(2)*(1+PJ_0(n)*h_Jm(m,n)))-B(m,n);
    end
end
maximize(eta)  
subject to
sum(a2(1,:))*(1/T)-eta>=0;
sum(a2(2,:))*(1/T)-eta>=0;
for m=1:2
    for k=1:2
        for n=1:T
            log(1+PJ_0(n)* h_Je(k,n)+PG_0(m,n)*y_e(k))/log(2)+y_e(k)*(PG(m,n)-PG_0(m,n))/(log(2)*(1+PJ_0(n)* h_Je(k,n)+PG_0(m,n)*y_e(k)))+h_Je(k,n)*(PJ(n)-PJ_0(n))/(log(2)*(1+PJ_0(n)* h_Je(k,n)+PG_0(m,n)*y_e(k)))-log(1+PJ(n)*h_Je(k,n))/log(2)-B(m,n)<=0;
        end
    end
end

PG<=PG_max*ones(2,T);
PG>=zeros(2,T);
sum(PG(:))*(1/T)-PG_avg<=0;
PJ<=PJ_max*ones(1,T);
PJ>=zeros(1,T);
sum(PJ)*(1/T)-PJ_avg<=0;
cvx_end
sub2_cishu=sub2
if (abs((cvx_optval-cvx2)/cvx_optval)<=0.01)
%        *************************************************************
%     xlswrite('PG.xlsx',PG)
%      xlswrite('PJ.xlsx',PJ)
     %**************************************************************
    break;
else
    cvx2=cvx_optval;
    PJ_0=PJ;
    PG_0=PG;
end
end

BCD_cishu=BCD_N %(鏄剧ず鎬婚棶棰樿凯浠ｆ鏁?
if(abs((cvx2-BCD)/cvx2)<=0.01)
 break;
else
    BCD=cvx2;
    MMP(BCD_cishu)=BCD;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 R_e_m=[];
  R_e=[];
   R_m=[];
    R_sec=[];
for m=1:2
    for k=1:2
        for n=1:T
           R_e_m(m,k,n)=log2(1+PG(m,n)*y*(norm(HG-HE_T(:,k)))^(-fai)/(1+PJ(n)*y/(VJ_0(n)^(2)+(norm(HJ(:,n)-HE_T(:,k)))^(2)))); 
        end
    end
end

for m=1:2
        for n=1:T
            if R_e_m(m,1,n)>=R_e_m(m,2,n)
                R_e(m,n)=R_e_m(m,1,n);
            else
                 R_e(m,n)=R_e_m(m,2,n);
            end
        end
end
for m=1:2
    for n=1:T
        R_m(m,n)=log2(1+exp(-k1)*PG(m,n)*y*(norm(HG-HU(:,m)))^(-fai)/(1+PJ(n)*y/(VJ_0(n)^(2)+(norm(HJ(:,n)-HU(:,m)))^(2)))); 
    end
end
for m=1:2
    for n=1:T
        R_sec(m,n)=R_m(m,n)-R_e(m,n);
    end
end
R_avg=sum( R_sec(:))*(1/T)*0.5
% xlswrite('R_avg_wulubang_90.xlsx',R_avg)
% xlswrite('wulubang_T90.xlsx',HJ)
xlswrite('R_avg_wulubang_PG_2.7_0.0012.xlsx',R_avg)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% xlswrite('R_avg_wulubang_PJ_0.0027.xlsx',R_avg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HJ0_y0=-100*ones(1,T);
HJ0_x0=[];
for i = 1:T
    a = -200+(400/T)*i;
     HJ0_x0=[HJ0_x0,a];
end
HJ0=[HJ0_x0;HJ0_y0];
figure(1)
plot(HJ_0(1,:),HJ_0(2,:),'-kv', ...%绘制优化后的轨迹（蓝色下三角，黑线）
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
HE1=[HE(1,1);HE(2,1)];
plot(HE1(1,1),HE1(2,1),'v', ...%窃听者一（黑色下三角）
'MarkerSize',7, ...
'MarkerEdgeColor','k', ...
'MarkerFaceColor','k')
hold on
HE2=[HE(1,2);HE(2,2)];
plot(HE2(1,1),HE2(2,1),'v', ...%窃听者二（黑色下三角）
'MarkerSize',7, ...
'MarkerEdgeColor','k', ...
'MarkerFaceColor','k')
hold on
plot(HG(1,1),HG(2,1),'h', ...%基站（绿色六角星）
'MarkerSize',7, ...
'MarkerEdgeColor','k', ...
'MarkerFaceColor','g')
title('trajectory')
xlabel('m')
ylabel('m')
legend('HJ0','HJ_0','HE1','HE2','HU1','HU2','HG')
% text(10,0,'user1')
% text(10,0,'user2')
% text(-170,150,'EVE1 fangcha=5(m)')
% text(150,-150,'EVE2 fangcha=30(m)')

PU=[];
for i=1:T
    a=PJ_0(i);
    PU(:,i)=[i;a];
end
figure(2)
plot(PU(1,:),PU(2,:),'.b-')
xlabel('Time(s)')
ylabel('Power(w)')
legend('PU')
PG1=[];
for i=1:T
    a=PG_0(1,i);
    PG1(:,i)=[i;a];
end
figure(3)
plot(PG1(1,:),PG1(2,:),'.b-')
hold on

PG2=[];
for i=1:T
    a=PG_0(2,i);
    PG2(:,i)=[i;a];
end
plot(PG2(1,:),PG2(2,:),'.r-')
xlabel('Time(s)')
ylabel('Power(w)')
legend('PG1','PG2')
% CNM=[];
% for i=1:(BCD_N)
%     CNM(i)=i;
% end
% OPT=[CNM;MMP]
% figure(4)
% plot(OPT(1,:),OPT(2,:),'.b-')
% xlabel('Number of iterations')
% ylabel('Average secrecy rate (bps/Hz)')

% VC=[];
% for i=1:90
%     a=VJ_0(i);
%     VC(:,i)=[i;a];
% end
% figure(5)
% plot(VC(1,:),VC(2,:),'.b-')
% xlabel('Time(s)')
% ylabel('Height(m)')

FF=sum(PJ)
% R1=[];
% for n=1:90
%     a=(log(1+PJ(n)*h_Jm(1,n)+y_m(1)*PG(1,n))/log(2)-log(1+PJ_0(n)*h_Jm(1,n))/log(2)-(PJ(n)-PJ_0(n))*h_Jm(1,n)/(log(2)*(1+PJ_0(n)*h_Jm(1,n)))-B(1,n));
%     R1(:,n)=[n;a];
% end
% figure(6)
% plot(R1(1,:),R1(2,:),'.b-')
% hold on
% 
% R2=[];
% for n=1:90
%     a=log(1+PJ(n)*h_Jm(2,n)+y_m(2)*PG(2,n))/log(2)-log(1+PJ_0(n)*h_Jm(2,n))/log(2)-(PJ(n)-PJ_0(n))*h_Jm(2,n)/(log(2)*(1+PJ_0(n)*h_Jm(2,n)))-B(2,n);
%     R2(:,n)=[n;a];
% end
% plot(R2(1,:),R2(2,:),'.r-')
% xlabel('Time(s)')
% ylabel('Average secrecy rate (bps/Hz)')
% legend('R1','R2')