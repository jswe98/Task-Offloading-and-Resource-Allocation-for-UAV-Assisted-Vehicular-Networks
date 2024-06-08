clc
close all
clear
% 总初始化
T=70; M=8;H=100; det=1;roadwidth=30;roadlength=800; %共有70个时隙 %solt length   TTT=T*det;%time span
VUAV_max=15;%无人机速度约束
W=10^(6);
theta=0.2;%乘在E_VR之前，无人机的会比较大，这是一个系数
Delta=1e-6;%背景噪声是
Pm=0.3*ones(T,M);
R=zeros(T,M);
[UAVposition] = InitializationUAVposition(T,roadlength,roadwidth);
UAVcvx=UAVposition;
c2_l=zeros(T,M);
c3_l=zeros(T,M);
c1=zeros(T,M);
c2=zeros(T,M);
omega_l=zeros(T,M);
rho_l=zeros(T,M);
distance_l=zeros(T,M);
q_l_SCA=UAVposition;
for bcd=1:2 %BCD开始  可以加while形成循环
UAVcvx=q_l_SCA;%结果返回了
[GVR,GVU,distanceVR,distanceVU,CARposition,L]=CARinfo(H,T,M,det,roadlength,roadwidth,UAVcvx);%使用了汽车信息函数
[X,Y] = Xdecision(Pm,GVU,GVR,T,M);%使用了决策函数

for t = 1:T %Qcvx里所必需的定量值
    for m = 1:M
        c2_l(t,m) = Pm(t,m)*GVR(m,t)/Delta;
        c3_l(t,m)=Pm(t,m)/Delta;
        omega_l(t,m)=(-1*c3_l(t,m)*log2(exp(1)))/(((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,2,t))^(2)+H^(2))*...
         ((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,1,t))^(2)+H^(2)+c3_l(t,m)));
         rho_l(t,m)=log2(1+c3_l(t,m)/((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,2,t))^(2)+H^(2)));
         distance_l(t,m)=((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,2,t))^(2));
        %第l次迭代，第t个时隙UAV与第m个浮标水平距离平方
    end
 end
distance_l_SCA=distance_l;
L_VR=sum(sum(W*X'.*log2(1+Pm.*c2_l)));
E_VR=sum(sum(Pm.*X'));%VRde能量消耗
E_VU=sum(sum(Pm.*Y'));%VUde能量消耗

epsilon2=0.1;%这个是丁收敛的阈值
flag=0; %开始循环的标志
chi1=0;%这个是除法之后的系数
ding_u=1;%看看循环了多少次 %求解Pm的dinkelbach轮次 ding-p表示求的是Pm
while (flag==0)                  %Dinkelbach
        flag=1;                %结束循环的标志
        chi1_rec(ding_u)=chi1;%chi1的record变量
[q,distance2] =Qcvx(theta,E_VU,E_VR,W,Y,L_VR,rho_l,omega_l,distance_l_SCA,CARposition,VUAV_max,det,T,M);%使用了Qcvx函数
q_l_SCA=q;%q是cvx里面的 variable   %记录本轮优化结果，供下轮使用
            for t=1:T
                for m=1:M
                    distance2(t,m)=square_pos(norm([q(t,1) q(t,2)]-[CARposition(m,1,t) CARposition(m,2,t)]));
                end
            end
            L_A=sum(sum(W*Y'.*(rho_l+omega_l.*(distance2-distance_l_SCA))));%凸+凹？
            L_total=L_A+L_VR;
        F_chi1=L_total -chi1*1*(theta*sum(sum(Pm.*X'))+(1-theta)*sum(sum(Pm.*Y'))) ;         % 这个是减法的
        chi1=L_total/(1*(theta*sum(sum(Pm.*X'))+(1-theta)*sum(sum(Pm.*Y'))));
        ding_u=ding_u+1;                %循环次数每次自己增加一
         if (F_chi1>epsilon2) %判断标志
             flag=0;
         end
end

chi2=0;%这个是除法之后的系数
ding_p=1;%看看循环了多少次 %求解Pm的dinkelbach轮次 ding-p表示求的是Pm
flag2=0; 
while (flag2==0)                  %Dinkelbach
        flag2=1;                %结束循环的标志
        chi2_rec(ding_p)=chi2;%chi1的record变量
        [P_m] = Pcvx(W,X,Y,Delta,GVR,GVU,theta,T,M,chi2) ;
        Pm=P_m;
      GVRR=GVR';
      GVUU=GVU';
    for t = 1:T %给出信道增益h(t,m)和c2(t,m)
        for m = 1:M
            c1(t,m) = GVRR(t,m)/Delta;
            c2(t,m) =  GVUU(t,m)/Delta;%原
        end
    end  %这里是为了求解轨迹规划后新的信道状态信息拿来优化Pm   1,1+P_b.*c1
        F_chi2=sum(sum(W*X'.*rel_entr(ones(T,M),1+P_m.*c1)))+sum(sum(W*Y'.*rel_entr(ones(T,M),1+P_m.*c2)))...
            -chi2*1*(theta*sum(sum(P_m.*X'))+(1-theta)*sum(sum(P_m.*Y'))) ;         % 这个是减法的
        chi2=(sum(sum(W*X'.*rel_entr(ones(T,M),1+P_m.*c1)))+sum(sum(W*Y'.*rel_entr(ones(T,M),1+P_m.*c2))))...
            /(1*(theta*sum(sum(P_m.*X'))+(1-theta)*sum(sum(P_m.*Y'))));
        ding_p=ding_p+1;                %循环次数每次自己增加一
         if (F_chi2>epsilon2) %判断标志
             flag2=0;
         end
         
end

%可以加入BCD何时结束，判定方法为abs()>0.1
end



figure;
hold on;
UAVX=zeros(T,1);UAVXX=UAVX;UAVY=UAVXX;UAVYY=UAVXX;
car_positions=zeros(M,T);
car_x=zeros(M,T);car_y=car_x;
legend_entries=cell(1,M+3);
for t=1:T
    UAVX(t)=UAVposition(t,1);
    UAVY(t)=UAVposition(t,2);
    UAVXX(t)=q(t,1);
    UAVYY(t)=q(t,2);
    for m = 1:M
    car_positions(m,1:2,t)=CARposition(m,1:2,t);
    car_x(m,t)=car_positions(m,1,t);
    car_y(m,t)=car_positions(m,2,t);
        if L(m)==-1;
           legend_entries{m} = ['leftVehicle ' num2str(m)]; 
        else
           legend_entries{m} = ['rightVehicle ' num2str(m)]; 
        end
    end
end
for m = 1:M
    if L(m)==-1;
    plot(car_x(m, :), car_y(m, :), 'LineWidth', 1,'Marker', '<','MarkerFaceColor', 'auto');  % 绘制每辆车的轨迹
    else
    plot(car_x(m, :), car_y(m, :), 'LineWidth', 1,'Marker', '>','MarkerFaceColor', 'auto');  % 绘制每辆车的轨迹  
    end
end
hold off;
% 绘制无人机未优化的初始轨迹
hold on;
plot(UAVX, UAVY, 'color', 'black','Marker', '+', 'LineWidth', 1.5)
hold off;
xlabel('马路长度');
ylabel('马路宽度');
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
y_center = mean(ylim);
% 绘制无人机轨迹优化
hold on;
plot(UAVXX, UAVYY, 'color', 'blue','Marker', '+', 'LineWidth', 1.5)
hold off;
% 添加归一化的水平箭头
x_start = 0.5;
y_position = 0.6;    % 箭头的垂直位置
x_length = 0.19;      % 箭头的水平长度
y_offset = 0.1;      % 箭头的垂直位移
hold on;
annotation('arrow', [x_start, x_start + x_length], [y_position + y_offset, y_position + y_offset], 'Units', 'normalized','color', 'black');
hold off;
 % 添加黄色双实线
hold on;
plot(xlim, [y_center+0.05, y_center+0.05],'-',  'color', 'yellow', 'LineWidth', 2);
plot(xlim, [y_center-0.05, y_center-0.05],'-', 'color', 'yellow', 'LineWidth', 2);
hold off;
legend_entries{M+1} = '无人机初始轨迹';
legend_entries{M+2} = 'Optimization UAV Trajectory';
legend_entries{M+3} = '黄色双实线';
% legend_entries{M+3} = '向右为正方向';
%legend_entries{M+4} = '黄色双实线';
legend(legend_entries, 'Location', 'best');
XX=1-X'


