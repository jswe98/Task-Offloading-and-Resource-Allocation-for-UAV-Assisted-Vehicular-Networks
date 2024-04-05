clc
close all
clear
% 总初始化
T=70; M=4;H=100; det=1;roadwidth=30;roadlength=800; %共有70个时隙 %solt length   TTT=T*det;%time span
VUAV_max=20;%无人机速度约束
W=10^(6);
theta=0.4;%乘在VR之前，无人的会比较大
Delta=1e-6;%背景噪声是
Pm=0.03*ones(T,M);
R=zeros(T,M);
[UAVposition] = InitializationUAVposition(T,roadlength,roadwidth);
UAVcvx=UAVposition;
c2_l=zeros(T,M);
c3_l=zeros(T,M);
omega_l=zeros(T,M);
rho_l=zeros(T,M);
distance_l=zeros(T,M);
q_l_SCA=UAVposition;
for bcd=1:2
[GVR,GVU,distanceVR,distanceVU,CARposition,L]=CARinfo(H,T,M,det,roadlength,roadwidth,UAVcvx);
[X,Y] = Xdecision(Pm,GVU,GVR,T,M);
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
[q,distance2] =Qcvx(theta,E_VU,E_VR,W,Y,L_VR,rho_l,omega_l,distance_l_SCA,CARposition,VUAV_max,det,T,M);
q_l_SCA=q;
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



