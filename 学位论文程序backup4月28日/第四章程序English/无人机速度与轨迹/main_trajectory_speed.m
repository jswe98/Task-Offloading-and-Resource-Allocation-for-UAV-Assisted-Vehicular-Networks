clc
close all 
clear
% 总初始化
T=90; M=6;H=100; det=1;roadwidth=30;roadlength=3000;VUAV_max=45;W=10^(3);Delta=1e-6;%背景噪声是%共有70个时隙 %solt length   TTT=T*det;%time span
fix_x = 2500;  % x坐标  
fix_y = 30; % y坐标  
cvx_quiet(true);
valuefix = zeros(T, 2);  % 创建一个与时间点数量相同的2列向量   % 模拟无人机在每个时间点的位置  
for i = 1:T  
    valuefix(i,1) = fix_x;  % 设置x坐标  
    valuefix(i,2) = fix_y;  % 设置y坐标   
end  
thetatheta=1; 
programmatic=1;
theta=0.4;%乘在E_VR之前，无人机的会比较大，这是一个系数
% for programmatic = 1:3  ;thetatheta=1; 
% for thetatheta=1:3 ;programmatic=1;
% theta=0.1+0.2*thetatheta;
[c2_l,c3_l,c1,c2,omega_l,rho_l,distance_l,R] =deal(zeros(T,M));
[UAVposition] = InitializationUAVposition(T,roadlength,roadwidth);
UAVcvx=UAVposition;
Pm=0.3*ones(T,M);
q_l_SCA=UAVposition;
flagBCD=0;%BCD开始  可以加while形成循环
BCDl=1;%BCD迭代轮数
EE_total_recod(1)=0;
      while (flagBCD==0)  %for bcd=1:5 %BCD开始  可以加while形成循环
             flagBCD = 1;
            if (BCDl>1)
               q_l_SCA=storage_BCD(BCDl-1).Q;
               Pm=storage_BCD(BCDl-1).Pb;
            end
UAVcvx=q_l_SCA;%结果返回了
[GVR,GVU,distanceVR,distanceVU,CARposition,L]=CARinfo(H,T,M,det,roadlength,roadwidth,UAVcvx,programmatic,thetatheta,BCDl);%使用了汽车信息函数
[X,Y] = Xdecision(Pm,GVU,GVR,T,M);%使用了决策函数
            for t = 1:T %Qcvx里所必需的定量值
                for m = 1:M
                    c2_l(t,m) = Pm(t,m)*GVR(m,t)/Delta;%车与基站
                    c3_l(t,m)=Pm(t,m)/Delta;
                    omega_l(t,m)=(-1*c3_l(t,m)*log2(exp(1)))/(((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,2,t))^(2)+H^(2))*...
                    ((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,1,t))^(2)+H^(2)+c3_l(t,m)));
                    rho_l(t,m)=log2(1+c3_l(t,m)/((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,2,t))^(2)+H^(2)));
                    distance_l(t,m)=((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,2,t))^(2));
                    %第l次迭代，第t个时隙UAV与第m个车水平距离平方
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
                  [q,distance2] =Qcvx(theta,E_VU,E_VR,W,Y,L_VR,rho_l,omega_l,distance_l_SCA,CARposition,VUAV_max,det,T,M,roadlength,roadwidth);%使用了Qcvx函数
                  q_l_SCA=q;%q是cvx里面的 variable   %记录本轮优化结果，供下轮使用
                  switch programmatic  
                         case 1  
                              q =q;  %优化轨迹
                         case 2  
                              q =q;  %固定轨迹
                         case 3  
                              qq = q;  %固定点
                         otherwise  
                              disp('Invalid value of programmatic');  
                  end 
                  for t=1:T
                      for m=1:M
                          distance2(t,m)=square_pos(norm([q(t,1) q(t,2)]-[CARposition(m,1,t) CARposition(m,2,t)]));
                      end
                  end
L_A=sum(sum(W*Y'.*(rho_l+omega_l.*(distance2-distance_l_SCA))));%凸+凹？
L_total=L_A+L_VR;
F_chi1=L_total -chi1*1*((1-theta)*sum(sum(Pm.*X'))+theta*sum(sum(Pm.*Y'))) ;         % 这个是减法的
chi1=L_total/(1*((1-theta)*sum(sum(Pm.*X'))+theta*sum(sum(Pm.*Y'))));       % 这个是除法的
ding_u=ding_u+1;  %循环次数每次自己增加一
                  if (F_chi1>epsilon2) %判断标志
                     flag=0;
                  end
             end
storage_BCD(BCDl).Q = UAVcvx;
q_l_SCA = storage_BCD(BCDl).Q;            
chi2=0;%这个是除法之后的系数
ding_p=1;%看看循环了多少次 %求解Pm的dinkelbach轮次 ding-p表示求的是Pm
flag2=0; 
             while (flag2==0)                  %Dinkelbach
                   flag2=1;                %结束循环的标志
                   chi2_rec(ding_p)=chi2;%chi1的record变量
                   [P_m] = Pcvx(W,X,Y,Delta,GVR,GVU,theta,T,M,chi2) ;
                   Pm=P_m;%将功率更新为最优值
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
storage_BCD(BCDl).Pb = Pm;        
L_VR_real = sum(W*X'.*log2(1+P_m.*c1));%   L_VR_real = sum(W*X'.*log2(1+P_m.*c1));   
L_VU_real = sum(W*Y'.*log2(1+P_m.*c2));
E_VR_real = sum(P_m.*X');
E_VU_real = sum(P_m.*Y');   
EE_total_recod(BCDl+1) = (sum(L_VR_real)+sum(L_VU_real))/...
(theta*sum(E_VR_real)+(1-theta)*sum(E_VU_real));
             if (abs(EE_total_recod(BCDl+1)-EE_total_recod(BCDl))/EE_total_recod(BCDl+1)>0.01) %可以加入BCD何时结束，判定方法为abs()>0.1
                 flagBCD=0;
                 BCDl=BCDl+1; 
             end
      end
bound1=1;
bound2=T;
[Ybound,Xbound,GVRbound,GVUbound,sinrr,sinru,E_VRbound,E_VUbound,Total_T] = deal(zeros(M, bound2-bound1+1));
[Pbound,qbound] = deal(zeros(bound2-bound1+1,M));
      for j=bound1:bound2
          for i=1:M
              Ybound(i,j-bound1+1)=Y(i,j-bound1+1);
              Xbound(i,j-bound1+1)=X(i,j-bound1+1);
              Pbound(j-bound1+1,i)=P_m(j-bound1+1,i);
              GVRbound(i,j-bound1+1)=GVR(i,j-bound1+1);
              GVUbound(i,j-bound1+1)=GVU(i,j-bound1+1);
              qbound(j-bound1+1,i)=distance2(j-bound1+1,i);
              sinrr(i,j-bound1+1)= 1/log(2)*(log(Pbound(j-bound1+1,i)*GVRbound(i,j-bound1+1)+Delta)-log(Delta)) ; 
              sinru(i,j-bound1+1)= 1/log(2)*(log(Pbound(j-bound1+1,i)*qbound(j-bound1+1,i)*GVUbound(i,j-bound1+1)+Delta)-log(Delta)) ; 
              E_VRbound(i,j-bound1+1)=Pbound(j-bound1+1,i).*Xbound(i,j-bound1+1);%VRde能量消耗
              E_VUbound(i,j-bound1+1)=Pbound(j-bound1+1,i).*Ybound(i,j-bound1+1);%VRde能量消耗
              Total_T(i,j-bound1+1)=(sinrr(i,j-bound1+1)+sinru(i,j-bound1+1))/((theta)*(E_VRbound(i,j-bound1+1)+(1-theta)*E_VUbound(i,j-bound1+1)));
              Total_T(j-bound1+1)=(sum(sinrr(:,j-bound1+1))+sum(sinru(:,j-bound1+1)))/((theta)*E_VRbound(i,j-bound1+1)+(1-theta)*E_VUbound(i,j-bound1+1));
          end
      end
% Total_theta(thetatheta,:)=0.01*sum(Total_T,1)
Total_theta(programmatic,:)=0.01*sum(Total_T,1);
% end %这个end表示最外面一个for循环的对应
figure;
zihao=24;
hold on;
grid on
box on
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
        if L(m)==-1
           legend_entries{m} = ['LeftVehicle ' num2str(m)]; 
           % legend_entries{m} = ['\fontname{Times New Roman}leftvehicle', '\fontname{Times New Roman}{' num2str(m) '}']; 
        else
           % legend_entries{m} = ['\fontname{Times New Roman}rightvehicle', '\fontname{Times New Roman}{' num2str(m) '}'];   %num2str(m)
           legend_entries{m} = ['RightVehicle ' num2str(m)]; 
        end
    end
end
set(gca,'FontName','Times New Roman')
for m = 1:M
    if L(m)==-1;
    plot(car_x(m, :), car_y(m, :), 'LineWidth', 1,'Marker', '<','MarkerFaceColor', 'auto');  % 绘制每辆车的轨迹
    else
    plot(car_x(m, :), car_y(m, :), 'LineWidth', 1,'Marker', '>','MarkerFaceColor', 'auto');  % 绘制每辆车的轨迹  
    end
end
hold off;
hold on;  % 绘制无人机未优化的初始轨迹
plot(UAVX, UAVY, 'color', 'black','Marker', '+', 'LineWidth', 1.5)
hold off;
xlabel('\fontname{Times New Roman}road length', 'FontSize', zihao);
ylabel('\fontname{Times New Roman}road widch', 'FontSize', zihao);
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
y_center = mean(ylim);
hold on; % 绘制无人机轨迹优化
plot(UAVXX, UAVYY, 'color', 'blue','Marker', '+', 'LineWidth', 1.5)
hold off;
x_start = 0.2; % 添加归一化的水平箭头
y_position = 0.6;    % 箭头的垂直位置
x_length = 0.19;      % 箭头的水平长度
y_offset = 0.1;      % 箭头的垂直位移
hold on;
annotation('arrow', [x_start, x_start + x_length], [y_position + y_offset, y_position + y_offset], 'Units', 'normalized','color', 'black');
hold off;
hold on;  % 添加黄色双实线
plot(xlim, [y_center+0.09, y_center+0.09],'-',  'color', 'yellow', 'LineWidth', 2);
plot(xlim, [y_center-0.09, y_center-0.09],'-', 'color', 'yellow', 'LineWidth', 2);
hold off;
legend_entries{M+1} = '\fontname{Times New Roman}Fix UAV trajectory';
legend_entries{M+2} = '\fontname{Times New Roman}Optimal UAV trajectory';
legend_entries{M+3} = '\fontname{Times New Roman}Double solid line';
legend_handle = legend(legend_entries, 'FontSize', zihao); 
legend(legend_entries, 'Location', 'best', 'FontSize', zihao-1);
set(gcf, 'Position', [100, 100, 2000, 1000]);
set(gca,'YTickLabel',{'0','5','10','15','20','25','30'}, 'FontSize', zihao-1)
set(gca,'FontName','Times New Roman')
YY=1-X'
%delete('cardatatemp.mat');
set(gca,'FontName','Times New Roman')
FinallyUAVpositions=q;
delta_t = 1; % 例如，每个时隙间隔1秒  
% 计算每个时隙的瞬时速度  
velocities = zeros(size(FinallyUAVpositions, 1) - 1, 2);
for t = 2:size(FinallyUAVpositions, 1)  
    % 计算当前位置与前一个时隙位置的差值  
    position_diff = FinallyUAVpositions(t, :) - FinallyUAVpositions(t-1, :);   
    % 计算瞬时速度（位置差除以时间间隔）  
    velocities(t-1, :) =  (position_diff) / delta_t;  
end  
Varin2=load('velocities.mat');
velocities = Varin2.velocities;  %调用了之前的数据S
FinallyUAVspeedx=  velocities(:, 1);
FinallyUAVspeedy=  velocities(:, 2);
vx = FinallyUAVspeedx;  
vy = FinallyUAVspeedy;  
vx(vx < 0) = -vx(vx < 0); 
vy(vy < 0) = -vy(vy < 0); 

% 创建一个新的图形窗口  
figure;   
% 绘制x方向速度随时间变化的曲线  
subplot(2, 1, 1); % 分割图形窗口为2行1列，当前激活第1个子图  
plot(vx, '-db', 'LineWidth', 1);  
title(' \fontname{Times New Roman}speed in the X diection');  
xlabel('\fontname{Times New Roman}T(s)');  
ylabel('\fontname{Times New Roman}speed \fontname{Times New Roman}(m/s)');  
grid on;  
set(gca,'FontName','Times New Roman')  
% 绘制y方向速度随时间变化的曲线  
subplot(2, 1, 2); % 激活第2个子图  
plot(vy, '-db', 'LineWidth', 1);  
title(' \fontname{Times New Roman}speed in the Y diection');  
xlabel('\fontname{Times New Roman}T(s)');  
ylabel('\fontname{Times New Roman}speed \fontname{Times New Roman}(m/s)');  
grid on;
set(gca,'FontName','Times New Roman')
delete('cardatatemp4.mat');
%{
figure;
hold on;
grid on;
box on;
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
xlabel('道路长度 ','FontSize',12);
ylabel('道路宽度','FontSize',12);
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
YY=1-X'
delete('cardatatemp4.mat');
%}



%{


% 绘制图形
figure;
% 选择要显示的列索引
% indices = (7:90);  
% indices = [70, 75, 80, 85, 90];
% indices = [ 75, 80, 85, 90,95];
indices = [ 56, 61, 66, 71,76];
% 提取相应列的数据
data1 = Total_theta(:, indices)';
% 在同一张图中绘制三个曲线
hold on;
plot(indices, Total_theta(1, indices), '-xb', 'LineWidth', 1.5);
plot(indices, Total_theta(2, indices), '-sr', 'LineWidth', 1.5);
plot(indices, Total_theta(3, indices), '-dg', 'LineWidth', 1.5);
% 图形标题和标签
title('能效');
xlabel('周期长度（s）');
ylabel('y');
% 设置 x 轴刻度标签
% set(gca, 'XTick', indices);
% set(gca, 'XTickLabel', {'70', '75', '80', '85', '90'});
% 图例
% legend('\theta=0.2', '\theta=0.4', '\theta=0.6');
legend('优化轨迹','固定轨迹','悬停')


figure;
indices = (7:100);  
% indices =[60 65 70 75 80];
% indices =[55 60 65 70 75];
% indices =[52 57 62 67 72];
%  indices =[50 55 60 65 70];
%  indices =[45 49 53 57 61];
% indices =[65 70 75 80 85];
indices =[70 75 80 85 90];
% indices =[75 80 85 90 95];
% indices =[ 80 85 90 95 100];
data1 = Total_theta(:, indices)';   
hold on;
plot(Total_theta(1,indices),'-xb', 'LineWidth', 1.5);    % plot(data1,'-*b'); %这个是一下子全都画出来
plot(Total_theta(2,indices),'-sr', 'LineWidth', 1.5); 
plot(Total_theta(3,indices),'-dg', 'LineWidth', 1.5); 
title('能效');  
xlabel('周期长度（s）');  
ylabel('y');
% set(gca,'XTickLabel',{'60','61','62','63','64','65','66','67'})
set(gca,'XTickLabel',{'70','75','80','85','90'})
% legend('悬停','固定轨迹','优化轨迹')
legend('\theta=0.2','\theta=0.4','\theta=0.6')


%}



% bound1=1;
% bound2=T;
% [Ybound,Xbound,GVRbound,GVUbound,sinrr,sinru,E_VRbound,E_VUbound,Total_T] = deal(zeros(M, bound2-bound1+1));
% [Pbound] = deal(zeros(bound2-bound1+1,M));
% for j=bound1:bound2
% for i=1:M
% Ybound(i,j-bound1+1)=Y(i,j-bound1+1);
% Xbound(i,j-bound1+1)=X(i,j-bound1+1);
% Pbound(j-bound1+1,i)=P_m(j-bound1+1,i);
% GVRbound(i,j-bound1+1)=GVR(i,j-bound1+1);
% GVUbound(i,j-bound1+1)=GVU(i,j-bound1+1);
% sinrr(i,j-bound1+1)= 1/log(2)*(log(Pbound(j-bound1+1,i)*GVRbound(i,j-bound1+1)+Delta)-log(Delta)) ; 
% sinru(i,j-bound1+1)= 1/log(2)*(log(Pbound(j-bound1+1,i)*GVUbound(i,j-bound1+1)+Delta)-log(Delta)) ; 
% E_VRbound(i,j-bound1+1)=Pbound(j-bound1+1,i).*Xbound(i,j-bound1+1);%VRde能量消耗
% E_VUbound(i,j-bound1+1)=Pbound(j-bound1+1,i).*Ybound(i,j-bound1+1);%VRde能量消耗
% Total_T(i,j-bound1+1)=(sinrr(i,j-bound1+1)+sinru(i,j-bound1+1))/(E_VRbound(i,j-bound1+1)+E_VUbound(i,j-bound1+1));
% Total_T(j-bound1+1)=(sum(sinrr(:,j-bound1+1))+sum(sinru(:,j-bound1+1)))/(E_VRbound(i,j-bound1+1)+E_VUbound(i,j-bound1+1));
% end
% end
% thetatheta=1
% Total_theta(thetatheta,:)=sum(Total_T,1);
% Total_TT = sum(Total_T, 1);
% 
% Total_theta=Total_TT
% 
% figure;
% indices = (62:75);  
% data1 = Total_theta(:, indices)';   
%  
% plot(data1,'-*b'); 
% hold on;
% plot(Total_theta(1,indices),'-*b', 'LineWidth', 2); 
% plot(Total_theta(2,indices),'-+r', 'LineWidth', 2); 
% plot(Total_theta(3,indices),'-*g', 'LineWidth', 2); 
% legend('\theta=0.2','\theta=0.4','\theta=0.6')
% title('能效');  
% xlabel('周期长度（s）');  
% ylabel('y');
% % set(gca,'XTickLabel',{'60','61','62','63','64','65','66','67'})
% set(gca,'XTickLabel',{'60','61','62','63','64','65','66','67'})

% figure
% indices = (50:67);  
% % plot((P_m(indices,1)),'-+r');
% grid on
% hold on
% plot((P_m(indices,2)),'-og');
% plot((P_m(indices,3)),'-*b');
% plot((P_m(indices,4)),'-sk');
% figure
% gfgg=Total_TT
% plot(gfgg)