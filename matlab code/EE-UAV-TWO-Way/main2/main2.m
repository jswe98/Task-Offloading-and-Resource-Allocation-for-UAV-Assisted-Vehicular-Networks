clc
close all
clear
% �ܳ�ʼ��
T=100; M=5;H=80; det=1;roadwidth=8;roadlength=800; %����70��ʱ϶ %solt length   TTT=T*det;%time span
VUAV_max=10;%���˻��ٶ�Լ��
W=10^(6);
Delta=1e-6;%����������
Pm=0.03*ones(T,M);
R=zeros(T,M);
[UAVposition] = InitializationUAVposition(T,roadlength,roadwidth);
UAVcvx=UAVposition;
[GVR,GVU,distanceVR,distanceVU,CARposition,L]=CARinfo(H,T,M,det,roadlength,roadwidth,UAVcvx);
[X,Y] = Xdecision(Pm,GVU,GVR,T,M);

q_l_SCA=UAVposition;
for t = 1:T %Qcvx��������Ķ���ֵ
    for m = 1:M
%         c2(t,m) = 1/(l_s2b(1,m)^(kappa)*af^(l_s2b(1,m))*Nf*deltaf);%ԭ
        c3_l(t,m)=Pm(t,m)/Delta;
        omega_l(t,m)=(-1*c3_l(t,m)*log2(exp(1)))/(((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,2,t))^(2)+H^(2))*...
         ((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,1,t))^(2)+H^(2)+c3_l(t,m)));
         rho_l(t,m)=log2(1+c3_l(t,m)/((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,2,t))^(2)+H^(2)));
         distance_l(t,m)=((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,2,t))^(2));
        %��l�ε�������t��ʱ϶UAV���m������ˮƽ����ƽ��
    end
 end
distance_l_SCA=distance_l;
[q,distance2] =Qcvx(W,Y,rho_l,omega_l,distance_l_SCA,CARposition,...
                          VUAV_max,det,T,M)









figure;
hold on;
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
           legend_entries{m} = ['leftCar ' num2str(m)]; 
        else
           legend_entries{m} = ['rightCar ' num2str(m)]; 
        end
    end
end
for m = 1:M
    if L(m)==-1;
    plot(car_x(m, :), car_y(m, :), 'LineWidth', 1,'Marker', '<','MarkerFaceColor', 'auto');  % ����ÿ�����Ĺ켣
    else
    plot(car_x(m, :), car_y(m, :), 'LineWidth', 1,'Marker', '>','MarkerFaceColor', 'auto');  % ����ÿ�����Ĺ켣  
    end
end
hold off;
hold on;
plot(UAVX, UAVY, 'color', 'black','Marker', '+', 'LineWidth', 1.5)
hold off;
xlabel('��·����');
ylabel('��·���');
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
y_center = mean(ylim);
hold on;
plot(xlim, [y_center+0.05, y_center+0.05],'-',  'color', 'yellow', 'LineWidth', 2);
plot(xlim, [y_center-0.05, y_center-0.05],'-', 'color', 'yellow', 'LineWidth', 2);
hold off;
% ��ӹ�һ����ˮƽ��ͷ
x_start = 0.5;
y_position = 0.2;    % ��ͷ�Ĵ�ֱλ��
x_length = 0.19;      % ��ͷ��ˮƽ����
y_offset = 0.1;      % ��ͷ�Ĵ�ֱλ��
annotation('arrow', [x_start, x_start + x_length], [y_position + y_offset, y_position + y_offset], 'Units', 'normalized','color', 'black');
legend_entries{M+1} = '���˻�';
legend_entries{M+2} = '��ɫ˫ʵ��';
% legend_entries{M+3} = '����Ϊ������';
hold on;
plot(UAVXX, UAVYY, 'color', 'blue','Marker', '+', 'LineWidth', 1.5)
hold off;
legend_entries{M+3} = '�Ż������˻��켣';
legend(legend_entries, 'Location', 'best');


