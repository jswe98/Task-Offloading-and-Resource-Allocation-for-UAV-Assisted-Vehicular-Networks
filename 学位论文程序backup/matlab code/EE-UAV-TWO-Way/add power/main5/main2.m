clc
close all
clear
% �ܳ�ʼ��
T=70; M=8;H=100; det=1;roadwidth=30;roadlength=800; %����70��ʱ϶ %solt length   TTT=T*det;%time span
VUAV_max=15;%���˻��ٶ�Լ��
W=10^(6);
theta=0.2;%����E_VR֮ǰ�����˻��Ļ�Ƚϴ�����һ��ϵ��
Delta=1e-6;%����������
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
for bcd=1:2 %BCD��ʼ  ���Լ�while�γ�ѭ��
UAVcvx=q_l_SCA;%���������
[GVR,GVU,distanceVR,distanceVU,CARposition,L]=CARinfo(H,T,M,det,roadlength,roadwidth,UAVcvx);%ʹ����������Ϣ����
[X,Y] = Xdecision(Pm,GVU,GVR,T,M);%ʹ���˾��ߺ���

for t = 1:T %Qcvx��������Ķ���ֵ
    for m = 1:M
        c2_l(t,m) = Pm(t,m)*GVR(m,t)/Delta;
        c3_l(t,m)=Pm(t,m)/Delta;
        omega_l(t,m)=(-1*c3_l(t,m)*log2(exp(1)))/(((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,2,t))^(2)+H^(2))*...
         ((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,1,t))^(2)+H^(2)+c3_l(t,m)));
         rho_l(t,m)=log2(1+c3_l(t,m)/((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,2,t))^(2)+H^(2)));
         distance_l(t,m)=((q_l_SCA(t,1)-CARposition(m,1,t))^(2)+(q_l_SCA(t,2)-CARposition(m,2,t))^(2));
        %��l�ε�������t��ʱ϶UAV���m������ˮƽ����ƽ��
    end
 end
distance_l_SCA=distance_l;
L_VR=sum(sum(W*X'.*log2(1+Pm.*c2_l)));
E_VR=sum(sum(Pm.*X'));%VRde��������
E_VU=sum(sum(Pm.*Y'));%VUde��������

epsilon2=0.1;%����Ƕ���������ֵ
flag=0; %��ʼѭ���ı�־
chi1=0;%����ǳ���֮���ϵ��
ding_u=1;%����ѭ���˶��ٴ� %���Pm��dinkelbach�ִ� ding-p��ʾ�����Pm
while (flag==0)                  %Dinkelbach
        flag=1;                %����ѭ���ı�־
        chi1_rec(ding_u)=chi1;%chi1��record����
[q,distance2] =Qcvx(theta,E_VU,E_VR,W,Y,L_VR,rho_l,omega_l,distance_l_SCA,CARposition,VUAV_max,det,T,M);%ʹ����Qcvx����
q_l_SCA=q;%q��cvx����� variable   %��¼�����Ż������������ʹ��
            for t=1:T
                for m=1:M
                    distance2(t,m)=square_pos(norm([q(t,1) q(t,2)]-[CARposition(m,1,t) CARposition(m,2,t)]));
                end
            end
            L_A=sum(sum(W*Y'.*(rho_l+omega_l.*(distance2-distance_l_SCA))));%͹+����
            L_total=L_A+L_VR;
        F_chi1=L_total -chi1*1*(theta*sum(sum(Pm.*X'))+(1-theta)*sum(sum(Pm.*Y'))) ;         % ����Ǽ�����
        chi1=L_total/(1*(theta*sum(sum(Pm.*X'))+(1-theta)*sum(sum(Pm.*Y'))));
        ding_u=ding_u+1;                %ѭ������ÿ���Լ�����һ
         if (F_chi1>epsilon2) %�жϱ�־
             flag=0;
         end
end

chi2=0;%����ǳ���֮���ϵ��
ding_p=1;%����ѭ���˶��ٴ� %���Pm��dinkelbach�ִ� ding-p��ʾ�����Pm
flag2=0; 
while (flag2==0)                  %Dinkelbach
        flag2=1;                %����ѭ���ı�־
        chi2_rec(ding_p)=chi2;%chi1��record����
        [P_m] = Pcvx(W,X,Y,Delta,GVR,GVU,theta,T,M,chi2) ;
        Pm=P_m;
      GVRR=GVR';
      GVUU=GVU';
    for t = 1:T %�����ŵ�����h(t,m)��c2(t,m)
        for m = 1:M
            c1(t,m) = GVRR(t,m)/Delta;
            c2(t,m) =  GVUU(t,m)/Delta;%ԭ
        end
    end  %������Ϊ�����켣�滮���µ��ŵ�״̬��Ϣ�����Ż�Pm   1,1+P_b.*c1
        F_chi2=sum(sum(W*X'.*rel_entr(ones(T,M),1+P_m.*c1)))+sum(sum(W*Y'.*rel_entr(ones(T,M),1+P_m.*c2)))...
            -chi2*1*(theta*sum(sum(P_m.*X'))+(1-theta)*sum(sum(P_m.*Y'))) ;         % ����Ǽ�����
        chi2=(sum(sum(W*X'.*rel_entr(ones(T,M),1+P_m.*c1)))+sum(sum(W*Y'.*rel_entr(ones(T,M),1+P_m.*c2))))...
            /(1*(theta*sum(sum(P_m.*X'))+(1-theta)*sum(sum(P_m.*Y'))));
        ding_p=ding_p+1;                %ѭ������ÿ���Լ�����һ
         if (F_chi2>epsilon2) %�жϱ�־
             flag2=0;
         end
         
end

%���Լ���BCD��ʱ�������ж�����Ϊabs()>0.1
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
    plot(car_x(m, :), car_y(m, :), 'LineWidth', 1,'Marker', '<','MarkerFaceColor', 'auto');  % ����ÿ�����Ĺ켣
    else
    plot(car_x(m, :), car_y(m, :), 'LineWidth', 1,'Marker', '>','MarkerFaceColor', 'auto');  % ����ÿ�����Ĺ켣  
    end
end
hold off;
% �������˻�δ�Ż��ĳ�ʼ�켣
hold on;
plot(UAVX, UAVY, 'color', 'black','Marker', '+', 'LineWidth', 1.5)
hold off;
xlabel('��·����');
ylabel('��·���');
xlim = get(gca, 'XLim');
ylim = get(gca, 'YLim');
y_center = mean(ylim);
% �������˻��켣�Ż�
hold on;
plot(UAVXX, UAVYY, 'color', 'blue','Marker', '+', 'LineWidth', 1.5)
hold off;
% ��ӹ�һ����ˮƽ��ͷ
x_start = 0.5;
y_position = 0.6;    % ��ͷ�Ĵ�ֱλ��
x_length = 0.19;      % ��ͷ��ˮƽ����
y_offset = 0.1;      % ��ͷ�Ĵ�ֱλ��
hold on;
annotation('arrow', [x_start, x_start + x_length], [y_position + y_offset, y_position + y_offset], 'Units', 'normalized','color', 'black');
hold off;
 % ��ӻ�ɫ˫ʵ��
hold on;
plot(xlim, [y_center+0.05, y_center+0.05],'-',  'color', 'yellow', 'LineWidth', 2);
plot(xlim, [y_center-0.05, y_center-0.05],'-', 'color', 'yellow', 'LineWidth', 2);
hold off;
legend_entries{M+1} = '���˻���ʼ�켣';
legend_entries{M+2} = 'Optimization UAV Trajectory';
legend_entries{M+3} = '��ɫ˫ʵ��';
% legend_entries{M+3} = '����Ϊ������';
%legend_entries{M+4} = '��ɫ˫ʵ��';
legend(legend_entries, 'Location', 'best');
XX=1-X'


