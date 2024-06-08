clc ;                                                                                                                                                                         %如果爱有天意，我希望望舒可以回心转意，当时只道是寻常。    2021/5/26                                    
clear ;
close all;
%% 乱七八糟的东西设置一下
M=3;%用户个数
delta=1;%最小时隙时间
Euav=150;%这个还要再调整无人机电池容量
mu=9.999999999999972e-10;%IRS能耗-60dBm 单个元件的
N=20;%时隙个数
K=50;%IRS元件个数90
eta=0.8;%能量吸收效率
P0=10;%无人机发射功率
Vmax=50;%无人机飞行最大速度
kappa=2;%莱斯指数3dB
alpha=2.8;%IRS与用户之间的路损
beta=3.5;%UAV与用户之间的路损
lambda=0.1;%波长
d=lambda/2;%元件间距
hv=30;%无人机高度
hr=15;%IRS高度
Qu=zeros(2,M);%用户坐标
% Qu(:,1)=[50,250]';
% Qu(:,2)=[150,450]';
% Qu(:,3)=[400,100]';
Qu(:,1)=[300,450]';
Qu(:,2)=[400,300]';
%Qu(:,2)=[200,450]';
Qu(:,3)=[100,50]';
Qr=[0,350]';%IRS坐标
Qv=zeros(2,N);%UAV初始化坐标
Qv(:,1)=[0,0]';
Qv(:,N)=[500,500]';
for n=2:N-1
    Qv(:,n)=Qv(:,n-1)+[500/(N-1),500/(N-1)]';
end
figure(1);%出个setting图
plot(Qr(1),Qr(2),'ob','markersize',10,'MarkerFaceColor','r')
hold on
plot(Qu(1,1),Qu(2,1),'pb','markersize',10,'MarkerFaceColor','b')
hold on
plot(Qu(1,2),Qu(2,2),'pb','markersize',10,'MarkerFaceColor','b')
hold on
plot(Qu(1,3),Qu(2,3),'pb','markersize',10,'MarkerFaceColor','b')
hold on
plot(Qv(1,:),Qv(2,:),'-^c','linewidth',1);
axis([0 500 0 500])% 设置图像边界
legend('IRS','User1','User2','User3');
xlabel('x/m');
ylabel('y/m');
beta0=10^(-5);%单位信道增益beta0=10^(-5)的结果是0.3左右，-4的结果是30左右
sigma2=3.9810717055349565e-21;%噪声功率-174dBm -25

%% UAV到用户的信道
Dvu=zeros(M,N);
for n=1:N
    for m=1:M
        Dvu(m,n)=norm([Qv(1,n) Qv(2,n) hv]-[Qu(1,m) Qu(2,m) 0]);
    end
end
pathloss3=zeros(M,N);
for n=1:N
    for m=1:M
        pathloss3(m,n)=sqrt(beta0/pow_pos(Dvu(m,n),beta));
    end
end
array3=zeros(M,N);
array_33=load('array3.mat');%定义计算无人机位置的坐标时A的值
array3=array_33.array3;
Hvu=zeros(M,N);
for n=1:N
    for m=1:M
        Hvu(m,n)=pathloss3(m,n)*array3(m,n);
    end
end
array_vu=load('array3.mat');%定义计算无人机位置的坐标时A的值
array_vu=array_vu.array3;
for n=1:N
    for m=1:M
        array_vu(m,n)=abs(array_vu(m,n));
    end
end
for n=1:N
   for m=1:M
       A(m,n)=sqrt(beta0)*array_vu(m,n);
   end
end
a0=Dvu;
%% 大的BCD迭代次数为L。无人机轨迹的迭代次数为Z
L=5;%BCD大循环次数
h1=zeros(M,N,L);
Rit=zeros(1,L);%用于得到迭代的结果
for n=1:N
    for m=1:M
         h1(m,n,1)=abs(Hvu(m,n));
    end
end
for l=1:L
%进行一些必要的初始化
%cvx_solver SDPT3
cvx_solver SeDuMi
cvx_precision best
cvx_begin 
           variable Em(M,N)                %自传能量定义
           variable tm(M,N)                %自传时间定义
           variable te(N)                  %能量吸收时间定义
           expression Eth(M,N)             %临时每个时隙用户的吸收能量定义
           for n=1:N
               for m=1:M
                   Eth(m,n)=eta*P0*te(n)*pow_p(h1(m,n,l),2);
               end
           end
           maximize(sum(sum(-rel_entr(tm,tm+Em.*pow_p(h1(:,:,l),2)/sigma2)))/log(2));
           subject to
%          for n=1:N
%                for m=1:M 
%                    sum(Em(m,:))<=sum(Eth(m,:));
%                   Em(m,n)<=Eth(m,n);
%                end
%          end
                sum(Em,2)<=sum(Eth,2);
           for n=1:N
               te(n)+sum(tm(:,n))<=delta;
           end
%          for n=1:N
          % mu*sum(sum(tm))<=eta*beta0*P0*sum(te.*pow_p(d2(n,l),-2));
           %mu*sum(tm(:,n))<=eta*P0*te(n)*beta0*pow_p(d2(n,l),-2);
           %mu*sum(tm(:,n))<=eta*P0*te(n)*beta0*inv_pos(pow_pos(d2(n,l),2));
%          end
           Em>=0;
           tm>=0;
           te>=0;   
cvx_end  

pm=Em./tm;

%% 轨迹优化
for i=1:3
%R=0;
%定义M0
for n=1:N
    for m=1:M
        M0(m,n)=1+pm(m,n)/sigma2*(pow_p(A(m,n),2)/pow_p(a0(m,n),beta));
    end
end
%定义M1
for n=1:N
    for m=1:M
        M1(m,n)=-pm(m,n)/sigma2*(beta*pow_p(A(m,n),2)/pow_p(a0(m,n),beta+1));
    end
end

cvx_solver SDPT3
cvx_begin 
           variable q_cvx(2,N)            %无人机位置定义
           variable a(M,N)                %无人机到用户距离的松弛变量
           expression R(M,N);
           expression Gth(M,N)
 for n=1:N
    for m=1:M
        R(m,n)=log(M0(m,n))/log(2)+ M1(m,n)/(M0(m,n)*log(2))*(a(m,n)-a0(m,n));
       %R=R+log(M0(m,n))+ M1(m,n)/(M0(m,n)*log(2))*(a(m,n)-a0(m,n))+M2(m,n)/(M0(m,n)*log(2))*(b(n)-b0(n));
    end
 end
R=tm.*R;
 for n=1:N
    for m=1:M
        Gth(m,n)=te(n)*(pow_p(abs(array_vu(m,n)),2)*(beta0/pow_p(a0(m,n),beta)-beta0*beta/pow_p(a0(m,n),beta+1)*(a(m,n)-a0(m,n))));
    end
 end
           maximize(sum(sum(R)));
          %maximize(R/log(2));
           subject to
           %松弛变量a与距离的约束
           for n=1:N
               for m=1:M
                   square_pos(norm([q_cvx(1,n) q_cvx(2,n) hv]-[Qu(1,m) Qu(2,m) 0]))+square_pos(a0(m,n))-2*a0(m,n)*a(m,n)<=0;
               end
           end
           %两个时刻之间的约束
           for n=1:N-1
               norm([q_cvx(1,n) q_cvx(2,n)]-[q_cvx(1,n+1) q_cvx(2,n+1)])<=Vmax*delta;
           end
           sum(Em,2)<=eta*P0*sum(Gth,2);
           q_cvx(1,1)==0;
           q_cvx(2,1)==0;
           q_cvx(1,N)==500;
           q_cvx(2,N)==500;
cvx_end  
           for n=1:N
               for m=1:M
                  a0(m,n)=a(m,n);
                 %a0(m,n)=norm([q_cvx(1,n) q_cvx(2,n) hv]-[Qu(1,m) Qu(2,m) 0])
               end
           end
end
%% refresh
%信道更新   
%UAV到用户的信道
for n=1:N%距离
    for m=1:M
        Dvu(m,n)=norm([q_cvx(1,n) q_cvx(2,n) hv]-[Qu(1,m) Qu(2,m) 0]);
    end
end
for n=1:N%路损
    for m=1:M
        pathloss3(m,n)=sqrt(beta0/pow_pos(Dvu(m,n),beta));
    end
end
for n=1:N%信道
    for m=1:M
        h1(m,n,l+1)=pathloss3(m,n)*array_vu(m,n);
    end
end        


figure(1);%出个setting图
plot(Qr(1),Qr(2),'ob','markersize',10,'MarkerFaceColor','r')
hold on
plot(Qu(1,1),Qu(2,1),'pb','markersize',10,'MarkerFaceColor','b')
hold on
plot(Qu(1,2),Qu(2,2),'pb','markersize',10,'MarkerFaceColor','b')
hold on
plot(Qu(1,3),Qu(2,3),'pb','markersize',10,'MarkerFaceColor','b')
hold on
plot(Qv(1,:),Qv(2,:),'-^c','linewidth',1);
hold on
linehandle=plot( q_cvx(1,:), q_cvx(2,:),'-or','linewidth',1);
set( linehandle, 'linesmoothing', 'on' );
axis([0 500 0 500])% 设置图像边界
legend('IRS','User1','User2','User3');
xlabel('x/m');
ylabel('y/m');
Rit(l)=sum(sum(tm.*log(1+Em.*pow_p(h1(:,:,l),2).*inv_pos(sigma2*tm))));
end
save('q_cvx.mat','q_cvx') 
figure(4);
plot(Rit,'-^g','linewidth',1.5);
xlabel('Iteration');
ylabel('sum-rate (bps/Hz)');