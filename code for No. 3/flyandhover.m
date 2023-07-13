clc ;                                                                                                                                                                         %����������⣬��ϣ��������Ի���ת�⣬��ʱֻ����Ѱ����    2021/5/26                                    
clear ;
close all;
%% ���߰���Ķ�������һ��
M=3;%�û�����
delta=1;%��Сʱ϶ʱ��
Euav=150;%�����Ҫ�ٵ������˻��������
mu=9.999999999999972e-10;%IRS�ܺ�-60dBm ����Ԫ����
N=20;%ʱ϶����
K=50;%IRSԪ������90
eta=0.8;%��������Ч��
P0=10;%���˻����书��
Vmax=50;%���˻���������ٶ�
kappa=2;%��˹ָ��3dB
alpha=2.8;%IRS���û�֮���·��
beta=3.5;%UAV���û�֮���·��
lambda=0.1;%����
d=lambda/2;%Ԫ�����
hv=30;%���˻��߶�
hr=15;%IRS�߶�
Qu=zeros(2,M);%�û�����
% Qu(:,1)=[50,250]';
% Qu(:,2)=[150,450]';
% Qu(:,3)=[400,100]';
Qu(:,1)=[300,450]';
Qu(:,2)=[400,300]';
%Qu(:,2)=[200,450]';
Qu(:,3)=[100,50]';
Qr=[0,350]';%IRS����
Qv=zeros(2,N);%UAV��ʼ������
Qv(:,1)=[0,0]';
Qv(:,N)=[500,500]';
for n=2:N-1
    Qv(:,n)=Qv(:,n-1)+[500/(N-1),500/(N-1)]';
end
figure(1);%����settingͼ
plot(Qr(1),Qr(2),'ob','markersize',10,'MarkerFaceColor','r')
hold on
plot(Qu(1,1),Qu(2,1),'pb','markersize',10,'MarkerFaceColor','b')
hold on
plot(Qu(1,2),Qu(2,2),'pb','markersize',10,'MarkerFaceColor','b')
hold on
plot(Qu(1,3),Qu(2,3),'pb','markersize',10,'MarkerFaceColor','b')
hold on
plot(Qv(1,:),Qv(2,:),'-^c','linewidth',1);
axis([0 500 0 500])% ����ͼ��߽�
legend('IRS','User1','User2','User3');
xlabel('x/m');
ylabel('y/m');
beta0=10^(-5);%��λ�ŵ�����beta0=10^(-5)�Ľ����0.3���ң�-4�Ľ����30����
sigma2=3.9810717055349565e-21;%��������-174dBm -25

%% UAV���û����ŵ�
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
array_33=load('array3.mat');%����������˻�λ�õ�����ʱA��ֵ
array3=array_33.array3;
Hvu=zeros(M,N);
for n=1:N
    for m=1:M
        Hvu(m,n)=pathloss3(m,n)*array3(m,n);
    end
end
array_vu=load('array3.mat');%����������˻�λ�õ�����ʱA��ֵ
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
%% ���BCD��������ΪL�����˻��켣�ĵ�������ΪZ
L=5;%BCD��ѭ������
h1=zeros(M,N,L);
Rit=zeros(1,L);%���ڵõ������Ľ��
for n=1:N
    for m=1:M
         h1(m,n,1)=abs(Hvu(m,n));
    end
end
for l=1:L
%����һЩ��Ҫ�ĳ�ʼ��
%cvx_solver SDPT3
cvx_solver SeDuMi
cvx_precision best
cvx_begin 
           variable Em(M,N)                %�Դ���������
           variable tm(M,N)                %�Դ�ʱ�䶨��
           variable te(N)                  %��������ʱ�䶨��
           expression Eth(M,N)             %��ʱÿ��ʱ϶�û���������������
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

%% �켣�Ż�
for i=1:3
%R=0;
%����M0
for n=1:N
    for m=1:M
        M0(m,n)=1+pm(m,n)/sigma2*(pow_p(A(m,n),2)/pow_p(a0(m,n),beta));
    end
end
%����M1
for n=1:N
    for m=1:M
        M1(m,n)=-pm(m,n)/sigma2*(beta*pow_p(A(m,n),2)/pow_p(a0(m,n),beta+1));
    end
end

cvx_solver SDPT3
cvx_begin 
           variable q_cvx(2,N)            %���˻�λ�ö���
           variable a(M,N)                %���˻����û�������ɳڱ���
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
           %�ɳڱ���a������Լ��
           for n=1:N
               for m=1:M
                   square_pos(norm([q_cvx(1,n) q_cvx(2,n) hv]-[Qu(1,m) Qu(2,m) 0]))+square_pos(a0(m,n))-2*a0(m,n)*a(m,n)<=0;
               end
           end
           %����ʱ��֮���Լ��
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
%�ŵ�����   
%UAV���û����ŵ�
for n=1:N%����
    for m=1:M
        Dvu(m,n)=norm([q_cvx(1,n) q_cvx(2,n) hv]-[Qu(1,m) Qu(2,m) 0]);
    end
end
for n=1:N%·��
    for m=1:M
        pathloss3(m,n)=sqrt(beta0/pow_pos(Dvu(m,n),beta));
    end
end
for n=1:N%�ŵ�
    for m=1:M
        h1(m,n,l+1)=pathloss3(m,n)*array_vu(m,n);
    end
end        


figure(1);%����settingͼ
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
axis([0 500 0 500])% ����ͼ��߽�
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