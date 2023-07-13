clc ;                                                                                                                                                                         %����������⣬��ϣ��������Ի���ת�⣬��ʱֻ����Ѱ����    2021/5/26                                    
clear ;
close all;
%% ���߰���Ķ�������һ��
M=3;%�û�����
delta=1;%��Сʱ϶ʱ��
Euav=150;%�����Ҫ�ٵ������˻��������
mu=9.999999999999972e-10;%IRS�ܺ�-60dBm ����Ԫ����
N=50;%ʱ϶����
K=90;%IRSԪ������90
eta=0.8;%��������Ч��
P0=10;%���˻����书��
Vmax=50;%���˻���������ٶ�
kappa=2;%��˹ָ��3dB
alpha=2.8;%IRS���û�֮���·��
beta=3.2;%UAV���û�֮���·��
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
%�Աȱ���
Dtaylor=zeros(8);
Doriginal=zeros(8);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ʼ�ŵ�����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% UAV��IRS���ŵ�����
%UAV��IRS��ʼ����
Drv=zeros(1,N);
for n=1:N
    Drv(n)=norm([Qv(1,n) Qv(2,n) hv]-[Qr(1) Qr(2) hr]);
end
%UAV��IRS·��
pathloss1=zeros(1,N);
for n=1:N
    pathloss1(n)=sqrt(beta0/square_pos(Drv(n)));%*exp(-1j*2*pi/lambda*Drv(n));ȥ��������ɵ���λ
end
%UAV��IRS������Ӧ
array1=zeros(K,N);
phi_rv=zeros(1,N);%UAV��IRS���������ֵ
for n=1:N
    phi_rv(n)=(Qr(1)-Qv(1,n))/Drv(n);
end
for n=1:N
    for k=1:K
        array1(k,n)=exp(-1j*2*pi/lambda*d*(k-1)*phi_rv(n));
    end
end
%�洢һ������
%save('array1.mat','array1') 
%UAV��IRS�ŵ�
Hrv=zeros(K,N);
for n=1:N
    Hrv(:,n)=pathloss1(n)*array1(:,n);
end
%test
% a=norm(Hrv(1,3));
% b=norm(Hrv(2,3));
%% UAV��IRS���ŵ�����
%IRS���û��ľ���
Dru=zeros(1,M);
for m=1:M
    Dru(m)=norm([Qu(1,m) Qu(2,m) 0]-[Qr(1) Qr(2) hr]);
end
%IRS��User_m·��
pathloss2=zeros(1,M);
for m=1:M
    pathloss2(m)=sqrt(beta0/pow_pos(Dru(m),alpha));%*exp(-1j*2*pi/lambda*Dru(m));ȥ��������ɵ���λ
end
% array2=zeros(K,M);%IRS���û�������Ӧ
 %phi_ru=zeros(1,M);%IRS���û����������ֵ
 %for m=1:M
%     phi_ru(m)=(Qu(1,m)-Qr(1))/Dru(m);
 %end
% for m=1:M
%     for k=1:K
%        %ȥ���˾�����ɵ���λ 
%        %array2(k,m)=sqrt(kappa/(1+kappa))*exp(-1j*2*pi/lambda*Dru(m))*exp(-1j*2*pi/lambda*d*(k-1)*phi_ru(m))+sqrt(1/(1+kappa))*sqrt(1/2)*(randn(1,1)+1j*randn(1,1));
%         array2(k,m)=sqrt(kappa/(1+kappa))*exp(-1j*2*pi/lambda*d*(k-1)*phi_ru(m))+sqrt(1/(1+kappa))*sqrt(1/2)*(randn(1,1)+1j*randn(1,1));
%     end
% end
%% %%��ά����
array2=zeros(K,M,N);%IRS���û�������Ӧ
phi_ru=zeros(1,M);%IRS���û����������ֵ
for m=1:M
    phi_ru(m)=(Qu(1,m)-Qr(1))/Dru(m);
end
%% ����array2
% for n=1:N
%     for m=1:M
%         for k=1:K
%            %ȥ���˾�����ɵ���λ 
%             %array2(k,m)=sqrt(kappa/(1+kappa))*exp(-1j*2*pi/lambda*Dru(m))*exp(-1j*2*pi/lambda*d*(k-1)*phi_ru(m))+sqrt(1/(1+kappa))*sqrt(1/2)*(randn(1,1)+1j*randn(1,1));
%            array2(k,m,n)=sqrt(kappa/(1+kappa))*exp(-1j*2*pi/lambda*d*(k-1)*phi_ru(m))+sqrt(1/(1+kappa))*sqrt(1/2)*(randn(1,1)+1j*randn(1,1));
%         end
%     end
% end
% save('array2.mat','array2')
array_22=load('array2.mat');%����������˻�λ�õ�����ʱA��ֵ
array2=array_22.array2;

%%
%IRS���û��ŵ�
% Hru=zeros(K,M);
% for m=1:M
%     Hru(:,m)=pathloss2(m)*array2(:,m,);
% end
Hru=zeros(K,M,N);
for n=1:N
    for m=1:M
        for k=1:K
            %Hru(:,m,n)=pathloss2(m)*array2(:,m,n);
            Hru(k,m,n)=pathloss2(m)*array2(k,m,n);
        end
    end
end
save('Hru.mat','Hru') 
%text
 %a=norm(Hru(3,2));
 %b=norm(Hru(4,2));
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
%% ����array3
array3=zeros(M,N);
% for n=1:N
%     for m=1:M
%         array3(m,n)=sqrt(1/2)*(randn(1,1)+1j*randn(1,1));
%     end
% end
% save('array3.mat','array3') 
array_33=load('array3.mat');%����������˻�λ�õ�����ʱA��ֵ
array3=array_33.array3;
Hvu=zeros(M,N);
for n=1:N
    for m=1:M
        Hvu(m,n)=pathloss3(m,n);%*array3(m,n);
    end
end
%% ���BCD��������ΪL�����˻��켣�ĵ�������ΪZ
%����һЩ��Ҫ�ĳ�ʼ��
L=10;%BCD��ѭ������
Rit=zeros(1,L);%���ڵõ������Ľ��
h1=zeros(M,N);%����ÿ�θ���direct���棬�����ǳ�ʼ��
h1(:,:)=Hvu(:,:);
% h2=zeros(K,M,N);%����ÿ�θ���R-U���棬�����ǳ�ʼ��
% h2(:,:,:)=Hru(:,:,:);
Pu=zeros(M,N,L);
for m=1:M
    Pu(m,:,1)=0.001;
end
C1=zeros(M,N,L);
C2=zeros(M,N,L);
d1=zeros(M,N,L);
d1(:,:,1)=Dvu(:,:);
d2=zeros(N,L);
d2(:,1)=Drv(:);
hz1=zeros(M,N,L);%���ڵڶ���ģ�����ʱ�����
hz2=zeros(M,N,L);%���ڵڶ���ģ�����ʱ�����
hz0=zeros(M,N,L);%���ڵڶ���ģ�����ʱ����ȵ�ƽ��
a0=Dvu;
b0=Drv;
array_vu=load('array3.mat');%����������˻�λ�õ�����ʱA��ֵ
array_vu=array_vu.array3;
for n=1:N
    for m=1:M
        array_vu(m,n)=abs(array_vu(m,n));
    end
end
Hru_ru=load('Hru.mat');%����������˻�λ�õ�����ʱB��ֵ
Hru_ru=Hru_ru.Hru;
h2=zeros(K,M,N);%����ÿ�θ���R-U���棬�����ǳ�ʼ��
h2(:,:,:)=Hru(:,:,:);
for n=1:N
    for m=1:M
        for k=1:K
            Hru_ru(k,m,n)=abs(Hru_ru(k,m,n));
        end
    end
end
%����A
for n=1:N
   for m=1:M
       A(m,n)=sqrt(beta0)*array_vu(m,n);
   end
end
%����B
for n=1:N
   for m=1:M
       B(m,n)=sqrt(beta0)*sum(Hru_ru(:,m,n));
   end
end
%% BCD�㷨
for l=1:L
%% ��һ���Ż�������Ϣ������λ
for n=1:N
    for m=1:M
         %C1(m,n,l)=sqrt(beta0)*abs(h1(m,n));
         C1(m,n,l)=sqrt(beta0)*abs(array_vu(m,n));
    end
end

for n=1:N
    for m=1:M
        for k=1:K
             C2(m,n,l)=C2(m,n,l)+sqrt(beta0)*abs(h2(k,m,n));
        end
    end
end
%% �ڶ����Ż�����ʱ��
%��ʱ�ŵ�����hz
for n=1:N
    for m=1:M
        hz1(m,n,l)=C1(m,n,l)/pow_pos(d1(m,n,l),beta/2);
    end
end
for n=1:N
    for m=1:M
         hz2(m,n,l)=C2(m,n,l)/d2(n,l);
    end
end
for n=1:N
    for m=1:M
         hz0(m,n,l)=pow_pos((hz1(m,n,l)+hz2(m,n,l)),2);
    end
end
cvx_solver SDPT3%��ʼ�ܴ���ѡ�������ʼ��С��ѡ������
%cvx_solver SeDuMi
%cvx_precision best
cvx_begin 
           variable Em(M,N)                %�Դ���������
           variable tm(M,N)                %�Դ�ʱ�䶨��
           variable te(N)                  %��������ʱ�䶨��
           expression Eth(M,N)             %��ʱÿ��ʱ϶�û���������������
           for n=1:N
               for m=1:M
                   Eth(m,n)=eta*P0*te(n)*pow_p(abs(h1(m,n)),2);
               end
           end
           maximize(sum(sum(-rel_entr(tm,tm+Em.*hz0(:,:,l)/sigma2)))/log(2));
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
           mu*sum(sum(tm))<=eta*beta0*P0*sum(te.*pow_p(d2(n,l),-2));
           %mu*sum(tm(:,n))<=eta*P0*te(n)*beta0*pow_p(d2(n,l),-2);
           %mu*sum(tm(:,n))<=eta*P0*te(n)*beta0*inv_pos(pow_pos(d2(n,l),2));
%          end
           Em>=0;
           tm>=0;
           te>=0;   
cvx_end  
%text
pm=Em./tm;
zz=eta*beta0*P0*sum(te.*pow_p(d2(n,l),-2))-mu*sum(sum(tm));
gg=te(n)+sum(tm(:,n));
xx=sum(Em,2)-sum(Eth,2);
ZZ=sum(sum(tm.*log(1+Em.*hz0(:,:,l).*inv_pos(sigma2*tm))));
%TEST��ʽ�Ƿ���ȷ
% pp=Em./tm;
% AA=sum(sum(-rel_entr(tm,tm+Em.*hz0(:,:,l)/sigma2)));
% BB=sum(sum(tm.*log(1+Em.*hz0(:,:,l).*inv_pos(sigma2*tm))));
% CC=0;
% for m=1:M
%     for n=1:N
%         CC=CC+tm(m,n)*log(1+Em(m,n)*hz0(m,n,l)*inv_pos(sigma2*tm(m,n)));
%     end
% end
%% �������Ż����˻����й켣
%���¼������UAV-IRS��UAV-USERS
%UAV-IRS Drv_tΪ��ʱ����
for i=1:3
%R=0;
%����M0
for n=1:N
    for m=1:M
        M0(m,n)=1+pm(m,n)/sigma2*(pow_p(A(m,n),2)/pow_p(a0(m,n),beta)+pow_p(B(m,n),2)/pow_p(b0(n),2)+2*A(m,n)*B(m,n)/(pow_p(a0(m,n),beta/2)*b0(n)));
    end
end
%����M1
for n=1:N
    for m=1:M
        M1(m,n)=-pm(m,n)/sigma2*(beta*pow_p(A(m,n),2)/pow_p(a0(m,n),beta+1)+beta*A(m,n)*B(m,n)/(pow_p(a0(m,n),beta/2+1)*b0(n)));
    end
end
%����M2
for n=1:N
    for m=1:M
        M2(m,n)=-pm(m,n)/sigma2*(2*pow_p(B(m,n),2)/pow_p(b0(n),3)+2*A(m,n)*B(m,n)/(pow_p(a0(m,n),beta/2)*pow_p(b0(n),2)));
    end
end
%cvx_clear
cvx_solver SDPT3
cvx_begin 
           variable q_cvx(2,N)            %���˻�λ�ö���
           variable a(M,N)                %���˻����û�������ɳڱ���
           variable b(N)                  %���˻���IRS������ɳڱ���
           expression R(M,N);
           expression Gth(M,N)
           expression Gth2(N)
 for n=1:N
    for m=1:M
        R(m,n)=log(M0(m,n))/log(2)+ M1(m,n)/(M0(m,n)*log(2))*(a(m,n)-a0(m,n))+M2(m,n)/(M0(m,n)*log(2))*(b(n)-b0(n));
       %R=R+log(M0(m,n))+ M1(m,n)/(M0(m,n)*log(2))*(a(m,n)-a0(m,n))+M2(m,n)/(M0(m,n)*log(2))*(b(n)-b0(n));
    end
 end
R=tm.*R;
 for n=1:N
    for m=1:M
        Gth(m,n)=te(n)*(pow_p(abs(array_vu(m,n)),2)*(beta0/pow_p(a0(m,n),beta)-beta0*beta/pow_p(a0(m,n),beta+1)*(a(m,n)-a0(m,n))));
    end
 end
  for n=1:N
      Gth2(n)=te(n)*(beta0/pow_p(b0(n),2)-beta0*2/pow_p(b0(n),3)*(b(n)-b0(n)));
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
           %�ɳڱ���b������Լ��
           for n=1:N
               square_pos(norm([q_cvx(1,n) q_cvx(2,n) hv]-[Qr(1) Qr(2) hr]))+square_pos(b0(n))-2*b0(n)*b(n)<=0;
           end          
           %����ʱ��֮���Լ��
           for n=1:N-1
               norm([q_cvx(1,n) q_cvx(2,n)]-[q_cvx(1,n+1) q_cvx(2,n+1)])<=Vmax*delta;
           end
           sum(Em,2)<=eta*P0*sum(Gth,2);
           mu*sum(sum(tm))<=sum(Gth2);                    
%            %ÿ��n�ڵ�Լ�����������ʵ��
%             for n=1:N
%                 for m=1:M 
%                     Em(m,n)<=eta*P0*te(n)*abs(array3(m,n))*(beta0/pow_p(a0(m,n),beta)-beta0*beta/pow_p(a0(m,n),beta+1)*(a(m,n)-a0(m,n)));
%                 end
%             end            
           q_cvx(1,1)==0;
           q_cvx(2,1)==0;
           q_cvx(1,N)==500;
           q_cvx(2,N)==500;
cvx_end  
% %̩����ԭ�����Ա�
% Dtaylor(i)=sum(sum(R));
% Doriginal(i)=sum(sum(tm.*log(1+Em.*hz0(:,:,l).*inv_pos(sigma2*tm))));

%�����ɳڱ�������͹�ƽ���ĸ���
           for n=1:N
               for m=1:M
                   a0(m,n)=a(m,n);
               end
           end
           for n=1:N
                   b0(n)=b(n);
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

end
%% refresh �����������ö����������������������
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
        h1(m,n)=pathloss3(m,n)*array_vu(m,n);
    end
end
for n=1:N
    for m=1:M
        d1(m,n,l+1)=Dvu(m,n);
    end
end
for n=1:N
    d2(n,l+1)=norm([q_cvx(1,n) q_cvx(2,n) hv]-[Qr(1) Qr(2) hr]);
end

d1(:,:,1)
figure(2);%����settingͼ
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
ccc=zeros(M,N);
for n=1:N%�ŵ�
    for m=1:M
        ccc(m,n)=norm([q_cvx(1,n) q_cvx(2,n) hv]-[Qu(1,m) Qu(2,m) 0]);
    end
end
Rit(l)=sum(sum(tm.*log(1+Em.*hz0(:,:,l).*inv_pos(sigma2*tm))));
%̩����ԭ�����Ա�
Dtaylor(l)=sum(sum(R));
Doriginal(l)=sum(sum(tm.*log(1+Em.*hz0(:,:,l).*inv_pos(sigma2*tm))));
end

figure(3);%����settingͼ
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

figure(4);
plot(Rit,'-^g','linewidth',1.5);
xlabel('Iteration');
ylabel('sum-rate (bps/Hz)');
%Rit90=[1.24142902216963,22.2968087951863,24.2076365686919,24.3257009310573,24.3424981714509];%N=50
%Rit50=[1.24139952969553,13.8092605492158,20.2063584187745,20.4828618620443,20.4897106540524]
%Rnophase90=[1.24121462667579,12.9964157715700,18.6694812352935,19.8024329072787,19.6855841001474]%N=50
%Rnophase50=[1.24105015912072,12.9335495319730,17.7082130947848,18.4359893594171,18.3765949668826]

figure(5);
v=zeros(1,N);
for n=1:N-1
    v(n)= norm([q_cvx(1,n) q_cvx(2,n)]-[q_cvx(1,n+1) q_cvx(2,n+1)]);
end
plot(v,'-ob','linewidth',1.5);
hold on
plot([0 1],[0 49.9680392014224],'-ob','linewidth',1.5 );
xlabel('Time slot');
ylabel('UAV flying speed (m/s)');
grid on
axis([0 50 0 55])
% p=zeros(M,N);%IRSλ�õĹ��������
% p=Em./tm;
% figure(6);
% plot(p(1,:),'-^r','linewidth',1.5);
% hold on
% plot(p(2,:),'-^b','linewidth',1.5);
% hold on
% plot(p(3,:),'-^g','linewidth',1.5);
