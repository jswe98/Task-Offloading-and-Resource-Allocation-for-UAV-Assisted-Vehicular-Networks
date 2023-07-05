clc ;                                                                                           
clear ;
close all;
global K
global M1
global M2
global N
% global eta
% global repeat
%% parameters setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(5)% fix random results
K=1;%the number of users
J=1;%the number of eves
M1=40;%the size of IRS1
M2=M1;%the size of IRS2
N=5;%the number of AP's antennas
P0=db2pow(15);%transmit power of AP
%P = db2pow(-5:5:25); % 发射功率15dBm
eta=1;%背景噪声
sigmaK2 = db2pow(-88); 
repeat=1000;%Gaussian_random times
%% channel setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L0=30;%log distance path loss 30dB 10^(-L0/10)
C0 = db2pow(-30); % 参考距离时的路损
%AP coordinate setup
x_ap=0;
y_ap=0;
z_ap=0;
%IRS coordinate setup
x_irs1=10;
y_irs1=-5;
z_irs1=15;
x_irs2=40;
y_irs2=-5;
z_irs2=15;
%user coordinate setup
for k=1:K
    x_user(k)=50;
    y_user(k)=0;
    z_user(k)=0;
end
% y_user=linspace(49,51,K);
%  x_user=[47.5 48.75 50 51.25 52.5 ];
% y_user=[49.2929 49.6173 ];
% x_user=[0.7071 0.9239 0.3827 0 1 ];
% y_user=[49.2929 49.6173 49.0761 49 50];
%Eve coordinate  setup
for j=1:J
    x_eve(j)=45;
%     y_eve(1)=45;
    y_eve(j)=-5;
    z_eve(j)=0;
end
% x_eve(1)=1;
% x_eve(2)=1.5;
% x_eve(3)=1.25;
[G1,G2,h_iu1,h_iu2,Q,g1,g2]=double_Channel_setup1(M1,M2,K,N,J,x_ap,y_ap,...
    z_ap,x_irs1,y_irs1,z_irs1,x_irs2,y_irs2,z_irs2,L0,C0,x_user,y_user,...
    z_user,x_eve,y_eve,z_eve);
%% optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=zeros(K,100);
% R_L=zeros(100,1);
R_e=zeros(K,100);
R_u=zeros(K,100);
W_L=zeros(N,N,K,100);
Z_L=zeros(N,N,K,100);
v1_L=zeros(M1,1,100);
v2_L=zeros(M2,1,100);
v1_wuHL=zeros(1,M1,100);
v2_wuHL=zeros(1,M2,100);
V1_L=zeros(M1+1,M1+1,100);
V2_L=zeros(M2+1,M2+1,100);

l=1;
R_L(l)=0;
% w=1/N*ones(N,K);
% z=1/N*ones(N,K);
% w=rand(N,K);
% z=rand(N,K);
% w=ones(N,K);
% z=ones(N,K);
%  for k=1:K
%         W(:,:,k)=w(:,k)*w(:,k)';
%         Z(:,:,k)=z(:,k)*z(:,k)';
%         W_L(:,:,k,l)=W(:,:,k);
%         Z_L(:,:,k,l)=Z(:,:,k);
%     end
          load('pqfile.mat','W_out','Z_out')
          W_L(:,:,k,l)=W_out(1:N,1:N);%最大为10
          Z_L(:,:,k,l)=Z_out(1:N,1:N);
    theta=2*pi*rand(1,500);%phase shift
    v1_wuHL(:,:,l)=exp(theta(:,1:M1)*1j);%initialization passive beamforming vector v.  Note that v_up(K,1)=[e^(jθ1) ,...,e^(jθN))]^T
    v2_wuHL(:,:,l)=exp(theta(:,1:M2)*1j); 
    v1_L(:,:,l)=v1_wuHL(:,:,l)';
    v2_L(:,:,l)=v2_wuHL(:,:,l)';
    V1_L(:,:,l)=[v1_L(:,:,l)',1]'*[v1_L(:,:,l)',1];
    V2_L(:,:,l)=[v2_L(:,:,l)',1]'*[v2_L(:,:,l)',1];
h=1;
    while 1 
    l
%% v1 optimization
    [v1_out]=double_iteration_v1_optimization1(J,repeat,l,M1,M2,N,K,G1,G2,...
        h_iu1,h_iu2,Q,g1,g2,V1_L,v1_L,v2_L,v1_wuHL,v2_wuHL,eta,W_L,Z_L);
    v1_L(:,:,l+1)=v1_out;
    v1_wuHL(:,:,l+1)=v1_out';
    V1_L(:,:,l+1)=[v1_L(:,:,l+1)',1]'*[v1_L(:,:,l+1)',1];
%% v2 optimization
    [v2_out]=double_iteration_v2_optimization1(J,repeat,l,M1,N,K,G1,G2,h_iu1,...
        h_iu2,Q,g1,g2,V2_L,v1_L,v2_L,v1_wuHL,v2_wuHL,eta,W_L,Z_L);
    v2_L(:,:,l+1)=v2_out;
    v2_wuHL(:,:,l+1)=v2_out';
    V2_L(:,:,l+1)=[v2_L(:,:,l+1)',1]'*[v2_L(:,:,l+1)',1];
%% beamforming and AN optimization    
    [W_out,Z_out]=double_iteration_Beamforming_AN_optimization1(J,l,M1,M2,N,...
        K,G1,G2,h_iu1,h_iu2,Q,g1,g2,v1_L,v2_L,v1_wuHL,v2_wuHL,P0,eta,W_L,Z_L);
    W_L(:,:,:,l+1)=W_out;
    Z_L(:,:,:,l+1)=Z_out;
 %% output
  u_k_L=zeros(K,N);
  M_L=zeros(N,N,K);
  for k=1:K
     u_k_L(k,:)=h_iu2(:,k)'*diag(v2_wuHL(:,:,l+1))*Q*diag(v1_wuHL(:,:,l+1))*G1+...
                    h_iu1(:,k)'*diag(v1_wuHL(:,:,l+1))*G1+...
                    h_iu2(:,k)'*diag(v2_wuHL(:,:,l+1))*G2;
     M_L(:,:,k)=u_k_L(k,:)'*u_k_L(k,:);
 end

  u_e_L=zeros(J,N);
  M_e_L=zeros(N,N,J);
  for j=1:J
  u_e_L(j,:)=g2(:,j)'*diag(v2_wuHL(:,:,l+1))*Q*diag(v1_wuHL(:,:,l+1))*G1+...
                 g1(:,j)'*diag(v1_wuHL(:,:,l+1))*G1+...
                 g2(:,j)'*diag(v2_wuHL(:,:,l+1))*G2;
  M_e_L(:,:,j)=u_e_L(j,:)'*u_e_L(j,:);
  end

M_W_L=zeros(1,K);
M_Z_L=zeros(1,K);
Me_W_L=zeros(1,K);
Me_Z_L=zeros(1,K);
D1=zeros(1,K);
D2=zeros(J,K);
D3=zeros(1,K);
for k=1:K
for j=1:J
    for i=1:K
        M_W_L(i)=trace(M_L(:,:,k)*W_L(:,:,i,l+1));
        M_Z_L(i)=trace(M_L(:,:,k)*Z_L(:,:,i,l+1));
        Me_W_L(i)=trace(M_e_L(:,:,j)*W_L(:,:,i,l+1));
        Me_Z_L(i)=trace(M_e_L(:,:,j)*Z_L(:,:,i,l+1));
     end
    D1(k)=log2(sum(M_W_L)+sum(M_Z_L)+1)-log2(sum(M_W_L)+sum(M_Z_L)+1-trace(M_L(:,:,k)*W_L(:,:,k,l+1))); 
    D2(j,k)=log2(sum(Me_W_L)+sum(Me_Z_L)+1)-log2(sum(Me_W_L)+sum(Me_Z_L)+1-trace(M_e_L(:,:,j)*W_L(:,:,k,l+1)));
end
D3(k)=D1(k)-max(D2(:,k));
end
D1
D2
D3

% R_L(l+1)=real(D1-max(D2));
R_u_plot(l+1)=real(sum(D1))
R_e_plot(l+1)=real(sum(max(D2)))
R_L_sum(l+1)=real(sum(D3))
R_L_min(l+1)=min(real(D3))
h=h+1;
if R_L_sum(l+1)-R_L_sum(l)<=1e-4
%     R_L_sum(l+1)=R_L_sum(l)
% end
% if h>=30
break
end
l=l+1;

end
 %% plot
figure(1)
x=[0:1:l];
plot(x,R_L_sum,'r-');
xlabel('number of iterations');
ylabel('sum secrecy rate');
figure(2)
x=[0:1:l];
plot(x,R_u_plot,'r-',x,R_e_plot,'b-')
xlabel('number of iterations');
% ylabel('min rate');