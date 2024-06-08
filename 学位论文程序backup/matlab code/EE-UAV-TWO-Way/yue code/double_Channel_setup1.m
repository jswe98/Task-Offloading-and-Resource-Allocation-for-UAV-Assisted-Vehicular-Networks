% function [G1,G2,h_iu1,h_iu2,Q,g1,g2]=double_Channel_setup1(M1,M2,K,N,J,x_ap,y_ap,z_ap,x_irs1,y_irs1,z_irs1,...
%     x_irs2,y_irs2,z_irs2,C0,x_user,y_user,z_user,x_eve,y_eve,z_eve)
% %computing AP-IRS distance
% d_ai1=sqrt((x_ap-x_irs1)^2+(y_ap-y_irs1)^2+(z_ap-z_irs1)^2);
% d_ai2=sqrt((x_ap-x_irs2)^2+(y_ap-y_irs2)^2+(z_ap-z_irs2)^2);
% %computing IRS-USER distance
% d_iu1=zeros(M1,K);
% d_iu2=zeros(M2,K);
% for k=1:K
%     d_iu1(:,k)=sqrt((x_irs1-x_user(k))^2+(y_irs1-y_user(k))^2+(z_irs1-z_user(k))^2);
%     d_iu2(:,k)=sqrt((x_irs2-x_user(k))^2+(y_irs2-y_user(k))^2+(z_irs2-z_user(k))^2);
% end
% %computing IRS1-IRS2 distance
% d_irs=sqrt((x_irs1-x_irs2)^2+(y_irs1-y_irs2)^2+(z_irs1-z_irs2)^2);
% %computing IRS-Eve distance
% d_ie1=zeros(M1,J);
% d_ie2=zeros(M2,J);
% for j=1:J
%     d_ie1(:,j)=sqrt((x_irs1-x_eve(j))^2+(y_irs1-y_eve(j))^2+(z_irs1-z_eve(j))^2);
%     d_ie2(:,j)=sqrt((x_irs2-x_eve(j))^2+(y_irs2-y_eve(j))^2+(z_irs2-z_eve(j))^2);
% end
% 
% %% channel setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C0 = db2pow(-30); % 参考距离时的路损
% D0=1; % 参考距离
% sigmaK2 = db2pow(-80); % 噪声功率
% L = @(d, alpha)C0*(d/D0)^(-alpha); % 路损模型
% 
% % 路损参数
% alpha_G1= 2;
% alpha_G2= 2;
% alpha_Q= 2;
% alpha_Iu1 = 2.8;
% alpha_Iu2 = 2.8;
% alpha_Ie1 = 3.5; 
% alpha_Ie2 = 3.5; 
% beta_IU = 0; % IRS到User考虑瑞利衰落信道，AP和IRS之间为纯LoS信道
% 
% G1 = zeros(M1,N);
% G1 = sqrt(L(d_ai1,alpha_G1)/sigmaK2)*ones(M1,N);% LoS信道,除以噪声功率是为了进行噪声功率归一化，因为G和hr是级联的，对一个信道进行归一化即可
% % G1=sqrt(10^(-L0/10)/pow_pos(d_ai1,2));
% % temp11=(rand(M1,120)+1j*randn(M1,120));
% % G1=G1.*(sqrt(1/2)*temp11(:,1:N));
% 
% G2=zeros(M2,N);
% G2 = sqrt(L(d_ai2,alpha_G2)/sigmaK2)*ones(M2,N);
% % G2(:,:)=sqrt(10^(-L0/10)/pow_pos(d_ai2,2));
% % temp12=(rand(M2,120)+1j*randn(M2,120));
% % G2=G2.*(sqrt(1/2)*temp12(:,1:N));
% 
% %IRS-USER loss exponents=2.8
% h_iu1=zeros(M1,K);
% for k=1:K
%     h_iu1(:,k) = sqrt(L(d_iu1(k),alpha_Iu1)/2)*(randn(M1,1)+1i*randn(M1,1)); % 瑞利信道, IRS-User
% end
% % for k=1:K
% %     h_iu1(:,k)=sqrt(10^(-L0/10)/pow_pos(d_iu1(k),2));
% % end
% % temp21=(rand(500,K)+1j*randn(500,K));
% % h_iu1=h_iu1.*(sqrt(1/2)*temp21(1:M1,:));
% 
% h_iu2=zeros(M2,K);
% for k=1:K
%     h_iu2(:,k) = sqrt(L(d_iu2(k),alpha_Iu2)/2)*(randn(M2,1)+1i*randn(M2,1)); % 瑞利信道, IRS-User
% end
% % for k=1:K
% %     h_iu2(:,k)=sqrt(10^(-L0/10)/pow_pos(d_iu2(k),2));
% % end
% % temp22=(rand(500,K)+1j*randn(500,K));
% % h_iu2=h_iu2.*(sqrt(1/2)*temp22(1:M2,:));
% 
% % %IRS-Eve loss exponents=4
% g1=zeros(M1,J);
% for j=1:J
%     g1(:,j) = sqrt(L(d_ie1(j),alpha_Ie1)/2)*(randn(M1,1)+1i*randn(M1,1)); % 瑞利信道, IRS-User
% end
% % for j=1:J
% %     g1(:,j)=sqrt(10^(-L0/10)/pow_pos(d_ie1(j),4));
% % end
% % g1=g1.*(sqrt(1/2)*(rand(M1,J)+1j*randn(M1,J)));
% 
% g2=zeros(M2,J);
% for j=1:J
%     g2(:,j) = sqrt(L(d_ie2(j),alpha_Ie2)/2)*(randn(M2,1)+1i*randn(M2,1));
% end
% % for j=1:J
% %     g2(:,j)=sqrt(10^(-L0/10)/pow_pos(d_ie2(j),4));
% % end
% % g2=g2.*(sqrt(1/2)*(rand(M2,J)+1j*randn(M2,J)));
% 
% % %IRS1-IRS2 loss exponents=2
% Q=zeros(M2,M1);
% Q = sqrt(L(d_irs,alpha_Q))*ones(M2,M1); 
% % Q(:,:)=sqrt(10^(-L0/10)/pow_pos(d_irs,2));
% % Q=Q.*(sqrt(1/2)*(rand(M2,M1)+1j*randn(M2,M1)));

% end

function [G1,G2,h_iu1,h_iu2,Q,g1,g2]=double_Channel_setup1(M1,M2,K,N,J,x_ap,y_ap,z_ap,x_irs1,y_irs1,z_irs1,...
    x_irs2,y_irs2,z_irs2,L0,C0,x_user,y_user,z_user,x_eve,y_eve,z_eve)
sigmaK2 = db2pow(-80); % 噪声功率
%computing AP-IRS distance
d_ai1=sqrt((x_ap-x_irs1)^2+(y_ap-y_irs1)^2+(z_ap-z_irs1)^2);
d_ai2=sqrt((x_ap-x_irs2)^2+(y_ap-y_irs2)^2+(z_ap-z_irs2)^2);
%computing IRS-USER distance
d_iu1=zeros(M1,K);
d_iu2=zeros(M2,K);
for k=1:K
    d_iu1(:,k)=sqrt((x_irs1-x_user(k))^2+(y_irs1-y_user(k))^2+(z_irs1-z_user(k))^2);
    d_iu2(:,k)=sqrt((x_irs2-x_user(k))^2+(y_irs2-y_user(k))^2+(z_irs2-z_user(k))^2);
end
%computing IRS1-IRS2 distance
d_irs=sqrt((x_irs1-x_irs2)^2+(y_irs1-y_irs2)^2+(z_irs1-z_irs2)^2);
%computing IRS-Eve distance
d_ie1=zeros(M1,J);
d_ie2=zeros(M2,J);
for j=1:J
    d_ie1(:,j)=sqrt((x_irs1-x_eve(j))^2+(y_irs1-y_eve(j))^2+(z_irs1-z_eve(j))^2);
    d_ie2(:,j)=sqrt((x_irs2-x_eve(j))^2+(y_irs2-y_eve(j))^2+(z_irs2-z_eve(j))^2);
end
% d_ie1(:)=sqrt((x_irs1-x_eve)^2+(y_irs1-y_eve)^2+(z_irs1-z_eve)^2);
% 
% d_ie2(:)=sqrt((x_irs2-x_eve)^2+(y_irs2-y_eve)^2+(z_irs2-z_eve)^2);

%% channel setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G1=zeros(M1,N);
G1(:,:)=sqrt(10^(-L0/10)/pow_pos(d_ai1,2));
temp11=(rand(M1,120)+1j*randn(M1,120));
G1=G1.*(sqrt(1/2)*temp11(:,1:N))/sqrt(sigmaK2);

% /sqrt(sigmaK2);

G2=zeros(M2,N);
G2(:,:)=sqrt(10^(-L0/10)/pow_pos(d_ai2,2));
temp12=(rand(M2,120)+1j*randn(M2,120));
G2=G2.*(sqrt(1/2)*temp12(:,1:N))/sqrt(sigmaK2);

% /sqrt(sigmaK2);

%IRS-USER loss exponents=2
h_iu1=zeros(M1,K);
for k=1:K
    h_iu1(:,k)=sqrt(10^(-L0/10)/pow_pos(d_iu1(k),2));
end
temp21=(rand(500,K)+1j*randn(500,K));
h_iu1=h_iu1.*(sqrt(1/2)*temp21(1:M1,:));
%h_iu1(:,2)=0; 

h_iu2=zeros(M2,K);
for k=1:K
    h_iu2(:,k)=sqrt(10^(-L0/10)/pow_pos(d_iu2(k),2));
end
temp22=(rand(500,K)+1j*randn(500,K));
h_iu2=h_iu2.*(sqrt(1/2)*temp22(1:M2,:));
%h_iu2(:,2)=0;


% %IRS-Eve loss exponents=2
g1=zeros(M1,J);
for j=1:J
    g1(:,j)=sqrt(10^(-L0/10)/pow_pos(d_ie1(j),2));
end
g1=g1.*(sqrt(1/2)*(rand(M1,J)+1j*randn(M1,J)));

g2=zeros(M2,J);
for j=1:J
    g2(:,j)=sqrt(10^(-L0/10)/pow_pos(d_ie2(j),2));
end
g2=g2.*(sqrt(1/2)*(rand(M2,J)+1j*randn(M2,J)));

% %IRS1-IRS2 loss exponents=2
Q=zeros(M2,M1);
Q(:,:)=sqrt(10^(-L0/10)/pow_pos(d_irs,2));
Q=Q.*(sqrt(1/2)*(rand(M2,M1)+1j*randn(M2,M1)));

% %% IRS1-IRS2 LOS
% beta0=10^(-5);%单位信道增益beta0=10^(-5)的结果是0.3左右，-4的结果是30左右
% alpha=2.8;%IRS与用户之间的路损
% beta=3.5;%UAV与用户之间的路损
% %IRS1到IRS2距离
% D_12=norm([x_irs1 y_irs1 z_irs1]-[x_irs2 y_irs2 z_irs2]);
% %IRS1到IRS2路损
% pathloss_12=sqrt(beta0/square_pos(D_12));%beta0
% %IRS1到IRS2天线相应
% array_12=zeros(M2,M1);%待
% %arrive
% phi_12ar=pi/2;
% theta_12ar=pi/2;
% for m= 0:sqrt(M2)-1
%     for n= 0:sqrt(M2)-1
%         array_12ar(m*(sqrt(M2))+n+1) = exp( 1i* pi* ( m*sin(phi_12ar)*sin(theta_12ar) + n*cos(theta_12ar) ) );
%     end
% end
% array_12ar = array_12ar.'/sqrt(M2);
% phi_12at=pi/2;
% theta_12at=pi/2;
% for m= 0:sqrt(M1)-1
%     for n= 0:sqrt(M1)-1
%         array_12at(m*(sqrt(M1))+n+1) = exp( 1i* pi* ( m*sin(phi_12at)*sin(theta_12at) + n*cos(theta_12at) ) );
%     end
% end
% array_12at = array_12at.'/sqrt(M1);
% array_12=array_12ar*array_12at;
% Q=pathloss_12*array_12;
% 
end

