%% part 1:parameters initialization
clc
close all
clear
WR=[0,0];%coordinate of RIS
WG=[0,70];%coordinate of ground user
zG=150e3;
zU=80;%hight of UAV
zR=40;%hight of RIS
vmax=25;
d=0.5;%anntenna separation
lamda=1;%carrier wavelength
det=1;%solt length
D=vmax*det;%maximum distance
M=90;%number of RIS elements
sigma=db2pow(-80)*1e-3;%AWGN
P=1;%transmit power
gamma0=P/sigma;
epsilon=1e-3;%convergence accuracy
beita=db2pow(3);%Rician factor
alpha=2.8;
k=3.5;%pass loss exponent
rou=db2pow(-20);%pass loss at reference distance=1;
N=100;% number of time slots
T=N*det;%time span
Q0=zeros(N,2);% trajectory of UAV
Q0(1,:)=[-500,20];%initial horizontal location
Q0(N,:)=[500,20];%final horizontal location
for q=2:N-1   %initialization trajectory
    Q0(q,:)=[Q0(1,1)+q*abs(Q0(N,1)-Q0(1,1))/N,20];
end
b=0.158;
mm=19.4;
omega=1.29;
%% part 2: channel genelization
dUR=zeros(1,N);%distance between UAV-RIS
dUG=zeros(1,N);%distance between UAV-Ground user
hUR=zeros(N,M);%channel between UAV-RIS
hUG=zeros(1,N);%channel between UAV-Ground user
h=ShadowedRicianRVGenerator(1,b,mm,omega);
for n=1:N
    dUR(n)=sqrt((zU-zR)^2+norm(WR-Q0(n,:))^2);
    dUG(n)=sqrt((zU-zG)^2+norm(Q0(n,:)-WG)^2);
    hUR(n,:)=sqrt(rou*dUR(n)^(-2))*exp(-1i*2*pi*(0:M-1)*d*((WR(1,1)-Q0(n,1))/dUR(n))/lamda);
    hUG(n)=sqrt(rou*dUG(n)^(-k))*h;    
end
dRG=sqrt((zR-zG)^2+norm(WR-WG)^2);
hRG=sqrt(rou*dRG^(-alpha))*ShadowedRicianRVGenerator(M,b,mm,omega).';
A=sqrt(rou)*norm(h);
B=sqrt(rou)*sum(abs(hRG));
%% part 3: optimization
u0=dUG+10;
v0=dUR+10;
A0=zeros(1,N);
B0=zeros(1,N);
C0=zeros(1,N);
phi=zeros(N,M);
for n=1:N
    for m=1:M
        phi(n,m)=exp(1i*(angle(h)+angle(hRG(m))-angle(hUR(n,m))));
    end
    A0(n)=1+gamma0*(A^2/u0(n)^k+B^2/v0(n)^2+2*A*B/(u0(n)^(k/2)*v0(n)));
    B0(n)=-gamma0*(k*A^2/(u0(n)^(k+1))+k*A*B/(v0(n)*u0(n)^(k/2+1)));
    C0(n)=-gamma0*(2*B^2/v0(n)^3+2*A*B/(u0(n)^(k/2)*v0(n)^2));
end
ite=0;%initializtion of iteration number for AO
obj=zeros(1,100);%record objective value in each iteration
while 1
    %-----------------step1:optimize trajectory----------------------------
        cvx_begin
        cvx_solver sedumi
        variable Q(N,2)
        variable u(1,N)
        variable v(1,N)
        expression Rate_N(1,N)
        for n=1:N
            Rate_N(n)=B0(n)*u(n)/A0(n)/log(2)+C0(n)*v(n)/A0(n)/log(2);
        end
        maximize 1e10*sum(Rate_N)/N
        subject to
        for n=1:N
            (zU-zG)^2+sum_square(Q(n,:)-WG)+u0(n)^2-2*u0(n)*u(n)<=0;
            (zU-zR)^2+sum_square(Q(n,:)-WR)+v0(n)^2-2*v0(n)*v(n)<=0;
        end
        for n=1:N-1
            sum_square(Q(n+1,:)-Q(n,:))<=D^2;
        end
        sum_square(Q(N,:)-[500,20])<=D^2;
        Q(1,:)==[-500,20];
        cvx_end
        %------------------step 1.2:update---------------------------------
        Q0=Q;
        u0=u;
        v0=v;
        for n=1:N
            A0(n)=1+gamma0*(A^2/u0(n)^k+B^2/v0(n)^2+2*A*B/(u0(n)^(k/2)*v0(n)));
            B0(n)=-gamma0*(k*A^2/(u0(n)^(k+1))+k*A*B/(v0(n)*u0(n)^(k/2+1)));
            C0(n)=-gamma0*(2*B^2/v0(n)^3+2*A*B/(u0(n)^(k/2)*v0(n)^2));
        end
    %-----------------step2:optimize RIS phase shifts----------------------
    for n=1:N
        dUR(n)=sqrt((zU-zR)^2+norm(WR-Q0(n,:))^2);
        hUR(n,:)=sqrt(rou*dUR(n)^(-2))*exp(-1i*2*pi*(0:M-1)*d*((WR(1,1)-Q0(n,1))/dUR(n))/lamda);
    end
    for n=1:N
        for m=1:M
            phi(n,m)=exp(1i*(angle(h)+angle(hRG(m))-angle(hUR(n,m))));
        end
    end
    %----------------step3:judge convergence-------------------------------
    ite=ite+1;    
    for n=1:N
        obj(ite)=obj(ite)+(1/N)*log2(1+gamma0*norm(hUG(n)+hUR(n,:)*diag(phi(n,:))*hRG')^2);
    end
    if ite>1 && (abs(obj(ite)-obj(ite-1))/abs(obj(ite-1))<=epsilon || ite>=20)
        break
    end
end
hold on;axis on;grid on
% plot(WR(1),WR(2),'sk','MarkerSize',10);
% plot(WG(1),WG(2),'Ob','MarkerSize',10);
plot(Q0(1:N,1),Q0(1:N,2),'-^','LineWidth',1.2,'MarkerSize',5);













