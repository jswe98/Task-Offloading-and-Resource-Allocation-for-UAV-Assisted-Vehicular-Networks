function [no_p] = no_p(~)
clc
close all
clear
V1=[20,32,22,38,30];%���ͳ����ٶ�
V2=[0,0,0,0,0];%���ջ�վ���ٶ�
v=zeros(5,5);%�����վ֮�������ٶȾ���
h=zeros(5,5);%�����վ֮�����ǰʱ�̵�ֵ
for i=1:5
    for j=1:5
        v(i,j)=abs(V1(i)-V2(j));
        h(i,j)=1+0.3*rand(1);
    end
end
h=[1.19329543905811,1.11358281479808,1.24347413748474,1.15984767663984,1.10521813107307;1.28170046859997,1.26278284344790,1.16504690286953,1.18674252580037,1.17611341135943;1.06232268781991,1.09037389908385,1.14127700455528,1.06914644806347,1.25329263780862;1.05842928687011,1.06777653429172,1.05121241414436,1.06829928934497,1.13070960523117;1.09333068599512,1.27701389263097,1.12906221739888,1.05544489603724,1.27146429060397];
Delta=1e-6;%��������
f5 = 5.8; f4 = f5; f3 = f4; f2 = f3; f1 = f2;
t_imax=0.25; %ִ��ʱ�����������̵�ʱ��
Tc=0.001;%transmission latency between RSU to Cloud is set to fix
f_ba=5;%��Ե�Ĺ̶�����
d_up=0.5;%�ϴ���������
c_exe=2;%ִ�е�������
tau=100;%���ݰ�����
f_total=30;%�Ʒָ���Ե�ļ�����Դ
varsigma=200;%���ݰ�����ǿ��
Dmax=0.2;%ʱ����ֵ
e2=0.1;%�ڶ�����ֵ
B=1;%����
T=0.0005;%��վ�ɼ���֮ͨ�ŵķ��䳵�����ŵ���·��CSI����(V2I)
c=3*1e+8;%����
fc=5.9*1e+9;%������Ƶ������Ƶ��
j1=2*pi*fc*T/c*v;%0�ױ���������������V2I��
epsi=besselj(0,j1);%����������ֵ��V2I��
a=(epsi.^2).*h;%�����վͨ���ƶ��ŵ���˥���ģ���е�CSI�������֣���ǰʱ�̵��ŵ�״̬��Ϣ��
l=0.0001;%�����վͨ����˥��ϵ��
G1=l*a;%V2Iͨ��˥�䣨��ǰʱ�̵��ŵ����棨Ҳ�ǲ�õ�ƽ���ŵ����棩��
G2=l*(1-epsi.^2);%v2iͨ���ŵ�˥�����ֵ�ƽ���ŵ�����
G=G1+G2;%�����վͨ���ܵ�ƽ���ŵ�������󣨸þ���������ȡƽ������ȣ��������ʺͣ�
  for i=1:5
      G(i,i)=100*G(i,i);%��Ч�ŵ�����ı������󣨴���ͨ�ţ�
  end 
MUk=0.5;
alfak=l*(1-epsi.^2);
betak=l*(1-epsi.^2);
diagG=-diag(G1)/9;%SINR��ֵ��Ϊ9
G11=zeros(5,5);
 for i=1:5
     G11(i,i)=diagG(i,1);
 end
 for i=1:5
      G1(i,i)=0;
 end
 G1=G11+G1;
X=G1+MUk*alfak+betak;%xi����
SIGMA=0.5;
E=0.1;%��һ����ֵ
H=X+SIGMA*sqrt(-2*log(E))*alfak; 
I=H;%����˹̹����Լ������
global N
N=40;%��������
arfa=zeros(N,5);
SINR=zeros(N,5);
Pmax=-9*ones(N,5);
C1=linspace(1,1,5);
S=zeros(N,5);
s=zeros(N,5);
Ith=2e-6;%������ֵGamaҪ����Ӧ�Ĵ洢�ռ�
Temp=zeros(N,5);%ÿ�ε������������Ź����ݴ��� 40ӦΪZ
Z=1;%��������ÿ�μ�һ
I1=zeros(N,5);
Total_T=zeros(N,1);
P=zeros(N,5);
Lamda=zeros(N,5);
Miu=zeros(N,5);
Ksai=zeros(N,5);
Fai=zeros(N,5);
fi=zeros(N,5);
sumfai=zeros(N,5);
sumfi=zeros(N,5);
u_iexe=zeros(N,5);
P=-8*ones(N,5);
Lamda(1,:)=linspace(15,15,5);
Miu(1,:)=linspace(20,20,5);

    for i=1:5
    I1(i)=exp(P(1,:))*G(:,i)+Delta-G(i,i)*exp(P(1,i));
    SINR(1,i)=G(i,i)*exp(P(1,i))/I1(i);      
    end
for u=1:N
    u_iexe(u,:)=max([0,0,0,0,0],(t_imax*C1-(Tc+(c_exe./(f_ba*C1+fi(u,:)))))/t_imax*(1/d_up/log(2)));%����ִ��ʱ���Ч�ö����һ���ϴ�ǰ׺
    fi(1,:)=[f1,f2,f3,f4,f5];
    Ksai(1,:)=linspace(3000,3000,5);
    Fai(1,:)=linspace(7,7,5);
    for i=1:5
    I1(i)=exp(P(1,:))*G(:,i)+Delta-G(i,i)*exp(P(1,i));
    SINR(1,i)=G(i,i)*exp(P(1,i))/I1(i);      
    arfa(1,i)=SINR(1,i)/(1+SINR(1,i));
    Total_T(u)=dot(log(C1+SINR(1,:)),u_iexe(u,:)); %���ǵ�һ��������ĺ͵���ʽ
    end
    for i=1:5
        RRR=(exp(P(1,:))*G(i,i))/(I1(i)+Delta);
        Ri=(log(C1+RRR))/log(2);
        D1=C1./(tau*Ri-varsigma);
        BR=Ri.*c_exe/(t_imax*d_up); 
        BRR=sum(BR)+c_exe.*Ksai(u,:);
        sumfai(u,i)=sum(Fai(u,:));
        xiajie(u,i)=c_exe/(t_imax-Tc)-f_ba;
        fi(u+1,:)=max(xiajie(u,i),sqrt(BRR/sumfai(u,i))-f_ba);
    end    
   f1= fi(u+1,1);
   f2= fi(u+1,2);
   f3= fi(u+1,3);
   f4= fi(u+1,4);
   f5= fi(u+1,5);
   sumfi(u,1)=sum(fi(u,:));       
    for i=1:5       
        Ksai(u+1,i)=Ksai(u,i)+1000/u*(c_exe/(fi(u,i)+f_ba)+D1(1,i)-Dmax);
    end    
    for i=1:5       
        sumfi(u,i)=sum(fi(u,:));
        Fai(u+1,i)=max(0,Fai(u,i)+1/u*(sumfi(u,i)-f_total));
    end  
end     
R=zeros(N,1);
U=zeros(N,1);
no_p=Total_T(:,1);
%     for i=1:5
%         R=(log(1+(exp(P(:,i))*G(j,j))/(I1(i)+Delta)))/log(2);  
%         U=1-(Tc+c_exe./(fi(:,i)+f_ba))/t_imax;
%     end
%     R(N+1,1)=R(N,1);
%     Total_T=R.*U/d_up;
% figure
% plot((fi(:,1)),'-+r');
% hold on
% plot((fi(:,2)),'-og');
% plot((fi(:,3)),'-*b');
% plot((fi(:,4)),'-sk');
% plot((fi(:,5)),'-dc');
% legend('F1','F2','F3','F4','F5')
% xlabel('iteration');
% % set(gca,'ylim',[3,10]); 
% set(gca,'xlim',[20,50]);
% ylabel('computing');
% 
% figure
% plot((Total_T(:,1)),'-+r');
% legend('F1')
% xlabel('iteration');
% ylabel('EE');