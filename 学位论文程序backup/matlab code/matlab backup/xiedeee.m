clc
close all
clear
%-------------------ϵͳ����-------------------%�о�����ĳ��ʱ���5��������ͨ�ţ�ÿ��������һ������CM�����ͷCH���ճ�����ͨ��
%-------------------ϵͳ����-------------------%���ͳ�CM�ٶ�32��30��34��34��30.��Ӧ��CH���ճ��ٶ�30��34��35��36��33
  L=zeros(5,5);
for F=1:6
V1=28*ones(1,5)+2*F*ones(1,5);%���ͳ�CM�ٶ�
% V1=[0,0,0,0,0];
V2=[30,30,30,30,30];%���ճ�CH�ٶ�
v=zeros(5,5);%���복֮�������ٶȾ���
h=zeros(5,5);%���복֮�����ǰʱ�̵�ֵ
g=zeros(5,5);%�����Ӳ���Ϊ1ָ���ֲ�
G=zeros(5,5);%���복֮���ŵ��������
I=zeros(5,5);%���복֮������ŵ�������󣨶Խ���ȫΪ0Ԫ�أ���Ϊ�Խ�������Ч�ŵ����棩
for i=1:5
    for j=1:5
        v(i,j)=abs(V1(i)-V2(j));
        h(i,j)=1+0.3*rand(1);
    end
end
% [6,1,4,0,3;0,7,2,6,9;8,1,6,2,1;10,3,8,4,1;0,7,2,6,9];����ٶȾ���
%h=[0.845418353066450,0.995224514302217,0.687960311535798,0.741632431164121,0.706496751675251;0.716599213586070,0.675355816604331,0.609143849450247,0.779761672779394,0.697456077357881;0.947490617895953,0.811444306977048,0.965654067256766,0.989572070823563,0.834170383464932;0.647590153069472,0.970613035572538,0.837424308992463,0.953446121969907,0.769790372435062;0.842902954353616,0.628305432990921,0.969908972473279,0.856831766456680,0.641799867159030];
h=[1.21026765186132,1.00186773622345,1.11230368038275,1.27044868578780,1.09550344842318;1.17912487696444,1.08933855474170,1.03750431438572,1.11650670059712,1.24530632319216;1.29435271127300,1.25859695058615,1.02514626462880,1.10131359710114,1.07083874991164;1.09534161637965,1.29533452164231,1.16447526737693,1.22477537780633,1.25255552041723;1.05006691209611,1.27092928278688,1.03153723834667,1.22352792893222,1.21881158135137];
Delta=1e-6;%��������
T=0.001;%���г��복�ƶ��ŵ�CSI�������ڣ�V2V��
c=3*1e+8;%����
fc=5.9*1e+9;%������Ƶ������Ƶ��
j1=2*pi*fc*T/c*v;%0�ױ���������������V2V��
epsi=besselj(0,j1);%����������ֵ��V2V��
a=(epsi.^2).*h;%����ͨ���ƶ��ŵ���˥���ģ���е�CSI�������֣���ǰʱ�̵��ŵ�״̬��Ϣ��
l=0.0001;%����ͨ����˥��ϵ��
G1=l*a;%����ͨ��˥�䣨��ǰʱ�̵��ŵ����棨Ҳ�ǲ�õ�ƽ���ŵ����棩��
G2=l*(1-epsi.^2);%����ͨ���ŵ�˥�����ֵ�ƽ���ŵ�����
G=G1+G2;%����ͨ���ܵ�ƽ���ŵ�������󣨸þ���������ȡƽ������ȣ��������ʺͣ�
  for i=1:5
      G(i,i)=100*G(i,i);%��Ч�ŵ�����ı������󣨴���ͨ�ţ�
  end 
MUk=0.5;
alfak=l*(1-epsi.^2);
betak=l*(1-epsi.^2);
X=G1+MUk*alfak+betak;%%xi����
SIGMA=0.5;
E=0.1;
H=X+SIGMA*sqrt(-2*log(E))*alfak; 
% ----------���챴��˹̹���Ƹ���Լ����ϵ������----------% 
I=H;%����˹̹����Լ������
    for i=1:5
      I(i,i)=0;%�������Լ���е�ϵ������
    end
    
global N
N =25;
arfa=zeros(N,5);
beta=zeros(N,5);
SINR=zeros(N,5);
Pmax=log(0.5)*ones(N,5);
B=1;
C1=linspace(1,1,5);
GF=linspace(1.5,1.5,5);%����ϵ��
PC=0.05;%��·�й̶����ĵĹ���
S=zeros(N,5);
s=zeros(N,5);
Gama=0.1;%ԭ����
Ith=2e-6;%������ֵGamaҪ����Ӧ�Ĵ洢�ռ�
Temp=zeros(40,5);%ÿ�ε������������Ź����ݴ��� 40ӦΪZ
Z=1;

while 1
P=zeros(N,5);
P(1,:)=[-8,-8,-8,-8,-8]; 
Lamda=zeros(N,5);
Lamda(1,:)=linspace(20,20,5);
for k=1:N
    for i=1:5
    I1(i)=exp(P(k,:))*G(:,i)+Delta-G(i,i)*exp(P(k,i));
    SINR(k,i)=G(i,i)*exp(P(k,i))/I1(i);     
    arfa(k,i)=SINR(k,i)/(1+SINR(k,i));
    Total_T(k)=B*dot(log(C1+SINR(k,:))/log(2),C1);
    end
    for i=1:5
        for j=1:5
            if j~=i
               S(k,j)=B/log(2)*arfa(k,j)*G(i,j)*SINR(k,j)/((exp(P(k,j))*G(j,j)))+Lamda(k,j)*I(i,j);%%������G(i,j)����G(j,i),Ӧ����G(i,j)����j�����ǵ�j�У�������Ǹ���CM��ĳ����ͷ�������Ϣ��
            else
               S(k,j)=0;
            end
         end
          s(k,i)=sum(S(k,:))+Gama(Z)*GF(i);
    end
%  Total_T(k+1)=B*dot(log(C1+SINR(k,:))/log(2),C1); 
P(k+1,:)=min(log(B.*arfa(k,:)/log(2))-log(s(k,:)),Pmax(k,:));
     if max(abs(exp(P(k+1,:))-exp(P(k,:))))<=1e-5
        P(k+1,:)=P(k,:);
     end
     for i=1:5
         Lamda(k+1,i)=max(0,Lamda(k,i)+1000000/k*(exp(P(k,:))*I(:,i)-Ith));
         w(k,i)=exp(P(k,:))*I(:,i)-Ith;
     end 
     if  max(abs(w(k,:)))<=3e-7    %����ֵ�ɸģ���ʾ������һ���ĳ̶��¾���Ϊ�Ѿ��ȶ���
         Lamda(k+1,:)=Lamda(k,:);
     end
end
D(Z)=Total_T(N);
Temp(Z,:)=P(N,:);
error =abs(D(Z)-Gama(Z)*(GF*exp(Temp(Z,:)')+PC));
Gama(Z+1)=D(Z)/(GF*exp(Temp(Z,:)')+PC);
      if  error<=1e-3
          L(F,:)=exp(Temp(Z,:));
          Q(F)=D(Z);
          A(F)=Gama(Z+1);
break;
%       else
%        Lamda=zeros(N,5);
%     Lamda(1,:)=linspace(4,6,5);
%         P=zeros(N,5);
%         P =Temp(Z,:); 
      end 
          Z=Z+1;
end
end
   
%{
plot(10*log(exp(P(:,1)))/log(10),'-+r');
hold on
plot(10*log(exp(P(:,2)))/log(10),'-og');
plot(10*log(exp(P(:,3)))/log(10),'-*b');
plot(10*log(exp(P(:,4)))/log(10),'-sk');
plot(10*log(exp(P(:,5)))/log(10),'-dc');
legend('P1','P2','P3','P4','P5')
xlabel('iteration');
ylabel('Power conversion value(dB)');
%}
  Q1=Q;
  A1=A;
plot(exp(P(:,1)),'-+r');
hold on
plot(exp(P(:,2)),'-og');
plot(exp(P(:,3)),'-*b');
plot(exp(P(:,4)),'-sk');
plot(exp(P(:,5)),'-dc');

legend('P1','P2','P3','P4','P5')
xlabel('iteration');
ylabel('Power conversion value(W)');

figure
plot(Lamda(:,1),'-+r');
hold on
plot(Lamda(:,2),'-og');
plot(Lamda(:,3),'-*b');
plot(Lamda(:,4),'-sk');
plot(Lamda(:,5),'-dc');
legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5')
xlabel('iteration');
ylabel('Multiplier convergence');

figure
plot(Gama,'-sb','linewidth',1.5) 
hold on
plot(D,'-^r','linewidth',1.5)

hold on

[AX,H1,H2] =plotyy(Z,Gama,Z,D,@plot);% ��ȡ�����ᡢͼ����
set(AX(1),'XColor','k','YColor','b','LineWidth',1.5);
set(AX(2),'XColor','k','YColor','r','LineWidth',1.5)

set(get(AX(1),'ylabel'),'color','blue','string', 'Average EE(bits/Joule/Hertz)','fontsize',12);
set(get(AX(2),'ylabel'),'color','red','string', 'Average SE(bits/sec/Hertz)','fontsize',12);

xlabel('iteration','fontsize',12);
legend('EE','SE')


    FF=[0,2,4,6,8,10];
    figure
 [AX,H1,H2] = plotyy(FF,A1,FF,Q1,'plot');

set(AX(1),'XColor','k','YColor','b','LineWidth',1.5);
set(AX(2),'XColor','k','YColor','r','LineWidth',1.5);

HH1=get(AX(1),'Ylabel');
set(HH1,'String','Average EE(bits/Hz/Joule)','fontsize',12);
set(HH1,'color','b');

HH2=get(AX(2),'Ylabel');
set(HH2,'String','Average SE(bits/s/Hz)','fontsize',12);
set(HH2,'color','r');

set(H1,'marker','s')
set(H1,'LineStyle','-','LineWidth',1.5);
set(H1,'color','b');

set(H2,'LineStyle','-','LineWidth',1.5);
set(H2,'color','r');
set(H2,'marker','^')
% hold on 
% EE=[1.53634734393827,1.53579347534572,1.53541361465141,1.53510426222248,1.53482959098909];
% SE=[0.231773787999855,0.231780191199803,0.231784581226184,0.231788155496835,0.231791328385694];

legend([H1,H2],{'EE';'SE'});
xlabel('\Delta\upsilon (m/s)','fontsize',12);
legend('EE,C=C^*','SE,C=C^*')


% set(H1,'Linestyle','--');
% 
% set(H2,'Linestyle',':');

% set(gcf,'color','white')
% set(gca,'xcolor','r');
% set(gca,'ycolor','r','linewidth',1.5) %���ñ߿���  
