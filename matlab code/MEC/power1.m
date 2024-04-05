clc
close all
clear
V1=[24,32,28,34,30];%���ͳ����ٶ�
V2=[0,0,0,0,0];%���ջ�վ���ٶ�
v=zeros(5,5);%�����վ֮�������ٶȾ���
h=zeros(5,5);%�����վ֮�����ǰʱ�̵�ֵ
for i=1:5
    for j=1:5
        v(i,j)=abs(V1(i)-V2(j));
        h(i,j)=1+0.3*rand(1);
    end
end
%%h=[1.21026765186132,1.00186773622345,1.11230368038275,1.27044868578780,1.09550344842318;1.17912487696444,1.08933855474170,1.03750431438572,1.11650670059712,1.24530632319216;1.29435271127300,1.25859695058615,1.02514626462880,1.10131359710114,1.07083874991164;1.09534161637965,1.29533452164231,1.16447526737693,1.22477537780633,1.25255552041723;1.05006691209611,1.27092928278688,1.03153723834667,1.22352792893222,1.21881158135137];
%%%%h=[1.02789667806120,1.13904677432867,1.00279975360833,1.27450777352074,1.19282252173993;1.00042571746072,1.00911558174848,1.06254106699160,1.13648984350063,1.03817980999754;1.00259430636343,1.21812388034456,1.10623493981698,1.23413378382881,1.13099699249059;1.13096643471168,1.01476395410695,1.01488957045582,1.02733005262392,1.17821110943324;1.07232521655071,1.25241073059164,1.25716382922718,1.28908366028541,1.14666993584818];
Delta=1e-6;%��������
f5 = 6; f4 = f5; f3 = f4; f2 = f3; f1 = f2;
t_imax=0.3; %ִ��ʱ�����������̵�ʱ��
Tc=0.01;%transmission latency between RSU to Cloud is set to fix
f_ba=5;%��Ե�Ĺ̶�����
d_up=0.5;%�ϴ���������
c_exe=2;%ִ�е�������
tau=100;%���ݰ�����
f_total=25;%�Ʒָ���Ե�ļ�����Դ
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
      G(i,i)=10*G(i,i);%��Ч�ŵ�����ı������󣨴���ͨ�ţ�
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
N=90;%��������
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
Total_T=zeros(N,5);
P=zeros(N,5);
Lamda=zeros(N,5);
Miu=zeros(N,5);
Ksai=zeros(N,5);
Fai=zeros(N,5);
fi=zeros(N,5);
sumfai=zeros(N,5);
sumfi=zeros(N,5);
u_iexe=zeros(N,5);
P(1,:)=[-8,-8,-8,-8,-8]; 
Lamda(1,:)=linspace(15,15,5);
Miu(1,:)=linspace(20,20,5);
for bcd=1:90
for k=1:N
    fik(k,:)=[f1,f2,f3,f4,f5];%����ÿ���Ƹ���Ե���������
    u_iexe(k,:)=(t_imax*C1-(Tc+(c_exe./(f_ba*C1+fik(k,:)))))/t_imax*(1/d_up/log(2));%����ִ��ʱ���Ч�ö����һ���ϴ�ǰ׺
    for i=1:5
    I1(i)=exp(P(k,:))*G(:,i)+Delta-G(i,i)*exp(P(k,i));
    SINR(k,i)=G(i,i)*exp(P(k,i))/I1(i);      
    arfa(k,i)=SINR(k,i)/(1+SINR(k,i));
    Total_T(k)=dot(log(C1+SINR(k,:)),u_iexe(k,:)); %���ǵ�һ��������ĺ͵���ʽ
    end
    for i=1:5
        for j=1:5
            if j~=i
               S(k,j)=u_iexe(k,j)*arfa(k,j)*G(i,j)*SINR(k,j)/((exp(P(k,j))*G(j,j)))+Lamda(k,j)*I(i,j);%%������G(i,j)����G(j,i),Ӧ����G(i,j)����j�����ǵ�j�У��������Ǹ���CM��ĳ����ͷ�������Ϣ��
            else
               S(k,j)=-Miu(k,j)*(log(1-e2)-G(j,j));
            end
        end
               s(k,i)=sum(S(k,:));
    end
%     P(k+1,:)=min(log(u_iexe(k,:).*arfa(k,:)/log(2))-log(s(k,:)),Pmax(k,:));
    P(k+1,:)=log(u_iexe(k,:).*arfa(k,:)/log(2))-log(s(k,:));
    for i=1:5
         Lamda(k+1,i)=max(0,Lamda(k,i)+1000000/k*(exp(P(k,:))*G1(:,i)+Delta));
    end 
    for i=1:5
        A=(varsigma/tau+1/tau/(Dmax-c_exe/fi(k,i)));%z���ǹ̶�f��ĳ�����
        test1=power(2,A)*(Ith+Delta);
        Miu(k+1,i)=max(0,Miu(k,i)+100/k*(test1+exp(P(k,i))*(log(1-e2)-G(j,j))));%%������Ҫ�����󵼺�Ľ��
    end  
end
for u=1:N
    fi(1,:)=[f1,f2,f3,f4,f5];
    Ksai(1,:)=linspace(1,1,5);
    Fai(1,:)=linspace(15,15,5);
    for i=1:5
        RRR=(exp(P(k,:))*G(i,i))/(Ith+Delta);
        Ri=(log(C1+RRR))/log(2);
        D1=C1./(tau*Ri-varsigma);
        BR=Ri.*c_exe*c_exe/(t_imax*d_up); 
        BRR=BR.*Ksai(u,:);
        sumfai(u,i)=sum(Fai(u,:));
        fi(u+1,:)=max(0,sqrt(BRR)-f_ba);
    end    
   f1= fi(u+1,1);
   f2= fi(u+1,2);
   f3= fi(u+1,3);
   f4= fi(u+1,4);
   f5= fi(u+1,5);
    for i=1:5       
        Ksai(u+1,i)=Ksai(u,i)+10/u*(c_exe/(fi(u,i)+f_ba)+D1(1,i)-Dmax);
    end    
%     for i=1:5       
%         sumfi(u,i)=sum(fi(u,:));
%         Fai(u+1,i)=Fai(u,i)+1/u*(sumfi(u,i)-f_total);
%     end  
end    
end