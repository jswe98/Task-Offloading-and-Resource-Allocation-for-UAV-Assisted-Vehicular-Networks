clc
close all
clear
V1=[30,32,28,34,30];%���ͳ����ٶ�
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
Delta=1e-6;%��������
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
MUk=0.5;
alfak=l*(1-epsi.^2);
betak=l*(1-epsi.^2);
X=G1+MUk*alfak+betak;%xi����
SIGMA=0.5;
E=0.1;
H=X+SIGMA*sqrt(-2*log(E))*alfak; 
I=H;%����˹̹����Լ������
for i=1:5
    I(i,i)=0;
end%�������Խ���Ϊ0
N=30;%��������
arfa=zeros(N,5);
beta=zeros(N,5);
SINR=zeros(N,5);
Pmax=log(0.5)*ones(N,5);
C1=linspace(1,1,5);
S=zeros(N,5);
s=zeros(N,5);
Gama=0.1;%ԭ����
Ith=2e-6;%������ֵGamaҪ����Ӧ�Ĵ洢�ռ�
Temp=zeros(N,5);%ÿ�ε������������Ź����ݴ��� 40ӦΪZ
Z=1;%��������ÿ�μ�һ
for sub1=1:N
P=zeros(N,5);
Lamda=zeros(N,5);
Miu=zeros(N,5);
P(1,:)=[-8,-8,-8,-8,-8]; 
Lamda(1,:)=linspace(20,20,5);
Miu(1,:)=linspace(20,20,5);


end%sub1 end