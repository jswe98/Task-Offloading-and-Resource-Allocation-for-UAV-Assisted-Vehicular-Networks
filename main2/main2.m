clc
close all
clear
T=60;%����70��ʱ϶
M=4;%��������
det=1;%solt length   TTT=T*det;%time span
roadwidth=8;
roadlength=800;
H=80;%���˻��߶�
VUAV_max=10;%���˻��ٶ�Լ��
%ʹ������Ƚ����ж�
X=zeros(M,T);
Delta=1e-6;%����������
C1=linspace(1,1,M);
%�����������һ��ʱ϶��
P=zeros(T,M);
R=zeros(T,M);
[CARposition,l,GVR]=CARinfo(T,M,det,roadlength,roadwidth);
[UAVposition] = InitializationUAVposition(T,roadlength,roadwidth);

