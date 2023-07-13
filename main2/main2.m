clc
close all
clear
T=60;%共有70个时隙
M=4;%车的数量
det=1;%solt length   TTT=T*det;%time span
roadwidth=8;
roadlength=800;
H=80;%无人机高度
VUAV_max=10;%无人机速度约束
%使用信噪比进行判断
X=zeros(M,T);
Delta=1e-6;%背景噪声是
C1=linspace(1,1,M);
%求解吞吐量第一个时隙的
P=zeros(T,M);
R=zeros(T,M);
[CARposition,l,GVR]=CARinfo(T,M,det,roadlength,roadwidth);
[UAVposition] = InitializationUAVposition(T,roadlength,roadwidth);

