clc
close all
clear
T=70;%����70��ʱ϶
t=1;%��һ��ʱ϶
lleight=1;%��������������
lleft=-1;%��������������
%UAV coordinate setup
VUAV_max=10;%���˻��ٶ�Լ��
VUAU_x
VUAU_y
for t=1:T
QUAV(t)=[VUAV_max*t,VUAV_max*t,100];
end
%vehicle coordinate setup
for i=1:4
    x_user(i)=6;
    y_user(i)=2;
end
x_user(1)=6;
x_user(2)=8;
x_user(3)=8;
y_user(1)=2;
y_user(2)=1.5;
y_user(3)=2;