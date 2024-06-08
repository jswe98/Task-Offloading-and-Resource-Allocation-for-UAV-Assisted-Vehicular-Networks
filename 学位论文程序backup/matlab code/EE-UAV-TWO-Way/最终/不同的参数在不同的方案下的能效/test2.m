clc
close all
clear
% 总初始化
import truncnorm
T=70; M=4;
Vvelocity=zeros(1,M);
for i=1:M
    Vvelocity(i)=randi([10, 20]);
end
% random_number = randi([10, 20]);
% Vvelocity = 10;  % 车车速度为40
% 输入范围和均值、标准差
vmin = 60;  % 最小速度
vmax = 120; % 最大速度
mu = 90;    % 均值
sigma = 10; % 标准差

% 计算截断高斯分布的参数
a = (vmin - mu) / sigma; % 下截断参数
b = (vmax - mu) / sigma; % 上截断参数

% 生成截断高斯分布随机数
n = 1000; % 生成的随机数数量
X_trunc = truncnormrnd(a, b, mu, sigma, [n,1]);

% 绘制直方图
histogram(X_trunc, 'Normalization', 'pdf')
xlabel('车辆速度')
ylabel('概率密度')
title('截断高斯分布')

% 显示均值和方差
mean(X_trunc)
var(X_trunc)