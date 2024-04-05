clc %https://zhuanlan.zhihu.com/p/312069817
clear all
close all
x = 2:4;
y = rand(length(x),4);
M=10
% vector = [1, zeros(1, M-1)];
C1=linspace(1,1,M)-[1, zeros(1, M-1)];

bar(x,y,'stacked')