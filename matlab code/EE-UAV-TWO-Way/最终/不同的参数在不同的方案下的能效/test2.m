clc
close all
clear
% �ܳ�ʼ��
import truncnorm
T=70; M=4;
Vvelocity=zeros(1,M);
for i=1:M
    Vvelocity(i)=randi([10, 20]);
end
% random_number = randi([10, 20]);
% Vvelocity = 10;  % �����ٶ�Ϊ40
% ���뷶Χ�;�ֵ����׼��
vmin = 60;  % ��С�ٶ�
vmax = 120; % ����ٶ�
mu = 90;    % ��ֵ
sigma = 10; % ��׼��

% ����ضϸ�˹�ֲ��Ĳ���
a = (vmin - mu) / sigma; % �½ضϲ���
b = (vmax - mu) / sigma; % �Ͻضϲ���

% ���ɽضϸ�˹�ֲ������
n = 1000; % ���ɵ����������
X_trunc = truncnormrnd(a, b, mu, sigma, [n,1]);

% ����ֱ��ͼ
histogram(X_trunc, 'Normalization', 'pdf')
xlabel('�����ٶ�')
ylabel('�����ܶ�')
title('�ضϸ�˹�ֲ�')

% ��ʾ��ֵ�ͷ���
mean(X_trunc)
var(X_trunc)