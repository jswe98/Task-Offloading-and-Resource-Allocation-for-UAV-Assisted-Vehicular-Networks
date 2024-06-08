% INPUTS:
%   H       : number of reflecting elements along the vertical direction
%   W       : number of reflecting elements along the horizontal direction
%   NantRX  : number of antennas at the RX
%   NantTX  : number of antennas at the TX
%   NrayRX  : number of sub-paths in the IRS-RX channel
%   NrayTX  : number of sub-paths in the TX-IRS channel
%
% OUTPUTS:
%   H_RX    : IRS-RX channel matrix
%   H_TX    : TX-IRS channel matrix

syms r ct real;%����rctΪʵ����
x4=r*exp(ct*i);%��ָ����ʽ�����ķ��ű��
subs(x4,{r,ct},{sqrt(2),3*pi/4})%�������ֵ
v = [2 1 -1 -2 -5]; %�Խ���Ԫ��
SEITA=
% [H_RX, H_TX] = channelGen(2, 2, 1, 1, 1, 1);
% 
% A=H_RX;
% B=H_TX