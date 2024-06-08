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

syms r ct real;%声明rct为实数型
x4=r*exp(ct*i);%复指数形式复数的符号表达
subs(x4,{r,ct},{sqrt(2),3*pi/4})%代入具体值
v = [2 1 -1 -2 -5]; %对角线元素
SEITA=
% [H_RX, H_TX] = channelGen(2, 2, 1, 1, 1, 1);
% 
% A=H_RX;
% B=H_TX