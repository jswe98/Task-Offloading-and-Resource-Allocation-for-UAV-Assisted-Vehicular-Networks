
function [t_up,e_up] = t_up_and_e_up(d,p)
d_up=linspace(d,d,5);%�ϴ���������
t_up=d_up./(rate(V1,p));%�ϴ���ʱ��
e_up=p.*t_up;%�ϴ����ܺ�

