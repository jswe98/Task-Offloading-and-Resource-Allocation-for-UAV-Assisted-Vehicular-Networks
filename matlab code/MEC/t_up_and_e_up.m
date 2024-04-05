
function [t_up,e_up] = t_up_and_e_up(d,p)
d_up=linspace(d,d,5);%上传的数据量
t_up=d_up./(rate(V1,p));%上传的时间
e_up=p.*t_up;%上传的能耗

