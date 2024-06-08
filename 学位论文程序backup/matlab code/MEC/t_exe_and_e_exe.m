function [t_exe,e_exe] = t_exe_and_e_exe(c,f)
c_exe=linspace(c,c,5);%基站处理的数据量
t_exe=c_exe./f;%基站处理的时间
k=1;%系数
e_exe=k*c_exe.*f.*f;%基站处理能耗