function [t_exe,e_exe] = t_exe_and_e_exe(c,f)
c_exe=linspace(c,c,5);%��վ�����������
t_exe=c_exe./f;%��վ�����ʱ��
k=1;%ϵ��
e_exe=k*c_exe.*f.*f;%��վ�����ܺ�