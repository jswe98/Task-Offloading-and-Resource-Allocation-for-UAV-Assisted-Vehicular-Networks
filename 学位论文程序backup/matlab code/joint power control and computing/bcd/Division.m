% ��������[2,3]
x = 2:0.1:3;
a = 2;
b = 3;
% ����ο�yֵ����y=0��ֱ��
ty = zeros(1,length(x));
% ����������[2,3]�ϵ����ĺ���fun
fun = @(x) 3*x.^3-9*x.^2+5.6*x-7.5;
y = fun(x);
% ����ʾ��ͼ
plot(x,y,'b.-',x,ty,'r--','LineWidth',3.5);
xlabel('x��');
ylabel('y��');
title('���ַ�Ѱ������');
% ������㾫��ep����ʱ����tmp
ep = 1e-6;
tmp = 1;
% ���ڲ�֪��ѭ�������ģ���whlie��ʵ��
while(tmp>ep)
    if fun(a)*fun(b) < 0
       mid = (a+b)/2;
       if fun(a)*fun(mid) < 0
          b = mid;
       end
       if fun(b)*fun(mid) < 0
          a = mid;
       end
       tmp = abs(a-b);
    else
        disp('�����䲻���ڸ�������');
        break;
    end
end
JG = fun(mid);
disp(['������[2,3]�Ľ��ƽ�Ϊ',num2str(JG)]);
hold on;
% ���ƽ��ƽ���
plot(mid,JG,'r.','MarkerSize',50);
hold off;


clc;
fun=@(x) exp(x)+10*x-2;
Division(fun,0.005,0,1);