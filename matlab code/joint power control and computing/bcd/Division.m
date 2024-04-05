% 定义区间[2,3]
x = 2:0.1:3;
a = 2;
b = 3;
% 定义参考y值，即y=0的直线
ty = zeros(1,length(x));
% 定义在区间[2,3]上单调的函数fun
fun = @(x) 3*x.^3-9*x.^2+5.6*x-7.5;
y = fun(x);
% 绘制示意图
plot(x,y,'b.-',x,ty,'r--','LineWidth',3.5);
xlabel('x轴');
ylabel('y轴');
title('二分法寻根测试');
% 定义计算精度ep和临时变量tmp
ep = 1e-6;
tmp = 1;
% 对于不知道循环次数的，用whlie来实现
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
        disp('此区间不存在根！！！');
        break;
    end
end
JG = fun(mid);
disp(['在区间[2,3]的近似解为',num2str(JG)]);
hold on;
% 绘制近似交点
plot(mid,JG,'r.','MarkerSize',50);
hold off;


clc;
fun=@(x) exp(x)+10*x-2;
Division(fun,0.005,0,1);