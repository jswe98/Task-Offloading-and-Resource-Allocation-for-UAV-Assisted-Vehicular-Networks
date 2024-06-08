% 创建一些数据  
x = 0:0.1:10;  
y1 = sin(x);  
y2 = cos(x);  
  
% 创建图形  
plot(x, y1, 'r-', x, y2, 'b--')  
  
% 添加图例并插入希腊字母theta  
legend(['\theta=0.2', '\theta cos(x)'], 'Interpreter', 'latex', 'Location', 'best')