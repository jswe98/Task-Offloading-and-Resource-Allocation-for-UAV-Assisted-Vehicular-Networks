% ����һЩ����  
x = 0:0.1:10;  
y1 = sin(x);  
y2 = cos(x);  
  
% ����ͼ��  
plot(x, y1, 'r-', x, y2, 'b--')  
  
% ���ͼ��������ϣ����ĸtheta  
legend(['\theta=0.2', '\theta cos(x)'], 'Interpreter', 'latex', 'Location', 'best')