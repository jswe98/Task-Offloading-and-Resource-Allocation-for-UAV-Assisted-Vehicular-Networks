for programmatic = 1:3  
    switch programmatic  
        case 1  
            qq = 'value1';  
        case 2  
            qq = 'value2';  
        case 3  
            qq =valuefix;  
        otherwise  
            disp('Invalid value of programmatic');  
    end  
    % ������ʹ�� q ��ֵ������������  
    disp(qq);  
end




T=90;
% �������˻��Ĺ̶�����  
fix_x = 5;  % x����  
fix_y = 10; % y����  
valuefix = zeros(T, 2);  % ����һ����ʱ���������ͬ��2������   % ģ�����˻���ÿ��ʱ����λ��  
for i = 1:T  
    valuefix(i,1) = fix_x;  % ����x����  
    valuefix(i,2) = fix_y;  % ����y����   
end  
for programmatic = 1:3  
    switch programmatic  
        case 1  
            qq = 'value1';  
        case 2  
            qq = 'value2';  
        case 3  
            qq =valuefix;  
        otherwise  
            disp('Invalid value of programmatic');  
    end  
    % ������ʹ�� q ��ֵ������������  
    disp(qq);  
end
  
% ��ʾλ������  
disp(valuefix);