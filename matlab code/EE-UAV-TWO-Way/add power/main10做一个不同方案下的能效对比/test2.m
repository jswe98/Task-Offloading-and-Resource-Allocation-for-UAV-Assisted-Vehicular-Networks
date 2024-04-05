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
    % 在这里使用 q 的值进行其他操作  
    disp(qq);  
end




T=90;
% 定义无人机的固定坐标  
fix_x = 5;  % x坐标  
fix_y = 10; % y坐标  
valuefix = zeros(T, 2);  % 创建一个与时间点数量相同的2列向量   % 模拟无人机在每个时间点的位置  
for i = 1:T  
    valuefix(i,1) = fix_x;  % 设置x坐标  
    valuefix(i,2) = fix_y;  % 设置y坐标   
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
    % 在这里使用 q 的值进行其他操作  
    disp(qq);  
end
  
% 显示位置向量  
disp(valuefix);