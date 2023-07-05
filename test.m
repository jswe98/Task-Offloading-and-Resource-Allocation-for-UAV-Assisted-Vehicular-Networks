% 车辆初始位置和速度
lleight=1;%车道索引车向右
lleft=-1;%车道索引车向左
initial_x = [0, 0, 0; 10, 0, 0; 20, 0, 0; 30, 0, 0]; % 四辆车的初始x坐标
velocity = 40; % 车辆的速度

% 模拟行驶过程
time = 0:0.1:10; % 时间间隔为0.1秒，总共模拟10秒的行驶过程
num_cars = size(initial_x, 1); % 车辆数量
Vposition(:,1) = [0; 0; 0];% 定义无人机初始位置

for t = time
    clf; % 清除之前的图像
    
    % 更新车辆位置
    positions = initial_x + velocity * t;
    
    % 绘制车辆
    hold on;
    for i = 1:num_cars
        plot3([positions(i,1), positions(i,1)], [positions(i,2), positions(i,2)], [positions(i,3), positions(i,3)+0.1], 'b'); % 表示车辆的线段，通过在z方向上增加一个很小的高度来表示车道
        text(positions(i,1), positions(i,2), positions(i,3)+0.1, num2str(i)); % 在车辆上方添加编号
    end
    hold off;
    
    % 设置坐标轴
    xlim([0, max(positions(:,1))+10]);
    ylim([-10, 10]);
    zlim([-1, 1]);
    
    % 添加标签和标题
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(['Time: ', num2str(t), ' seconds']);
    
    % 控制刷新频率
    pause(0.01);
end

random_number = randi([1, 2]); % 生成1到2之间的随机整数（1代表200，2代表400）

% 根据随机数选择并显示数字
if random_number == 1
    disp(200);
else
    disp(400);
end

random_number = randi([0, 1]); % 生成0到1之间的随机整数（0代表-1，1代表1）
if random_number == 0
    ll = lleight;
    random_number = randi([0, 1]);
    if random_number == 0
        Vposition_x=0;
    else
        Vposition_x=200;
    end        
else
    ll = lleft;
    if random_number == 0
        Vposition_x=400;
    else
        Vposition_x=800;
    end 
end
fff=UAVposition(1:2, 88) 
% 显示结果
disp(result);
fffs=[1 1]