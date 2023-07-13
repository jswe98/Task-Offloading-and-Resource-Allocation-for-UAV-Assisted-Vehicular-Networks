function [UAVposition] = InitializationUAVposition(T,roadlength,roadwidth)
Qv=zeros(T,2);%UAV初始化坐标
Qv(1,:)=[0,0];
Qv(T,:)=[roadlength,roadwidth];
for n=2:T-1
    Qv(n,:)=Qv(n-1,:)+[roadlength/(T-1),roadwidth/(T-1)];
end

UAVposition=Qv;
% 绘制轨迹图
% figure;
% plot(positions(:, 1), positions(:, 2), 'b-', 'LineWidth', 2);
% hold on;
% plot([initial_position(1), target_position(1)], [initial_position(2), target_position(2)], 'r--', 'LineWidth', 2);
% xlabel('X (m)');
% ylabel('Y (m)');
% legend('无人机轨迹', '终点');
% grid on;