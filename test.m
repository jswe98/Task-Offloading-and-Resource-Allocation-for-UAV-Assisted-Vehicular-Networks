% ������ʼλ�ú��ٶ�
lleight=1;%��������������
lleft=-1;%��������������
initial_x = [0, 0, 0; 10, 0, 0; 20, 0, 0; 30, 0, 0]; % �������ĳ�ʼx����
velocity = 40; % �������ٶ�

% ģ����ʻ����
time = 0:0.1:10; % ʱ����Ϊ0.1�룬�ܹ�ģ��10�����ʻ����
num_cars = size(initial_x, 1); % ��������
Vposition(:,1) = [0; 0; 0];% �������˻���ʼλ��

for t = time
    clf; % ���֮ǰ��ͼ��
    
    % ���³���λ��
    positions = initial_x + velocity * t;
    
    % ���Ƴ���
    hold on;
    for i = 1:num_cars
        plot3([positions(i,1), positions(i,1)], [positions(i,2), positions(i,2)], [positions(i,3), positions(i,3)+0.1], 'b'); % ��ʾ�������߶Σ�ͨ����z����������һ����С�ĸ߶�����ʾ����
        text(positions(i,1), positions(i,2), positions(i,3)+0.1, num2str(i)); % �ڳ����Ϸ���ӱ��
    end
    hold off;
    
    % ����������
    xlim([0, max(positions(:,1))+10]);
    ylim([-10, 10]);
    zlim([-1, 1]);
    
    % ��ӱ�ǩ�ͱ���
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(['Time: ', num2str(t), ' seconds']);
    
    % ����ˢ��Ƶ��
    pause(0.01);
end

random_number = randi([1, 2]); % ����1��2֮������������1����200��2����400��

% ���������ѡ����ʾ����
if random_number == 1
    disp(200);
else
    disp(400);
end

random_number = randi([0, 1]); % ����0��1֮������������0����-1��1����1��
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
% ��ʾ���
disp(result);
fffs=[1 1]