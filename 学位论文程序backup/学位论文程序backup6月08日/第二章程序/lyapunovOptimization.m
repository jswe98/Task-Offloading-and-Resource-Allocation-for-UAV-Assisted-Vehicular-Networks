function result = lyapunovOptimization(lyapunov)

    % 假设有一个非凸函数
%     flyapunov = @(x) x^2;

    % 设置迭代参数
    f=lyapunov;
    max_iter = 100;
    step_size = 0.01;

    % 初始化
    x = 0;

    % 迭代
    for iter = 1:max_iter
        % 计算梯度
        grad = 2 * x;

        % 更新变量
        x = x - step_size * grad;
    end

    % 最终结果
    result = f(x);

    % 打印结果
%     disp(['最小值为：', num2str(result)]);
end
