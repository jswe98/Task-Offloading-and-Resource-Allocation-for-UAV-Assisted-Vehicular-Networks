function result = dc_programming(g,f)

    % 假设有两个凸函数
%     f = @(x) x^2;
%     g = @(x) exp(x);

    % 设置迭代参数
    max_iter = 100;
    step_size = 0.01;

    % 初始化
    x = 0;

    % 迭代
    for iter = 1:max_iter
        % 计算次梯度
        subgrad = 2*x - exp(x);

        % 更新变量
        x = x - step_size * subgrad;
    end

    % 最终结果
    result = f(x) - g(x);

    % 打印结果
%     disp(['最小值为：', num2str(result)]);
end
