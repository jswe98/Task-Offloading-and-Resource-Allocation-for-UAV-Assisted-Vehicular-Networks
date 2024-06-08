function result = lyapunovOptimization(lyapunov)

    % ������һ����͹����
%     flyapunov = @(x) x^2;

    % ���õ�������
    f=lyapunov;
    max_iter = 100;
    step_size = 0.01;

    % ��ʼ��
    x = 0;

    % ����
    for iter = 1:max_iter
        % �����ݶ�
        grad = 2 * x;

        % ���±���
        x = x - step_size * grad;
    end

    % ���ս��
    result = f(x);

    % ��ӡ���
%     disp(['��СֵΪ��', num2str(result)]);
end
