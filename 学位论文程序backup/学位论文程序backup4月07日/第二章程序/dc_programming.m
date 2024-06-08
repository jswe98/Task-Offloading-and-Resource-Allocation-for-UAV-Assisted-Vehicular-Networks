function result = dc_programming(g,f)

    % ����������͹����
%     f = @(x) x^2;
%     g = @(x) exp(x);

    % ���õ�������
    max_iter = 100;
    step_size = 0.01;

    % ��ʼ��
    x = 0;

    % ����
    for iter = 1:max_iter
        % ������ݶ�
        subgrad = 2*x - exp(x);

        % ���±���
        x = x - step_size * subgrad;
    end

    % ���ս��
    result = f(x) - g(x);

    % ��ӡ���
%     disp(['��СֵΪ��', num2str(result)]);
end
