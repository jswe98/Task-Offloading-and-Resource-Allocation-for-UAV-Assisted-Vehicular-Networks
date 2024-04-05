function result = binarySearchForMin(fx)
left=0;
right=6;
epsilon=1e-6;
f=fx;
% ����
    while right - left > epsilon
        mid1 = (2 * left + right) / 3;
        mid2 = (left + 2 * right) / 3;

        if f(mid1) > f(mid2)
            right = mid2;
        else
            left = mid1;
        end
    end

    % �������ұ߽��ƽ��ֵ��Ϊ���ֵ�Ĺ���
    max_val = (left + right) / 2;
    result=(f(max_val));
%     result1=result
%     disp(['���������ֵ����Ϊ��', num2str(f(max_val)), '����Ӧ�� x Ϊ��', num2str(max_val)]);
end