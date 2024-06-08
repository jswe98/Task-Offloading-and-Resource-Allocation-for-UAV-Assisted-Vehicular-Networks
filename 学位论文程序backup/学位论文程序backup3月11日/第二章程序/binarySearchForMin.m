function result = binarySearchForMin(fx)
left=0;
right=6;
epsilon=1e-6;
f=fx;
% 精度
    while right - left > epsilon
        mid1 = (2 * left + right) / 3;
        mid2 = (left + 2 * right) / 3;

        if f(mid1) > f(mid2)
            right = mid2;
        else
            left = mid1;
        end
    end

    % 返回左右边界的平均值作为最大值的估计
    max_val = (left + right) / 2;
    result=(f(max_val));
%     result1=result
%     disp(['函数的最大值估计为：', num2str(f(max_val)), '，对应的 x 为：', num2str(max_val)]);
end