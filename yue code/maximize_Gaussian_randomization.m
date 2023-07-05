function [nu_opt] = maximize_Gaussian_randomization(K,V,l,F,W_L,repeat)
[U, S, ~] = svd(V);
len = size(U, 1);
obj_opt = -inf;
obj=0;
for round = 1 : repeat
    r = sqrt(1/2) * ( randn(len,1) + 1j*randn(len,1));
    nu_bar = U*sqrt(S)*r;
    nu = exp( 1j * ...                 
        angle( nu_bar(1:len-1)/nu_bar(len) ) ...
        );
    for k=1:K
        obj=obj+real(trace(F(:,:,k)'*W_L(:,:,k,l)*F(:,:,k)*[nu; 1]*[nu; 1]'));
    end
     if obj > obj_opt
        obj_opt = obj; 
        nu_opt = nu;
     end

end
end

