cvx_begin 
variable z;
variable x;
variable y;
minimize (x)
subject to
y.*log(1+x) >= 2;
sqrt(x) >= 3;
% z*x-x-z==0;
cvx_end
