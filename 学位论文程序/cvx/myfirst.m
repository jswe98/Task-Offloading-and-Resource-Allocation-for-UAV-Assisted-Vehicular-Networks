cvx_begin
variable z;
variable x;
variable y;
minimize (z)
subject to
log(1+x) >= 2;
sqrt(x) >= 3;
y==8./x;
norm(2*z, (x-y)) >= x+y;
cvx_end
