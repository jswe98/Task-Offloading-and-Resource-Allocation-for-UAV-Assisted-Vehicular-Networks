clc
close all
clear

cvx_begin
c = [1;-1];
A = [-10,7;1,1/2];
b = [-5;3];
d=length(c);
        variable x_Recon(d);
        
        minimize(c'*x_Recon)
        subject to
            A*x_Recon<=b;
            x_Recon>=0;
    cvx_end
display(x_Recon);


clc
close all
clear
m = 16; n = 8;
A = randn(m,n);
b = randn(m,1);

cvx_begin
    variable x(n)
    minimize( norm(A*x-b) )
cvx_end

clc
close all
clear
m = 16; n = 8;
A = randn(m,n);
b = randn(m,1);
bnds = randn(n,2);
l = min( bnds, [], 2 );
u = max( bnds, [], 2 );
cvx_begin
    variable x(n)
    minimize( norm(A*x-b) )
    subject to
        l <= x <= u
cvx_end
