function [v2_out] = double_iteration_v2_optimization1(J,repeat,l,M1,N,K,G1,G2,...
    h_iu1,h_iu2,Q,g1,g2,V2_L,v1_L,v2_L,v1_wuHL,v2_wuHL,eta,W_L,Z_L)
% repeat=5000;%Gaussian operation times
[M2, ~,~] = size(v2_L);
% T=3;
A=zeros(M2,N,K);
b=zeros(K,N);
Q_bar=zeros(M2,M1,K);
F=zeros(N,M2+1,K);
for k=1:K
    Q_bar(:,:,k)=diag(h_iu2(:,k)')*Q;
    A(:,:,k)=Q_bar(:,:,k)*diag(v1_wuHL(:,:,l+1))*G1+diag(h_iu2(:,k)')*G2;
    b(k,:)=v1_L(:,:,l+1)'*diag(h_iu1(:,k)')*G1;
    F(:,:,k)=[A(:,:,k)',b(k,:)'];
%     F(:,:,k)=F(:,:,k)/eta;
end

A_e=zeros(M2,N,J);
b_e=zeros(J,N);
Q_ebar=zeros(M2,M1,J);
F_e=zeros(N,M2+1,J);
for j=1:J
    Q_ebar(:,:,j)=diag(g2(:,j)')*Q;
    A_e(:,:,j)=Q_ebar(:,:,j)*diag(v1_wuHL(:,:,l+1))*G1+diag(g2(:,j)')*G2;
    b_e(j,:)=v1_L(:,:,l+1)'*diag(g1(:,j)')*G1;
    F_e(:,:,j)=[A_e(:,:,j)',b_e(j,:)'];
% F_e=F_e/eta;
end

F_W=zeros(1,K);
F_Z=zeros(1,K);
Fe_W=zeros(1,K);
Fe_Z=zeros(1,K);

t=1;
v2_T=zeros(M2,1,100);
V2_T=zeros(M2+1,M2+1,100);
v2_T(:,:,t)=v2_L(:,:,l);
V2_T(:,:,t)=[v2_T(:,:,t)',1]'*[v2_T(:,:,t)',1];
Alpha_T=zeros(1,1);
Alpha_T(t)=0;
while 1
    t1=zeros(1,K);
    t2=zeros(J,K);
for k=1:K
    for j=1:J
        for i=1:k
            F_W(:,i)=real(trace(F(:,:,k)'*W_L(:,:,i,l)*F(:,:,k)*V2_T(:,:,t)));
            F_Z(:,i)=real(trace(F(:,:,k)'*Z_L(:,:,i,l)*F(:,:,k)*V2_T(:,:,t)));
            Fe_W(:,i)=real(trace(F_e(:,:,j)'*W_L(:,:,i,l)*F_e(:,:,j)*V2_T(:,:,t)));
            Fe_Z(:,i)=real(trace(F_e(:,:,j)'*Z_L(:,:,i,l)*F_e(:,:,j)*V2_T(:,:,t)));
        end
        t1(k)=1/(sum(F_W)+sum(F_Z)+eta-trace(F(:,:,k)'*W_L(:,:,k,l)*F(:,:,k)*V2_T(:,:,t)));
        t2(j,k)=1/(sum(Fe_W)+sum(Fe_Z)+eta);
    end
end

cvx_solver SDPT3
cvx_precision best
cvx_begin quiet

    variable V2(M2+1,M2+1) complex hermitian
    variable t3(1,K) 
    expression F_Wi(1,K)
    expression F_Zi(1,K)
    expression Fe_Wi(1,K)
    expression Fe_Zi(1,K)
    expression fai1(1,K)
    expression fai2(J,K)

for k=1:K
    for j=1:J
        for i=1:k
            F_Wi(:,i)=real(trace(F(:,:,k)'*W_L(:,:,i,l)*F(:,:,k)*V2));
            F_Zi(:,i)=real(trace(F(:,:,k)'*Z_L(:,:,i,l)*F(:,:,k)*V2));
            Fe_Wi(:,i)=real(trace(F_e(:,:,j)'*W_L(:,:,i,l)*F_e(:,:,j)*V2));
            Fe_Zi(:,i)=real(trace(F_e(:,:,j)'*Z_L(:,:,i,l)*F_e(:,:,j)*V2));
        end 
        fai1(k)=log(sum(F_Wi)+sum(F_Zi)+eta)-real(t1(k)*(sum(F_Wi)+sum(F_Zi)+eta-real(trace(F(:,:,k)'*W_L(:,:,k,l)*F(:,:,k)*V2))))+real(log(t1(k)))+1;
        fai2(j,k)=t2(j,k)*(sum(Fe_Wi)+sum(Fe_Zi)+eta)-log(sum(Fe_Wi)+sum(Fe_Zi)+eta-real(trace(F_e(:,:,j)'*W_L(:,:,k,l)*F_e(:,:,j)*V2)))-log(t2(j,k))-1;
    end
end

    maximize sum(fai1-t3)
    subject to

    for k=1:K
        for j=1:J
            t3(k)>=fai2(j,k);
        end
    end
    t3>=0;
%     for k=1:K
%         fai1(k)>=0;
%     end
    for m=1:M2+1
        V2(m,m)==1;
    end
    V2==hermitian_semidefinite(M2+1);

cvx_end

V=V2;
v2=maximize_Gaussian_randomization(K,V,l,F,W_L,repeat);

V2_T(:,:,t+1)=V2;
v2_T(:,:,t+1)=v2;
Alpha_T(t+1)=cvx_optval

if Alpha_T(t+1)-Alpha_T(t)<=1e-3
    %     abs(Alpha_T(t+1)-Alpha_T(t))<=1e-1
break
end
t=t+1;
end

if Alpha_T(t+1)-Alpha_T(t)>=0
   v2_out=v2_T(:,:,t+1);
else
   v2_out=v2_T(:,:,t); 
end

% v2_out=v2_T(:,:,t+1);
end

