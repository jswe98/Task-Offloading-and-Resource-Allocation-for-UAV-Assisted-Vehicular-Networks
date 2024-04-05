function [v1_out] = double_iteration_v1_optimization1(J,repeat,l,M1,M2,N,K,G1,G2,...
    h_iu1,h_iu2,Q,g1,g2,V1_L,v1_L,v2_L,v1_wuHL,v2_wuHL,eta,W_L,Z_L)
% repeat=5000;%Gaussian operation times

[M1, ~,~] = size(v1_L);
% T=3;

A=zeros(M1,N,K);
b=zeros(K,N);
Q_bar=zeros(M2,M1,K);
F=zeros(N,M1+1,K);
temp1=zeros(M1,N,M2);
for k=1:K
Q_bar(:,:,k)=diag(h_iu2(:,k)')*Q;
for m=1:M2
    tempt1(:,:,m)=conj(v2_L(m,1,l))*diag(Q_bar(m,:,k))*G1;
end
A(:,:,k)=sum(tempt1,3)+diag(h_iu1(:,k)')*G1;
b(k,:)=v2_L(:,:,l)'*diag(h_iu2(:,k)')*G2;
F(:,:,k)=[A(:,:,k)',b(k,:)'];
end

A_e=zeros(M1,N,J);
b_e=zeros(J,N);
Q_ebar=zeros(M2,M1,J);
F_e=zeros(N,M1+1,J);
tempt2=zeros(M1,N,M2);
for j=1:J
Q_ebar(:,:,j)=diag(g2(:,j)')*Q;
for m=1:M2
    tempt2(:,:,m)=conj(v2_L(m,1,l))*diag(Q_ebar(m,:,j))*G1;
end
A_e(:,:,j)=sum(tempt2,3)+diag(g1(:,j)')*G1;
b_e(j,:)=v2_L(:,:,l)'*diag(g2(:,j)')*G2;
F_e(:,:,j)=[A_e(:,:,j)',b_e(j,:)'];
% F_e=F_e/eta;
end

F_W=zeros(1,K);
F_Z=zeros(1,K);
Fe_W=zeros(1,K);
Fe_Z=zeros(1,K);

t=1;
v1_T=zeros(M1,1,100);
V1_T=zeros(M1+1,M1+1,100);
v1_T(:,:,t)=v1_L(:,:,l);
V1_T(:,:,t)=[v1_T(:,:,t)',1]'*[v1_T(:,:,t)',1];
Alpha_T=zeros(1,1);
Alpha_T(t)=0;
while 1
    t1=zeros(1,K);
    t2=zeros(J,K);
    for k=1:K
        for j=1:J
            for i=1:k
                F_W(:,i)=real(trace(F(:,:,k)'*W_L(:,:,i,l)*F(:,:,k)*V1_T(:,:,t)));   
                F_Z(:,i)=real(trace(F(:,:,k)'*Z_L(:,:,i,l)*F(:,:,k)*V1_T(:,:,t)));
                Fe_W(:,i)=real(trace(F_e(:,:,j)'*W_L(:,:,i,l)*F_e(:,:,j)*V1_T(:,:,t)));
                Fe_Z(:,i)=real(trace(F_e(:,:,j)'*Z_L(:,:,i,l)*F_e(:,:,j)*V1_T(:,:,t)));
            end
            t1(k)=1/(sum(F_W)+sum(F_Z)+eta-trace(F(:,:,k)'*W_L(:,:,k,l)*F(:,:,k)*V1_T(:,:,t)));
            t2(j,k)=1/(sum(Fe_W)+sum(Fe_Z)+eta);     
        end
    end

cvx_solver SDPT3
cvx_precision best
cvx_begin  quiet

    variable V1(M1+1,M1+1) complex hermitian
    variable t3(1,K)
    expression F_Wi(1,K)
    expression F_Zi(1,K)
    expression Fe_Wi(1,K)
    expression Fe_Zi(1,K)
    expression fai1(1,K)
    expression fai2(J,K)

    for k=1:K
        for j=1:J
            for i=1:K
                F_Wi(:,i)=real(trace(F(:,:,k)'*W_L(:,:,i,l)*F(:,:,k)*V1));
                F_Zi(:,i)=real(trace(F(:,:,k)'*Z_L(:,:,i,l)*F(:,:,k)*V1));
                Fe_Wi(:,i)=real(trace(F_e(:,:,j)'*W_L(:,:,i,l)*F_e(:,:,j)*V1));
                Fe_Zi(:,i)=real(trace(F_e(:,:,j)'*Z_L(:,:,i,l)*F_e(:,:,j)*V1));
             end
             fai1(k)=log(sum(F_Wi)+sum(F_Zi)+eta)-real(t1(k)*(sum(F_Wi)+sum(F_Zi)+eta-real(trace(F(:,:,k)'*W_L(:,:,k,l)*F(:,:,k)*V1))))+real(log(t1(k)))+1;
             fai2(j,k)=t2(j,k)*(sum(Fe_Wi)+sum(Fe_Zi)+eta)-log(sum(Fe_Wi)+sum(Fe_Zi)+eta-real(trace(F_e(:,:,j)'*W_L(:,:,k,l)*F_e(:,:,j)*V1)))-log(t2(j,k))-1;
        end
    end

    maximize sum(fai1-t3)
    subject to
    for k=1:K
        for j=1:J
            t3(k)>=fai2(j,k);
        end
%     fai1(k)>1;
    end
    t3>=0;
%     for k=1:K
%         fai1(k)>=0;
%     end
    for m=1:M1+1
        V1(m,m)==1;
    end
    V1==hermitian_semidefinite(M1+1);

cvx_end


V=V1;
v1=maximize_Gaussian_randomization(K,V,l,F,W_L,repeat);

V1_T(:,:,t+1)=V1;
v1_T(:,:,t+1)=v1;
Alpha_T(t+1)=cvx_optval

if  Alpha_T(t+1)-Alpha_T(t)<=1e-3
%     abs(Alpha_T(t+1)-Alpha_T(t))<=1e-1
break
end
t=t+1;
end

if Alpha_T(t+1)-Alpha_T(t)>=0
    v1_out=v1_T(:,:,t+1);
    else
    v1_out=v1_T(:,:,t); 
end

% v1_out=v1_T(:,:,t+1);
end

