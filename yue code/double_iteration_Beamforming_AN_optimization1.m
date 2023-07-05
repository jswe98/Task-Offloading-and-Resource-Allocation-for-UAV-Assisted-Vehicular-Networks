function [W_out,Z_out]=double_iteration_Beamforming_AN_optimization1(J,l,M1,M2,N,K,G1,G2,h_iu1,h_iu2,Q,g1,g2,...
    v1_L,v2_L,v1_wuHL,v2_wuHL,P0,eta,W_L,Z_L)
% T=3;
%user equivalent link
u_k=zeros(K,N);
M=zeros(N,N,K);
for k=1:K
    u_k(k,:)=h_iu2(:,k)'*diag(v2_wuHL(:,:,l+1))*Q*diag(v1_wuHL(:,:,l+1))*G1+...
        h_iu1(:,k)'*diag(v1_wuHL(:,:,l+1))*G1+...
        h_iu2(:,k)'*diag(v2_wuHL(:,:,l+1))*G2;
    M(:,:,k)=u_k(k,:)'*u_k(k,:);
end
% Eve equivalent link
u_e=zeros(J,N);
M_e=zeros(N,N,J);
for j=1:J
    u_e(j,:)=g1(:,j)'*diag(v1_wuHL(:,:,l+1))*G1+...
        g2(:,j)'*diag(v2_wuHL(:,:,l+1))*G2+...
        g2(:,j)'*diag(v2_wuHL(:,:,l+1))*Q*diag(v1_wuHL(:,:,l+1))*G1;
    M_e(:,:,j)=u_e(j,:)'*u_e(j,:);
end

M_W=zeros(1,K);
M_Z=zeros(1,K);
Me_W=zeros(1,K);
Me_Z=zeros(1,K);
t=1;
W_T=zeros(N,N,K,100);
W_T(:,:,:,t)=W_L(:,:,:,l);
Z_T=zeros(N,N,K,100);
Z_T(:,:,:,t)=W_L(:,:,:,l);
Alpha_T=zeros(1,1);
Alpha_T(t)=0;
% for t=1:T
while 1
    t1=zeros(1,K);
    t2=zeros(J,K);
    for k=1:K
        for i=1:K
            M_W(i)=real(trace(M(:,:,k)*W_T(:,:,i,t))); 
            M_Z(i)=real(trace(M(:,:,k)*Z_T(:,:,i,t)));
        end
        t1(k)=1/(sum(M_W)+sum(M_Z)+eta-real(trace(M(:,:,k)*W_T(:,:,k,t))));
    end
    for k=1:K
        for j=1:J
            for i=1:K
                Me_W(i)=real(trace(M_e(:,:,j)*W_T(:,:,i,t)));
                Me_Z(i)=real(trace(M_e(:,:,j)*Z_T(:,:,i,t))); 
            end
            t2(j,k)=1/(sum(Me_W)+sum(Me_Z)+eta);
        end
%         t2(j,k)=1/(sum(Me_W)+sum(Me_Z)+1);
    end

    cvx_solver SDPT3
    cvx_precision best
    cvx_begin  quiet

    variable W(N,N,K) complex hermitian
    variable Z(N,N,K) complex hermitian
    variable t3(1,K)
    expression M_Wi(1,K)
    expression M_Zi(1,K)
    expression Me_Wi(1,K) 
    expression Me_Zi(1,K)
    expression Trac(1,K)
%     expression Sumz(1,K)%xin
    expression fai1(1,K)
    expression fai2(J,K)

    for k=1:K
        for j=1:J
            for i=1:K
                M_Wi(i)=real(trace(M(:,:,k)*W(:,:,i)));
                M_Zi(i)=real(trace(M(:,:,k)*Z(:,:,i)));
                Me_Wi(i)=real(trace(M_e(:,:,j)*W(:,:,i)));
                Me_Zi(i)=real(trace(M_e(:,:,j)*Z(:,:,i)));
            end
            fai1(k)=log(sum(M_Wi)+sum(M_Zi)+eta)-t1(k)*(sum(M_Wi)+sum(M_Zi)+eta-trace(M(:,:,k)*W(:,:,k)))+log(t1(k))+1;
            fai2(j,k)=t2(j,k)*(sum(Me_Wi)+sum(Me_Zi)+eta)-log(sum(Me_Wi)+sum(Me_Zi)+eta-trace(M_e(:,:,j)*W(:,:,k)))-log(t2(j,k))-1;
        end
        Trac(k)=trace(W(:,:,k)+Z(:,:,k));
%         Sumz(k)=trace(Z(:,:,k));%xin

    end
        
    maximize sum(fai1-t3)
    subject to
    sum(Trac)<=P0; %AP功率约束
%     sum(Sumz)>=0.5; %xin
    for k=1:K
        for j=1:J
            t3(k)>=fai2(j,k);
        end
    end
%     fai1>zeros(1,K);
%     for k=1:K
%         fai1(k)>=0;
%     end
    t3>=0;
    for k=1:K %semidefinite constraint
        W(:,:,k)==hermitian_semidefinite(N);
        Z(:,:,k)==hermitian_semidefinite(N); 
    end
    cvx_end

    W_T(:,:,:,t+1)=W;
    Z_T(:,:,:,t+1)=Z;
    Alpha_T(t+1)=cvx_optval

    if Alpha_T(t+1)-Alpha_T(t)<=1e-3
%         abs(Alpha_T(t+1)-Alpha_T(t))<=1e-1
        break
    end
    t=t+1;
end

if Alpha_T(t+1)-Alpha_T(t)>=0
    W_out=W_T(:,:,:,t+1);
    Z_out=Z_T(:,:,:,t+1);
    else
    W_out=W_T(:,:,:,t);
    Z_out=Z_T(:,:,:,t);
end


% W_out=W_T(:,:,:,t+1);
% Z_out=Z_T(:,:,:,t+1);
end