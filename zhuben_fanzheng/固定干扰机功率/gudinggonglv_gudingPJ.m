T=70;
MMP=[];
fai=3;
y=1e+8;
VH_max=10;%���˻��ٶ�Լ��
% VV_max=5;
t=1;
k1=0.557;
PG_avg=2.7;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PG_max=2*PG_avg;
% PJ_avg=0.0007;
% PJ_max=4*PJ_avg;
HG=[0;0];
HJ_I=[-200;-100];%���˻�ˮƽ��㡢�յ�λ��
HJ_F=[200;-100];
% VJ_I=70;%���˻���ֱ��㡢�յ�λ��
% VJ_F=70;
HU=[-150,200;200,200];%�Ϸ��û�λ��
HE=[-100,50;-350,-380];%�����߹���λ��
HE_T=[-90,90;-360,-382];
QE=[10,40];%�������
H_NFZ=[-100,50;-200,-250];%���н�������λ��
Q_NFZ=[45,30];%���н����뾶
BCD=-20;
for i=1:T %��ʼ�����˻�ˮƽ�켣
    a=[-200+(400/T)*i;-100];
    HJ_0(:,i)=a;
end
VJ_0=70*ones(1,T);%��ʼ�����˻���ֱ�켣
PJ_0=0.0012*ones(1,T);%��ʼ�����˻�����
PG_0=0.6*ones(2,T);%��ʼ����վ����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

M_E0=160000*ones(2,T);

for BCD_N=1:100
   % �Ż����˻��켣
for m=1:2
    for n=1:T
         c_m(m,n)=exp(-k1)*PG_0(m,n)*y/((norm(HG-HU(:,m)))^(fai));
    end
end
for m=1:2
    for k=1:2
        for n=1:T
             e_m_k(m,k,n)=PG_0(m,n)*y/((norm(HG-HE(:,k))-QE(k))^(fai));
        end
    end
end
 
cvx1=-20;
for sub1=1:100
cvx_begin
variable eta
variable HJ(2,T)
% variable VJ(1,90)
variable D(2,T)
variable M_Um(2,T)
variable M_E(2,T)
variable lamda(2,T)
variable F_E(2,T)
expression a1(2,T)

for m=1:2
    for n=1:T
    a1(m,n)=log(1+c_m(m,n)-c_m(m,n)*PJ_0(n)*y*inv_pos(y*PJ_0(n)+M_Um(m,n)))/log(2)-D(m,n);
    end
end
maximize(eta)  
subject to
sum(a1(1,:))*(1/T)-eta>=0;
sum(a1(2,:))*(1/T)-eta>=0;
for m=1:2
    for n=1:T
        VJ_0(n)^(2)+(norm(HJ_0(:,n)-HU(:,m)))^(2)+2*(HJ(:,n)-HJ_0(:,n))'*(HJ_0(:,n)-HU(:,m))-M_Um(m,n)>=0; %VJ_0(n)^(2)+(VJ(n)-VJ_0(n))*2*VJ_0(n)+
    end
end
for m=1:2
    for k=1:2
        for n=1:T
            log(1+(e_m_k(m,k,n)/(1+PJ_0(n)*y/M_E0(k,n))))/log(2)+(M_E(k,n)-M_E0(k,n))*y*PJ_0(n)*e_m_k(m,k,n)/(log(2)*(PJ_0(n)*y+M_E0(k,n)+e_m_k(m,k,n)*M_E0(k,n))*(PJ_0(n)*y+M_E0(k,n)))-D(m,n)<=0;
        end
    end
end

for k=1:2
    for n=1:T
        [lamda(k,n)-1,0,HJ(1,n)-HE(1,k);0,lamda(k,n)-1,HJ(2,n)-HE(2,k);HJ(1,n)-HE(1,k),HJ(2,n)-HE(2,k),(-1)*lamda(k,n)*QE(k)^(2)+F_E(k,n)]==semidefinite(3);%%��ȷ��
    end
end

for k=1:2
    for n=1:T
        lamda(k,n)>=0;
    end
end

for k=1:2
    for n=1:T
        -(square_pos(norm(HJ(:,n)-HE(:,k)))+VJ_0(n)^(2)-M_E(k,n))-F_E(k,n)>=0;% -(square_pos(norm(HJ(:,n)-HE(:,k)))+VJ(n)^(2)-M_E(k,n))-F_E(k,n)>=0;
    end
end

for n=1:(T-1)
    norm(HJ(:,n+1)-HJ(:,n))-VH_max*t<=0;
end
% for n=1:(T-1)
%     norm(VJ(n+1)-VJ(n))-VV_max*t<=0;
% end
norm(HJ(:,1)-HJ_I)-VH_max*t<=0;
HJ(:,T)==HJ_F;
% norm(VJ(1)-VJ_I)<=VV_max*t;
% norm(VJ(T)-VJ_F)<=VV_max*t;
% 60*ones(1,T)<=VJ<=120*ones(1,T);%���и߶�Լ��
for d=1:2
    for n=1:T
        (norm(HJ_0(:,n)-H_NFZ(:,d)))^2+2*(HJ_0(:,n)-H_NFZ(:,d))'*(HJ(:,n)-HJ_0(:,n))>=Q_NFZ(d)^2;%���н���Լ��
    end
end
cvx_end
sub1_cishu=sub1
if (abs((cvx_optval-cvx1)/cvx_optval)<=0.001)
    break;
else
    cvx1=cvx_optval;
      HJ_0=HJ;
%       VJ_0=VJ;
      M_E0=M_E;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:2
y_m(m)=exp(-k1)*y/(norm(HG-HU(:,m)))^(fai);
end
for k=1:2
    y_e(k)=y/(norm(HG-HE(:,k))-QE(k))^(fai);
end
for m=1:2
    for n=1:T
        h_Jm(m,n)=y/(VJ_0(n)^(2)+(norm(HJ_0(:,n)-HU(:,m)))^(2));
    end
end
for k=1:2
    for n=1:T
        h_Je(k,n)=y/(VJ_0(n)^(2)+(norm(HJ_0(:,n)-HE(:,k))+QE(k))^(2));
    end
end
cvx2=-20;
for sub2=1:100
cvx_begin
variable eta
variable PG(2,T)
% variable PJ(1,90)
variable B(2,T)
expression a2(2,T)

for m=1:2
    for n=1:T
        a2(m,n)=log(1+PJ_0(n)*h_Jm(m,n)+y_m(m)*PG(m,n))/log(2)-log(1+PJ_0(n)*h_Jm(m,n))/log(2)-B(m,n);
    end
end
maximize(eta)  
subject to
sum(a2(1,:))*(1/T)-eta>=0;
sum(a2(2,:))*(1/T)-eta>=0;
for m=1:2
    for k=1:2
        for n=1:T
           log(1+PG_0(m,n)*y_e(k)/(1+PJ_0(n)*h_Je(k,n)))/log(2)+(PG(m,n)-PG_0(m,n))*y_e(k)/(1+PJ_0(n)*h_Je(k,n)+PG_0(m,n)*y_e(k))-B(m,n)<=0;    %log(1+PJ_0(n)* h_Je(k,n)+PG_0(m,n)*y_e(k))/log(2)+y_e(k)*(PG(m,n)-PG_0(m,n))/(log(2)*(1+PJ_0(n)* h_Je(k,n)+PG_0(m,n)*y_e(k)))+h_Je(k,n)*(PJ(n)-PJ_0(n))/(log(2)*(1+PJ_0(n)* h_Je(k,n)+PG_0(m,n)*y_e(k)))-log(1+PJ(n)*h_Je(k,n))/log(2)-B(m,n)<=0;
        end
    end
end

PG<=PG_max*ones(2,T);
PG>=zeros(2,T);
sum(PG(:))*(1/T)-PG_avg<=0;
% PJ<=PJ_max*ones(1,T);
% PJ>=zeros(1,T);
% sum(PJ)*(1/T)-PJ_avg<=0;
cvx_end
sub2_cishu=sub2
if (abs((cvx_optval-cvx2)/cvx_optval)<=0.01)
    break;
else
    cvx2=cvx_optval;
%     PJ_0=PJ;
    PG_0=PG;
end
end

BCD_cishu=BCD_N %(显示总问题迭代次�?
if(abs((cvx2-BCD)/cvx2)<=0.01)
 break;
else
    BCD=cvx2;
    MMP(BCD_cishu)=BCD;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
R_e_m=[];
  R_e=[];
   R_m=[];
    R_sec=[];
for m=1:2
    for k=1:2
        for n=1:T
           R_e_m(m,k,n)=log2(1+PG(m,n)*y*(norm(HG-HE_T(:,k)))^(-fai)/(1+PJ_0(n)*y/(VJ_0(n)^(2)+(norm(HJ(:,n)-HE_T(:,k)))^(2)))); 
        end
    end
end

for m=1:2
        for n=1:T
            if R_e_m(m,1,n)>=R_e_m(m,2,n)
                R_e(m,n)=R_e_m(m,1,n);
            else
                 R_e(m,n)=R_e_m(m,2,n);
            end
        end
end
for m=1:2
    for n=1:T
        R_m(m,n)=log2(1+exp(-k1)*PG(m,n)*y*(norm(HG-HU(:,m)))^(-fai)/(1+PJ_0(n)*y/(VJ_0(n)^(2)+(norm(HJ(:,n)-HU(:,m)))^(2)))); 
    end
end
for m=1:2
    for n=1:T
        R_sec(m,n)=R_m(m,n)-R_e(m,n);
    end
end
R_avg=sum( R_sec(:))*(1/T)*0.5
% xlswrite('gudinPJ_R_avg_90.xlsx',R_avg)
% xlswrite('gudinPJ_T90.xlsx',HJ)
xlswrite('gudinPJ_R_avg_PG_2.7_0.0012.xlsx',R_avg)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% xlswrite('gudinPJ_R_avg_PJ_0.0027.xlsx',R_avg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
