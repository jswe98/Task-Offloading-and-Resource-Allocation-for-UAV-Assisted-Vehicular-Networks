%确定路损为fai，噪声为sigma，单位距离信道增益为lou，lou/sigma为y，无人机水平飞行�?��速度为VH_max，无人机垂直飞行�?��速度为VV_max.时隙为d，欧拉常数为k，无人机初始功率为PJ_0,基站初始功率为PD_0,无人机的初始水平轨迹为HJ_0,
%初始垂直轨迹为VJ_0,基站坐标为HG，用�?的位置为HU1，用�?的位置为HU2，窃听�?1的估计位置为HE1，窃听�?2的估计位置为HE2，窃听�?1的误差半径为Q1,
%窃听�?的误差半径为Q2，禁飞区1的半径为Q_NFZ1,其水平投影圆心坐标为H_NFZ1,禁飞�?的半径为Q_NFZ2,其水平投影圆心坐标为H_NFZ2,
%基站平均功率为PG_avg，最大功率为PG_max;无人机平均与�?��为PJ_avg,PJ_max
lou=1e-7;
sigma=1e-14;
fai=3;
y=1e+7;
VH_max=7;
VV_max=5;
t=1;
k=0.557;
PG_avg=4;
PG_max=4*PG_avg;
PJ_avg=0.01;
PJ_max=4*PJ_avg;
HG=[0;0];
HJ_I=[-210;-100];
HJ_F=[210;-100];
VJ_I=70;
VJ_F=70;
HU1=[-100;100];
HU2=[100;100];
HE1=[-100;-200];
HE2=[100;-200];
Q1=20;
Q2=10;
H_NFZ1=[-150;-150];
Q_NFZ1=10;
H_NFZ2=[150;-150];
Q_NFZ2=10;
cvx1=0;
cvx2=0;
cvx3=0;
BCD=0;
%初始化无人机轨迹
for i=1:120
    a1=[-210+3.5*i;-100];
    HJ_0(:,i)=a1;
end
VJ_0=80*ones(1,120);
%初始化无人机发射功率
PJ_0=0.008*ones(1,120);
%初始化基站的发射功率
PD_0=0.1*ones(4,120);

%%%%%%M2_E_0(2,120)的初始�?%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:120
    M2_E_0(1,n)=VJ_0(n)^(2)+(norm(HJ_0(:,n)-HE1)+Q1)^(2);%(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
end

for n=1:120
    M2_E_0(2,n)=VJ_0(n)^(2)+(norm(HJ_0(:,n)-HE2)+Q2)^(2);
end
 %%%%%%%%%%%%%%求解第一个子问题(固定PD,PJ,HJ,VJ)%%%%%%%%%%%%%%%%
for BCD_N=1:10000
 for sub1=1:1000
 d_GE1xin=norm(HG-HE1)-Q1;
 d_GE2xin=norm(HG-HE2)-Q2;
for n=1:120
    h_JE1xin(n)=lou/(VJ_0(n)^(2)+(norm(HJ_0(:,n)-HE1)+Q1)^(2));
end
for n=1:120
    h_JE2xin(n)=lou/(VJ_0(n)^(2)+(norm(HJ_0(:,n)-HE2)+Q2)^(2));
end

d_GU1=norm(HG-HU1);
d_GU2=norm(HG-HU2);
for n=1:120
h_JU1(n)=lou/(VJ_0(n)^(2)+(norm(HJ_0(:,n)-HU1))^(2));
end
for n=1:120
h_JU2(n)=lou/(VJ_0(n)^(2)+(norm(HJ_0(:,n)-HU2))^(2));
end
%在轨迹和功率确定的清况下寻找窃听速率�?��的窃听�?
%窃听�?在每个时隙下的每个子载波的最大窃听�?率（不�?虑A�?
for n=1:120
    for d=1:4
        R_E1_D(d,n)=log2(1+(PD_0(d,n)*lou*d_GE1xin^(-fai))/(sigma+PJ_0(n)*h_JE1xin(n)));
    end
end
%窃听�?在每个时隙下的每个子载波的最大窃听�?率（不�?虑A�?
for n=1:120
    for d=1:4
        R_E2_D(d,n)=log2(1+(PD_0(d,n)*lou*d_GE2xin^(-fai))/(sigma+PJ_0(n)*h_JE2xin(n)));
    end
end
%找出每个时隙下，每个子载波下的窃听�?率最大的窃听�?
for n=1:120
    for d=1:4
        if (R_E1_D(d,n)>=R_E2_D(d,n))
            R_K_MAX(d,n)=R_E1_D(d,n);
        else
            R_K_MAX(d,n)=R_E2_D(d,n);
        end
    end
end

%第一个合法接收�?在每个时隙下的，每个子载波下的传输�?�?
for n=1:120
    for d=1:4
        R_U1_D(d,n)=log2(1+(exp(-k)*PD_0(d,n)*lou*d_GU1^(-fai))/(sigma+PJ_0(n)*h_JU1(n)));
    end
end
%第二个合法接收�?在每个时隙下的，每个子载波下的传输�?�?
for n=1:120
    for d=1:4
        R_U2_D(d,n)=log2(1+(exp(-k)*PD_0(d,n)*lou*d_GU2^(-fai))/(sigma+PJ_0(n)*h_JU2(n)));
    end
end

cvx_begin      
variable eta
variable A(2,4,120)
expression a1(4,120)
expression a2(4,120)
expression a3(2,4,120)
expression a4(4,120)
    for n=1:120
        for d=1:4
            a1(d,n)=A(1,d,n)*(R_U1_D(d,n)-R_K_MAX(d,n));
        end
    end
    
    for n=1:120
        for d=1:4
            a2(d,n)=A(2,d,n)*(R_U2_D(d,n)-R_K_MAX(d,n));
        end
    end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for n=1:120
          for m=1:2
              for d=1:4
                  a3(m,d,n)=(1/120)*A(m,d,n)*PD_0(d,n);
              end
          end
      end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      for n=1:120
          for d=1:4
              a4(d,n)=sum(A(:,d,n));%%%(三维矩阵求任意维度的�?
          end
      end
      
      maximize(eta) 
      subject to 
      sum(a1(:))-eta>=0;
      sum(a2(:))-eta>=0;
      sum(squeeze(sum(sum(a3))))-PG_avg<=0;
      a4>=zeros(4,120);
      a4<=ones(4,120);
      A>=zeros(2,4,120);
      A<=ones(2,4,120);
      cvx_end
      A_0=A;
      sub1_cishu=sub1 %(显示子问�?迭代次数)
      if ((cvx_optval-cvx1)<=0.001)
    break;
else
    cvx1=cvx_optval;
      end
end
          
      
     
 %%%%%%%%%%%%%%%%%%求解第二个子问题(固定A,HJ,VJ)%%%%%%%%%%%%%%
 for sub2=1:1000
 y_U1=exp(-k)*lou/(d_GU1^(fai));
 y_U2=exp(-k)*lou/(d_GU2^(fai));
 y_E1=lou/(d_GE1xin^(fai));
 y_E2=lou/(d_GE2xin^(fai));
 cvx_begin
 variable eta
 variable PJ(1,120)
 variable PD(4,120)
 variable B(4,120)
 expression a5(4,120)
 expression a6(4,120)
 expression a7(2,4,120)
 for n=1:120
     for d=1:4
         a5(d,n)=A_0(1,d,n)*(log(y_U1*PD(d,n)+sigma+PJ(n)*h_JU1(n))/log(2)-log(sigma+PJ_0(n)*h_JU1(n))/log(2)-h_JU1(n)*(PJ(n)-PJ_0(n))/(log(2)*(sigma+PJ_0(n)*h_JU1(n)))-B(d,n));
     end
 end

for n=1:120
     for d=1:4
         a6(d,n)=A_0(2,d,n)*(log(y_U2*PD(d,n)+sigma+PJ(n)*h_JU2(n))/log(2)-log(sigma+PJ_0(n)*h_JU2(n))/log(2)-h_JU2(n)*(PJ(n)-PJ_0(n))/(log(2)*(sigma+PJ_0(n)*h_JU2(n)))-B(d,n));
     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for n=1:120
          for m=1:2
              for d=1:4
                  a7(m,d,n)=(1/120)*A_0(m,d,n)*PD(d,n);
              end
          end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 maximize(eta) 
 subject to 
 sum(a5(:))-eta>=0;
 sum(a6(:))-eta>=0;
 for n=1:120
    for d=1:4
     log(y_E1*PD_0(d,n)+sigma+PJ_0(n)*h_JE1xin(n))/log(2)+y_E1*(PD(d,n)-PD_0(d,n))/(log(2)*(y_E1*PD_0(d,n)+sigma+PJ_0(n)*h_JE1xin(n)))+h_JE1xin(n)*(PJ(n)-PJ_0(n))/(log(2)*(y_E1*PD_0(d,n)+sigma+PJ_0(n)*h_JE1xin(n)))-log(sigma+PJ(n)*h_JE1xin(n))/log(2)-B(d,n)<=0;
    end
 end                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

for n=1:120
    for d=1:4
     log(y_E2*PD_0(d,n)+sigma+PJ_0(n)*h_JE2xin(n))/log(2)+y_E2*(PD(d,n)-PD_0(d,n))/(log(2)*(y_E2*PD_0(d,n)+sigma+PJ_0(n)*h_JE2xin(n)))+h_JE2xin(n)*(PJ(n)-PJ_0(n))/(log(2)*(y_E2*PD_0(d,n)+sigma+PJ_0(n)*h_JE2xin(n)))-log(sigma+PJ(n)*h_JE2xin(n))/log(2)-B(d,n)<=0;
    end
end
 sum(squeeze(sum(sum(a7))))-PG_avg<=0;
 zeros(4,120)<=PD;
 PD<=PG_max*ones(4,120);
 zeros(1,120)<=PJ;
 PJ<=PJ_max*ones(1,120);
 1/120*sum(PJ)-PJ_avg<=0;
  cvx_end
  PJ_0=PJ;
  PD_0=PD;
  sub2_cishu=sub2 %(显示子问�?迭代次数)
  if ((cvx_optval-cvx2)<=0.001)
    break;
else
    cvx2=cvx_optval;
      end
 end
 
%%%%%%%%%%%%%%%%%%%%%求解第三个子问题（固定A,PJ,PD)%%%%%%%%%%%%%%%%%%
for sub3=1:1000
cvx_begin
variable eta
variable M1_U(2,120)
variable M2_E(2,120)
variable Dd(4,120)
variable HJ(2,120)
variable VJ(1,120)
variable lamda_E(2,120)
variable F_E(2,120)
expression st8(1,2)
expression st9(2,120)
expression st10(2,4,120)
expression st11(2,120)
expression st12(1,119)
expression st13(1,119)
expression st14(2,120)
for n=1:120
    for d=1:4
        a(d,n)=A_0(1,d,n)*(log(1+(exp(-k)*PD_0(d,n)*lou/(norm(HG-HU1))^(fai))/(sigma+(lou*PJ_0(n))/(M1_U(1,n))))/log(2)-Dd(d,n));
    end
end
st8(1,1)=sum(a(:))-eta;

for n=1:120
    for d=1:4
        a(d,n)=A_0(2,d,n)*(log(1+(exp(-k)*PD_0(d,n)*lou/(norm(HG-HU2))^(fai))/(sigma+(lou*PJ_0(n))/(M1_U(2,n))))/log(2)-Dd(d,n));
    end
end
st8(1,2)=sum(a(:))-eta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:120
    st9(1,n)=VJ_0(n)^(2)+2*VJ_0(n)*(VJ(n)-VJ_0(n))+(norm(HJ_0(:,n)-HU1))^(2)+2*(HJ(:,n)-HJ_0(:,n))'*(HJ_0(:,n)-HU1)-M1_U(1,n);
end

for n=1:120
    st9(2,n)=VJ_0(n)^(2)+2*VJ_0(n)*(VJ(n)-VJ_0(n))+(norm(HJ_0(:,n)-HU2))^(2)+2*(HJ(:,n)-HJ_0(:,n))'*(HJ_0(:,n)-HU2)-M1_U(2,n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:120
    for d=1:4
        st10(1,d,n)=log(1+(PD_0(d,n)*lou/(norm(HG-HE1)-Q1)^(fai))/(sigma+PJ_0(n)*lou/M2_E_0(1,n)))/log(2)+(M2_E(1,n)-M2_E_0(1,n))*(PD_0(d,n)*lou/(norm(HG-HE1)-Q1)^(fai)*PJ_0(n)*lou)/(log(2)*(PD_0(d,n)*lou/(norm(HG-HE1)-Q1)^(fai)*M2_E_0(1,n)+sigma*M2_E_0(1,n)+PJ_0(n)*lou)*(sigma*M2_E_0(1,n)+PJ_0(n)*lou))-Dd(d,n);
    end
end

for n=1:120
    for d=1:4
        st10(2,d,n)=log(1+(PD_0(d,n)*lou/(norm(HG-HE2)-Q2)^(fai))/(sigma+PJ_0(n)*lou/M2_E_0(2,n)))/log(2)+(M2_E(2,n)-M2_E_0(2,n))*(PD_0(d,n)*lou/(norm(HG-HE2)-Q2)^(fai)*PJ_0(n)*lou)/(log(2)*(PD_0(d,n)*lou/(norm(HG-HE2)-Q2)^(fai)*M2_E_0(2,n)+sigma*M2_E_0(2,n)+PJ_0(n)*lou)*(sigma*M2_E_0(2,n)+PJ_0(n)*lou))-Dd(d,n);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:120
    st11(1,n)=-((norm(HJ(:,n)-HE1))^(2)+VJ(n)^(2)-M2_E(1,n))-F_E(1,n);
end

for n=1:120
    st11(2,n)=-((norm(HJ(:,n)-HE2))^(2)+VJ(n)^(2)-M2_E(2,n))-F_E(2,n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:119
    st12(i)=norm(HJ(:,n+1)-HJ(:,n))-VH_max*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:119
    st13(i)=(VJ(n+1)-VJ(n))-VV_max*t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:120
    st14(1,n)=(norm(HJ_0(:,n)-H_NFZ1))^(2)+2*(HJ(:,n)-HJ_0(:,n))'*(HJ_(:,n)-H_NFZ1)-(Q_NFZ1)^(2);
end

for n=1:120
    st14(2,n)=(norm(HJ_0(:,n)-H_NFZ2))^(2)+2*(HJ(:,n)-HJ_0(:,n))'*(HJ_(:,n)-H_NFZ2)-(Q_NFZ2)^(2);
end

maximize(eta)  
 subject to 
 st8>=zeros(1,2);
 st9>=zeros(2,120);
 st10<=zeros(2,4,120);
 st11>=zeros(2,120);
 st12<=zeros(1,119);
 st13<=zeros(1,119);
 st14>=zeros(2,120);
 
 for n=1:120
     [lamda_E(1,n)-1,0,HJ(1,n)-HE1(1,1);0,lamda_E(1,n)-1,HJ(2,n)-HE1(2,1);HJ(1,n)-HE1(1,1),HJ(2,n)-HE1(2,1),(-1)*lamda_E(1,n)*Q1^(2)+F_E(1,n)]==semidefinite(3);
 end
 
  for n=1:120
     [lamda_E(2,n)-1,0,HJ(1,n)-HE2(1,1);0,lamda_E(2,n)-1,HJ(2,n)-HE2(2,1);HJ(1,n)-HE2(1,1),HJ(2,n)-HE2(2,1),(-1)*lamda_E(2,n)*Q2^(2)+F_E(2,n)]==semidefinite(3);
  end
  
lamda_E>=zeros(2,120);
norm(HJ(:,1)-HJ_I)-VH_max*t<=0;
norm(HJ(:,120)-HJ_F)-VH_max*t<=0;
VJ(1)-VJ_I<=VV_max*t;
VJ(120)-VJ_F<=VV_max*t;
20*ones(1,120)<=VJ<=120*ones(1,120);
cvx_end

HJ_0=HJ;
VJ_0=VJ;
M2_E_0=M2_E;
sub3_cishu=sub3 %(显示子问�?迭代次数)
 if ((cvx_optval-cvx3)<=0.001)
    break;
else
    cvx3=cvx_optval;
      end
end
 BCD_cishu=BCD_N %(显示总问题迭代次�?
if(eta-BCD<=0.001)
 break;
else
    BCD=eta;
end
end


