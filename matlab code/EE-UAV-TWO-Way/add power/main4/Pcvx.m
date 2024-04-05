function [P_m]  = Pcvx(W,X,Y,Delta,GVR,GVU,theta,T,M,chi2)  
 GVRR=GVR';
GVUU=GVU';
c1=zeros(T,M);
c2=zeros(T,M);
 SNRth=1;
 e1=0.1;%第1个阈值
 Pmax_m=0.5*ones(T,M);%车车的最大功率(W)            
% chi1=0;%这个以后要优化的丁克尔巴赫的系数
% cvx_clear                      
cvx_solver SeDuMi
cvx_precision best
cvx_begin
           variable P_m(T,M)
           

    for t = 1:T %给出信道增益h(t,m)和c2(t,m)
        for m = 1:M
            c1(t,m) = GVRR(t,m)/Delta;
            c2(t,m) =  GVUU(t,m)/Delta;%原
        end
    end  %这里是为了求解轨迹规划后新的信道状态信息拿来优化Pm   1,1+P_b.*c1
           
           L_VRR=sum(sum(W*X'.*rel_entr(ones(T,M),1+P_m.*c1)));
           L_AA=sum(sum(W*Y'.*rel_entr(ones(T,M),1+P_m.*c2)));

%            L_A=sum(sum(W*Y'.*(rho_l+omega_l.*(distance2-distance_l_SCA))));%凸+凹？
           E_VRR=sum(sum(P_m.*X'));%VRde能量消耗
           E_VUU=sum(sum(P_m.*Y'));%VUde能量消耗  -sum(sum(W*X'.*rel_entr(ones(T,M),1+P_m.*c1)))-sum(sum(W*Y'.*rel_entr(ones(T,M),1+P_m.*c2)))
            minimize(0+chi2*1*(theta* E_VRR+(1-theta)*E_VUU));  

           subject to
          GVRR.* P_m.*X' >= SNRth*Delta+ GVUU.*P_m.*Y'+log(1-e1)*P_m.*X';
             P_m >= 0;%<= Pmax_b;
             P_m <= Pmax_m;      
cvx_end   