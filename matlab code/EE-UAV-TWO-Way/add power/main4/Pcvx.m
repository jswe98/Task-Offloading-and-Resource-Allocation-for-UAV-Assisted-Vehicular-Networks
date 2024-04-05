function [P_m]  = Pcvx(W,X,Y,Delta,GVR,GVU,theta,T,M,chi2)  
 GVRR=GVR';
GVUU=GVU';
c1=zeros(T,M);
c2=zeros(T,M);
 SNRth=1;
 e1=0.1;%��1����ֵ
 Pmax_m=0.5*ones(T,M);%�����������(W)            
% chi1=0;%����Ժ�Ҫ�Ż��Ķ��˶��ͺյ�ϵ��
% cvx_clear                      
cvx_solver SeDuMi
cvx_precision best
cvx_begin
           variable P_m(T,M)
           

    for t = 1:T %�����ŵ�����h(t,m)��c2(t,m)
        for m = 1:M
            c1(t,m) = GVRR(t,m)/Delta;
            c2(t,m) =  GVUU(t,m)/Delta;%ԭ
        end
    end  %������Ϊ�����켣�滮���µ��ŵ�״̬��Ϣ�����Ż�Pm   1,1+P_b.*c1
           
           L_VRR=sum(sum(W*X'.*rel_entr(ones(T,M),1+P_m.*c1)));
           L_AA=sum(sum(W*Y'.*rel_entr(ones(T,M),1+P_m.*c2)));

%            L_A=sum(sum(W*Y'.*(rho_l+omega_l.*(distance2-distance_l_SCA))));%͹+����
           E_VRR=sum(sum(P_m.*X'));%VRde��������
           E_VUU=sum(sum(P_m.*Y'));%VUde��������  -sum(sum(W*X'.*rel_entr(ones(T,M),1+P_m.*c1)))-sum(sum(W*Y'.*rel_entr(ones(T,M),1+P_m.*c2)))
            minimize(0+chi2*1*(theta* E_VRR+(1-theta)*E_VUU));  

           subject to
          GVRR.* P_m.*X' >= SNRth*Delta+ GVUU.*P_m.*Y'+log(1-e1)*P_m.*X';
             P_m >= 0;%<= Pmax_b;
             P_m <= Pmax_m;      
cvx_end   