 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5.21,9:57%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ȫ�߹�ͨ,%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%δ�ܴ�ѭ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% E2RA
%% April 28 beginning
%
%% ��������
%{
l(M,1)%�������븡��ľ���
Wm %sink���� 
cmar %���Ӿ���connective martrix
cmi,cmj:���Ӿ���ļ�������
h(t,m) %���ߵ��ŵ�����
c2(t,m) &ˮ����·�ŵ�����
chi %���˶��ͺ�����ֵ
chi_rec %chi��record����
%}
clear
%% �ܳ�ʼ��
M=5;%sink����
T=90;%ʱ϶����
Ls=3000;%�����ģscale(m)
Vmax=80;
beta0=10^(-6);%��λ�ŵ�����
sigma2=10^(-14);%��������
deltaf=1000;%ˮ������
%af=2.59;%f=20kHz
af=4.1;%f=25kHz
kappa=1.5;%����
%Nf=10^(-14.2);%��������(7��2ͼ�ϵ�30����)
Nf=10^(-12.5);%��������(7��2ͼ�ϵ�50����)���ܺ�һЩ
W=10^(6);%���ߵ����1MHz
%Tt=100;%�������ڣ�s��
tau=1;%ʱ϶���ȣ�s��
H=100;%UAV���и߶ȣ�m��
Pr=8*10^(-4)*ones(T,M);%����������ݹ��ʣ�W��
%% ϵͳ�ɱ��ʼ��
theta=0.4;%������ֵȨ��
% L_buoy=300*1000*ones(1,M);%��ʱ�����и�������������&����ps����ɹ��Ĳ���
% L_sen=100*1000*ones(1,M);%��ʱ�����и�������������
L_buoy=300*1000*ones(1,M);%��ʱ�����и�������������
L_sen=100*1000*ones(1,M);%��ʱ�����и�������������
E_bub=0.02*ones(1,M);%buoy power budget����Ԥ��
E_bus=0.02*ones(1,M);%sensor power budget����Ԥ��
Pmax_b=0.1*ones(T,M);%buoy�����(W)
Pmax_s=0.1*ones(T,M);%sensor�����
L=1;%BCD�ܵ�������
epsilon2=0.1;%�����ж���ֵ
epsilon1=0.1;
%������Ϲ������
%{
theta=0.4;%������ֵȨ��
L_buoy=2000*1000*ones(1,M);%��ʱ�����и�������������
L_sen=100*1000*ones(1,M);%��ʱ�����и�������������
E_bub=0.01*ones(1,M);%buoy power budget����Ԥ��
E_bus=0.005*ones(1,M);%sensor power budget����Ԥ��
Pmax_b=0.1*ones(T,M);%buoy�����(W)
Pmax_s=0.1*ones(T,M);%sensor�����
L=1;%BCD�ܵ�������
epsilon=0.1;%�����ж���ֵ
%}
%ԭ����
%% �ڵ�ֲ�
W_m = [Ls/22.2,   Ls/3.98,  Ls/1.95,  Ls/1.41,  Ls/2.3;
       Ls/2.2,    Ls/1.28, Ls/1.23,  Ls/2.9,   Ls/6.2];
   geocent(1,1) = sum(W_m(1,:))/M;%��������geometry centre
   geocent(1,2) = sum(W_m(2,:))/M;
  
for m=1:M % ����ڵ�
    %{
%     W_m(1,m) = 800;
%     W_m(2,m) = 500;
    W_m(1,m) = Ls*rand;
    W_m(2,m) = Ls*rand;%sink�����������
    %}
    g2b(1,m)=sqrt((geocent(1,1)-W_m(1,m))^(2)+(geocent(1,2)-W_m(2,m))^(2));%geometry centre to buoys
    %����sinkλ��
    %{
    scatter (W_m(1,m),W_m(2,m),'mo');
    hold on
    %}
    %{
    gplot (cmar,q,'g<-');%�Ż��켣
    gplot (cmar,q_init,'g<-');
    %}
end
r_init = 0.8*max(g2b);
%% s2b ����
%l_s2b = 0.1;
%l_s2b = 1+rand(1,M);%�������ڵ㵽��Ӧ����ľ���(ˮ����km)
%l_s2b = [1.3816, 1.7655, 1.7952, 1.1869, 1.4898];
l_s2b = [2.3816, 1.9655, 2.2952, 2.1869, 2.4898];
%l_s2b = [2.5816, 2.9655, 2.3952, 2.7869, 2.4898];
%{
for m=1:M %�ṹ��
    sink(m).Wm_x = W_m(1,m);
    sink(m).Wm_y = W_m(2,m);%sinkλ�ô���ṹ��
    sink(m).se_bu = l_s2b(1,m);%sink--sensor�������ṹ��
end
%}
%% BCD��ʼ��
for m=1:M %��ʼ��alpha&beta������
    for t=1:T
        alpha_init(t,m)=1/M;%����ITʱ��
        beta_init(t,m)=0.09;%������ITʱ��
    end
end
for t=1:T %��ʼ���켣Q:q_init
    x_init = geocent(1,1)+r_init*cos(t*2*pi/T);
    y_init = geocent(1,2)+r_init*sin(t*2*pi/T);
    q_init(t,1)=x_init;%q�����������T*2��ʽ
    q_init(t,2)=y_init;
end
cmar=zeros(T,T);%�������Ӿ���
for cmi=1:T
    for cmj=1:T
        if (abs(cmi-cmj)==1||abs(cmi-cmj)==T-1)
            cmar(cmi,cmj)=1;
        end
    end
end
%���Ƴ�ʼ�켣
%{
gplot (cmar,q_init,'g<-');%��ʼ�켣��cmar�����Ӿ���connective martrix
gplot (cmar,q,'g<-');%�Ż��켣
%}
alpha1(:,:) = alpha_init;
beta1(:,:) = beta_init;
P_b = 1*10^(-4)*ones(T,M);
P_s = 1*10^(-3)*ones(T,M);

% P_b = 10^(-4)*ones(T,M);%5.18,16:51ԭ����,�������ϸĻأ���ʱQ����
% P_s = 10^(-3)*ones(T,M);
q_l=q_init;

%% BCD��ʼ
 
flagBCD=0;
l=1;%BCD��������
EE_total_recod(1)=0;
 while (flagBCD==0)
    flagBCD = 1;
    if (l>1)
        q_l=storage_BCD(l-1).Q;
        alpha1=storage_BCD(l-1).A;
        beta1=storage_BCD(l-1).B;
       % P_s=storage_BCD(l-1).Ps;
        P_b=storage_BCD(l-1).Pb;
    end
    %%%%%%%n_3%%%%%%%%%%%%%%%%%%%%%����P^s,P^b,A,B��������Q  
    %
    n_3=1;
    flag=0;
    q_l_SCA=q_l;
    while (flag==0)
        flag=1;
        for t = 1:T %Qcvx��������Ķ���ֵ
            for m = 1:M
                c2(t,m) = 1/(l_s2b(1,m)^(kappa)*af^(l_s2b(1,m))*Nf*deltaf);%ԭ
                c3_l(t,m)=P_b(t,m)/sigma2;
                omega_l(t,m)=(-1*c3_l(t,m)*beta0)/(((q_l_SCA(t,1)-W_m(1,m))^(2)+(q_l_SCA(t,2)-W_m(2,m))^(2)+H^(2))*log(2)*...
                    ((q_l_SCA(t,1)-W_m(1,m))^(2)+(q_l_SCA(t,2)-W_m(2,m))^(2)+H^(2)+c3_l(t,m)*beta0));
                rho_l(t,m)=log2(1+c3_l(t,m)*beta0/((q_l_SCA(t,1)-W_m(1,m))^(2)+(q_l_SCA(t,2)-W_m(2,m))^(2)+H^(2)));
                distance_l(t,m)=(q_l_SCA(t,1)-W_m(1,m))^(2)+(q_l_SCA(t,2)-W_m(2,m))^(2);
                %��l�ε�������t��ʱ϶UAV���m������ˮƽ����ƽ��
            end
        end
        distance_l_SCA=distance_l;
        L_Ul=sum(sum(deltaf*beta1.*log2(1+P_s.*c2)));
        E_tbl=sum(sum(P_b.*alpha1));
        E_rbl=sum(sum(Pr.*beta1));
        E_tsl=sum(sum(P_s.*beta1));
        %
        [q,distance2] = Qcvx(W,alpha1,rho_l,omega_l,distance_l_SCA,W_m,L_Ul,E_tbl,theta,...
            E_rbl,E_tsl,L_buoy,Vmax,tau,T,M) ;
        q_l_SCA=q;%��¼�����Ż������������ʹ��.Ǳ����в��ʵ�ڲ��в����distance2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%Q����ʷ���
        for t = 1:T %Qcvx��������Ķ���ֵ
            for m = 1:M
                h(t,m) = beta0/(H^(2)+(q_l_SCA(t,1)-W_m(1,m))^(2)+(q_l_SCA(t,2)-W_m(2,m))^(2));
                c1(t,m) = h(t,m)/sigma2;
            end
        end
        L_buoy_real = sum(W*alpha1.*log2(1+P_b.*c1));
        L_A_lb = sum(W*alpha1.*(rho_l+omega_l.*(distance2-distance_l_SCA)));%lower bound
        L_sen_real = sum(deltaf*beta1.*log2(1+P_s.*c2));
        E_buoy_real = sum(P_b.*alpha1+Pr.*beta1);
        E_sen_real = sum(P_s.*beta1);
        ob(n_3) = (sum(L_buoy_real)+L_Ul)/(theta*(E_tbl+E_rbl)+(1-theta)*E_tsl);%43��Fig.2��eta
        ob_lb(n_3) = (sum(L_A_lb)+L_Ul)/(theta*(E_tbl+E_rbl)+(1-theta)*E_tsl);
        if (abs(sum(ob(n_3)-ob_lb(n_3)))/sum(ob_lb(n_3))>0.01)
            flag=0;
        end
        n_3=n_3+1;
    end %��Ӧwhile
    storage_BCD(l).Q = q_l_SCA;
    q_l = storage_BCD(l).Q;
    %}
        
        
        
       
    
            %%%%%%%%%%%%%n_2%%%%%%chi2%%%%%%%%%%%%���P^s,P^b,Q��������A,B(linear program)
    %
    
    for t = 1:T %�����ŵ�����h(t,m)��c2(t,m)
        for m = 1:M
            h(t,m) = beta0/(H^(2)+(q_l(t,1)-W_m(1,m))^(2)+(q_l(t,2)-W_m(2,m))^(2));
            c1(t,m) = h(t,m)/sigma2;
            c2(t,m) = 1/(l_s2b(1,m)^(kappa)*af^(l_s2b(1,m))*Nf*deltaf);%ԭ
        end
    end
    %
    n_2=1; %���A,B��dinkelbach�ִ�
    flag=0;
    chi2=0;    
    while (flag==0)%Dinkelbach
        flag=1;
        chi2_rec(n_2)=chi2;
        [alpha1,beta1]=ABcvx(W,deltaf,c1,c2,chi2,theta,...
                                Pr,L_buoy,L_sen,E_bub,E_bus,P_b,P_s,T,M);%����CVX
        F_chi2=sum(sum(-W*alpha1.*rel_entr(1,1+P_b.*c1)))/0.6931...
                     +sum(sum(-deltaf*beta1.*rel_entr(1,1+P_s.*c2)))/0.6931...
                     -chi2*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1));%����F_chi
        chi2=(sum(sum(-1*W*alpha1.*rel_entr(1,1+P_b.*c1)))/0.6931...
                     +sum(sum(-1*deltaf*beta1.*rel_entr(1,1+P_s.*c2)))/0.6931)...
                     /sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1));
         n_2=n_2+1;
         if (F_chi2>epsilon2)
             flag=0;
         end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%A&B����ʷ���
            L_buoy_real = sum(W*alpha1.*log2(1+P_b.*c1));
            L_sen_real = sum(deltaf*beta1.*log2(1+P_s.*c2));
            E_buoy_real = sum(P_b.*alpha1+Pr.*beta1);
            E_sen_real = sum(P_s.*beta1);
            storage_BCD(l).A = alpha1;
            storage_BCD(l).B = beta1;
    %}
    

    
    %%%%%%%%%%chi1%%%%%F_chi1%%%%%%%%%%%%%%%%%%%%%����A,B,Q��������P^s,P^b
    %%%%%%%%%%%chi1%%%%%%%%%%%%%%%%%%%%%%%%%%����A,B,Q��������P^b
    %
    for t = 1:T %�����ŵ�����h(t,m)��c2(t,m)
        for m = 1:M
            h(t,m) = beta0/(H^(2)+(q_l(t,1)-W_m(1,m))^(2)+(q_l(t,2)-W_m(2,m))^(2));
            c1(t,m) = h(t,m)/sigma2;
            c2(t,m) = 1/(l_s2b(1,m)^(kappa)*af^(l_s2b(1,m))*Nf*deltaf);%ԭ
        end
    end
    %
    n_1=1;%��⹦�ʵ�dinkelbach�ִ�
    flag=0;
    chi1=0;
    while (flag==0)%Dinkelbach
        flag=1;
        chi1_rec(n_1)=chi1;%chi1��record����

        
        [P_b] = Pcvx_b(W,deltaf,alpha1,beta1,c1,c2,chi1,theta,...
                          Pr,L_buoy,L_sen,E_bub,E_bus,Pmax_s,Pmax_b,T,M,P_s);
%         [P_b,P_s] = Pcvx(W,deltaf,alpha1,beta1,c1,c2,chi1,theta,...
%                          Pr,L_buoy,L_sen,E_bub,E_bus,Pmax_s,Pmax_b,T,M);%����CVX��ԭ������
                      
         F_chi1=sum(sum(-W*alpha1.*rel_entr(1,1+P_b.*c1)))/0.6931...
                     +sum(sum(-deltaf*beta1.*rel_entr(1,1+P_s.*c2)))/0.6931...
                     -chi1*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1));%����F_chi
         chi1=(sum(sum(-1*W*alpha1.*rel_entr(1,1+P_b.*c1)))/0.6931...
                     +sum(sum(-1*deltaf*beta1.*rel_entr(1,1+P_s.*c2)))/0.6931)...
                     /sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1));
        n_1=n_1+1;
         if (F_chi1>epsilon1)
             flag=0;
         end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%P����ʷ���
            for t = 1:T %Qcvx��������Ķ���ֵ
                for m = 1:M
                    h(t,m) = beta0/(H^(2)+(q_l(t,1)-W_m(1,m))^(2)+(q_l(t,2)-W_m(2,m))^(2));
                    c1(t,m) = h(t,m)/sigma2;
                end
            end            
            L_buoy_real = sum(W*alpha1.*log2(1+P_b.*c1));
            L_sen_real = sum(deltaf*beta1.*log2(1+P_s.*c2));
            E_buoy_real = sum(P_b.*alpha1+Pr.*beta1);
            E_sen_real = sum(P_s.*beta1);
            storage_BCD(l).Pb = P_b;
%}
      

 
%%%%%%%%%%%%%%%%%%%%%%�ܼ����ʷ
    EE_total_recod(l+1) = (sum(L_buoy_real)+sum(L_sen_real))/...
        (theta*sum(E_buoy_real)+(1-theta)*sum(E_sen_real));
    if (abs(EE_total_recod(l+1)-EE_total_recod(l))/EE_total_recod(l+1)>0.01)
        flagBCD=0;
        l=l+1;
    end
    clear q q_l q_l_SCA distance_l distance_l_SCA distance2 alpha1 beta1 P_b %P_s
 end
%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�������ֺ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�������ֺ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
function [q,distance2] = Qcvx(W,alpha1,rho_l,omega_l,distance_l_SCA,W_m,L_Ul,E_tbl,theta,...
                          E_rbl,E_tsl,L_buoy,Vmax,tau,T,M) 
cvx_clear
%cvx_solver SeDuMi
cvx_begin
            variable q(T,2)
            for t=1:T
                for m=1:M
                    distance2(t,m)=square_pos(norm([q(t,1) q(t,2)]-[W_m(1,m) W_m(2,m)]));
                end
            end
            L_A=sum(sum(W*alpha1.*(rho_l+omega_l.*(distance2-distance_l_SCA))));%͹+����
 %           minimize (-1*L_A)
             minimize (-1*(L_A+L_Ul)/(theta*(E_tbl+E_rbl)+(1-theta)*E_tsl));%ԭ����
                
            subject to
                 sum(W*alpha1.*(rho_l+omega_l.*(distance2-distance_l_SCA))) >= L_buoy;
            
            for t=1:T-1
                norm([q(t,1) q(t,2)]-[q(t+1,1) q(t+1,2)]) <= Vmax*tau;
            end
            q(1,1)==q(T,1);
            q(1,2)==q(T,2);
cvx_end            

ABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABABAB
function [alpha1,beta1] = ABcvx(W,deltaf,c1,c2,chi2,theta,...
                                Pr,L_buoy,L_sen,E_bub,E_bus,P_b,P_s,T,M)
cvx_solver SeDuMi
cvx_precision best
cvx_begin
            variable alpha1(T,M)
            variable beta1(T,M)
%           minimize(-sum(sum(W*alpha1.*log(1+P_b.*eps(h)/eps(sigma2))/0.6931))...
%                      -sum(sum(deltaf*beta1.*log(1+P_s.*c2)/0.6931))...
%                      +chi2*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1)));
            minimize(-sum(sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1)))/0.6931...
                     -sum(sum(-deltaf*beta1.*rel_entr(ones(T,M),1+P_s.*c2)))/0.6931...
                     +chi2*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1)));
            subject to
            
           L_buoy <= sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1))/0.6931;%sum()��������
           L_sen <= sum(-deltaf*beta1.*rel_entr(ones(T,M),1+P_s.*c2))/0.6931;
           
           sum(P_b.*alpha1+Pr.*beta1) <= E_bub;
           sum(beta1.*P_s) <= E_bus;
           
           sum(alpha1,2) <= 1;
            0 <= beta1;% <= 1;
%            beta1 <= 1;
           sum(beta1,2) <= 1;
           alpha1 >= 0;
cvx_end

PsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPsPs
function [P_s] = Pcvx_s(W,deltaf,alpha1,beta1,c1,c2,chi4,theta,...
                          Pr,L_buoy,L_sen,E_bub,E_bus,Pmax_s,Pmax_b,T,M,P_b)      
% cvx_clear                      
cvx_solver SeDuMi
cvx_precision best
cvx_begin
%            variable P_b(T,M)
            variable P_s(T,M)
%             minimize(-sum(sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1)))/0.6931...
%                      -sum(sum(-deltaf*beta1.*rel_entr(ones(T,M),1+P_s.*c2)))/0.6931...
%                      +chi1*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1))); %ԭ  
                 minimize(-sum(sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1)))/0.6931...
                     -sum(sum(-deltaf*beta1.*rel_entr(ones(T,M),1+P_s.*c2)))/0.6931...
                     +chi4*1*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1)));  
            subject to
%              L_buoy <= sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1))/0.6931;%sum()��������
             L_sen <= sum(-deltaf*beta1.*rel_entr(ones(T,M),1+P_s.*c2))/0.6931;
             
%              sum(P_b.*alpha1+Pr.*beta1) <= E_bub;
             sum(P_s.*beta1) <= E_bus;
              
              P_s >= 0;%<= Pmax_s;
              P_s <= Pmax_s;
%              P_b >= 0;%<= Pmax_b;
%              P_b <= Pmax_b;           
cvx_end

PbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPbPb
function [P_b] = Pcvx_b(W,deltaf,alpha1,beta1,c1,c2,chi1,theta,...
                          Pr,L_buoy,L_sen,E_bub,E_bus,Pmax_s,Pmax_b,T,M,P_s)      
% cvx_clear                      
cvx_solver SeDuMi
cvx_precision best
cvx_begin
            variable P_b(T,M)
%             variable P_s(T,M)
%             minimize(-sum(sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1)))/0.6931...
%                      -sum(sum(-deltaf*beta1.*rel_entr(ones(T,M),1+P_s.*c2)))/0.6931...
%                      +chi1*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1))); %ԭ  
                 minimize(-sum(sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1)))/0.6931...
                     -sum(sum(-deltaf*beta1.*rel_entr(ones(T,M),1+P_s.*c2)))/0.6931...
                     +chi1*1*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1)));  
            subject to
             L_buoy <= sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1))/0.6931;%sum()��������
%              L_sen <= sum(-deltaf*beta1.*rel_entr(ones(T,M),1+P_s.*c2))/0.6931;
             
             sum(P_b.*alpha1+Pr.*beta1) <= E_bub;
%              sum(P_s.*beta1) <= E_bus;
              
%               P_s >= 0;%<= Pmax_s;
%               P_s <= Pmax_s;
             P_b >= 0;%<= Pmax_b;
             P_b <= Pmax_b;           
cvx_end   