 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5.21,9:57%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%全线贯通,%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%未能大循环%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% E2RA
%% April 28 beginning
%
%% 命名规则
%{
l(M,1)%传感器与浮标的距离
Wm %sink坐标 
cmar %连接矩阵，connective martrix
cmi,cmj:连接矩阵的计数变量
h(t,m) %无线电信道增益
c2(t,m) &水声链路信道增益
chi %丁克尔巴赫最优值
chi_rec %chi的record变量
%}
clear
%% 总初始化
M=5;%sink个数
T=90;%时隙个数
Ls=3000;%海域规模scale(m)
Vmax=80;
beta0=10^(-6);%单位信道增益
sigma2=10^(-14);%噪声功率
deltaf=1000;%水声带宽
%af=2.59;%f=20kHz
af=4.1;%f=25kHz
kappa=1.5;%常数
%Nf=10^(-14.2);%噪声功率(7文2图上的30纵轴)
Nf=10^(-12.5);%噪声功率(7文2图上的50纵轴)可能好一些
W=10^(6);%无线电带宽1MHz
%Tt=100;%飞行周期（s）
tau=1;%时隙长度（s）
H=100;%UAV飞行高度（m）
Pr=8*10^(-4)*ones(T,M);%浮标接收数据功率（W）
%% 系统可变初始化
theta=0.4;%能量价值权重
% L_buoy=300*1000*ones(1,M);%暂时让所有浮标的任务量相等&除了ps以外成功的参数
% L_sen=100*1000*ones(1,M);%暂时让所有浮标的任务量相等
L_buoy=300*1000*ones(1,M);%暂时让所有浮标的任务量相等
L_sen=100*1000*ones(1,M);%暂时让所有浮标的任务量相等
E_bub=0.02*ones(1,M);%buoy power budget功耗预算
E_bus=0.02*ones(1,M);%sensor power budget功耗预算
Pmax_b=0.1*ones(T,M);%buoy最大功率(W)
Pmax_s=0.1*ones(T,M);%sensor最大功率
L=1;%BCD总迭代次数
epsilon2=0.1;%收敛判定阈值
epsilon1=0.1;
%以上是瞎胡来的
%{
theta=0.4;%能量价值权重
L_buoy=2000*1000*ones(1,M);%暂时让所有浮标的任务量相等
L_sen=100*1000*ones(1,M);%暂时让所有浮标的任务量相等
E_bub=0.01*ones(1,M);%buoy power budget功耗预算
E_bus=0.005*ones(1,M);%sensor power budget功耗预算
Pmax_b=0.1*ones(T,M);%buoy最大功率(W)
Pmax_s=0.1*ones(T,M);%sensor最大功率
L=1;%BCD总迭代次数
epsilon=0.1;%收敛判定阈值
%}
%原来的
%% 节点分布
W_m = [Ls/22.2,   Ls/3.98,  Ls/1.95,  Ls/1.41,  Ls/2.3;
       Ls/2.2,    Ls/1.28, Ls/1.23,  Ls/2.9,   Ls/6.2];
   geocent(1,1) = sum(W_m(1,:))/M;%几何中心geometry centre
   geocent(1,2) = sum(W_m(2,:))/M;
  
for m=1:M % 浮标节点
    %{
%     W_m(1,m) = 800;
%     W_m(2,m) = 500;
    W_m(1,m) = Ls*rand;
    W_m(2,m) = Ls*rand;%sink坐标随机生成
    %}
    g2b(1,m)=sqrt((geocent(1,1)-W_m(1,m))^(2)+(geocent(1,2)-W_m(2,m))^(2));%geometry centre to buoys
    %画出sink位置
    %{
    scatter (W_m(1,m),W_m(2,m),'mo');
    hold on
    %}
    %{
    gplot (cmar,q,'g<-');%优化轨迹
    gplot (cmar,q_init,'g<-');
    %}
end
r_init = 0.8*max(g2b);
%% s2b 距离
%l_s2b = 0.1;
%l_s2b = 1+rand(1,M);%传感器节点到相应浮标的距离(水下是km)
%l_s2b = [1.3816, 1.7655, 1.7952, 1.1869, 1.4898];
l_s2b = [2.3816, 1.9655, 2.2952, 2.1869, 2.4898];
%l_s2b = [2.5816, 2.9655, 2.3952, 2.7869, 2.4898];
%{
for m=1:M %结构体
    sink(m).Wm_x = W_m(1,m);
    sink(m).Wm_y = W_m(2,m);%sink位置存入结构体
    sink(m).se_bu = l_s2b(1,m);%sink--sensor距离存入结构体
end
%}
%% BCD初始化
for m=1:M %初始化alpha&beta，均分
    for t=1:T
        alpha_init(t,m)=1/M;%浮标IT时间
        beta_init(t,m)=0.09;%传感器IT时间
    end
end
for t=1:T %初始化轨迹Q:q_init
    x_init = geocent(1,1)+r_init*cos(t*2*pi/T);
    y_init = geocent(1,2)+r_init*sin(t*2*pi/T);
    q_init(t,1)=x_init;%q的坐标必须是T*2形式
    q_init(t,2)=y_init;
end
cmar=zeros(T,T);%生成连接矩阵
for cmi=1:T
    for cmj=1:T
        if (abs(cmi-cmj)==1||abs(cmi-cmj)==T-1)
            cmar(cmi,cmj)=1;
        end
    end
end
%绘制初始轨迹
%{
gplot (cmar,q_init,'g<-');%初始轨迹，cmar是连接矩阵，connective martrix
gplot (cmar,q,'g<-');%优化轨迹
%}
alpha1(:,:) = alpha_init;
beta1(:,:) = beta_init;
P_b = 1*10^(-4)*ones(T,M);
P_s = 1*10^(-3)*ones(T,M);

% P_b = 10^(-4)*ones(T,M);%5.18,16:51原来的,不行马上改回，此时Q可以
% P_s = 10^(-3)*ones(T,M);
q_l=q_init;

%% BCD开始
 
flagBCD=0;
l=1;%BCD迭代轮数
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
    %%%%%%%n_3%%%%%%%%%%%%%%%%%%%%%接着P^s,P^b,A,B给定，求Q  
    %
    n_3=1;
    flag=0;
    q_l_SCA=q_l;
    while (flag==0)
        flag=1;
        for t = 1:T %Qcvx里所必需的定量值
            for m = 1:M
                c2(t,m) = 1/(l_s2b(1,m)^(kappa)*af^(l_s2b(1,m))*Nf*deltaf);%原
                c3_l(t,m)=P_b(t,m)/sigma2;
                omega_l(t,m)=(-1*c3_l(t,m)*beta0)/(((q_l_SCA(t,1)-W_m(1,m))^(2)+(q_l_SCA(t,2)-W_m(2,m))^(2)+H^(2))*log(2)*...
                    ((q_l_SCA(t,1)-W_m(1,m))^(2)+(q_l_SCA(t,2)-W_m(2,m))^(2)+H^(2)+c3_l(t,m)*beta0));
                rho_l(t,m)=log2(1+c3_l(t,m)*beta0/((q_l_SCA(t,1)-W_m(1,m))^(2)+(q_l_SCA(t,2)-W_m(2,m))^(2)+H^(2)));
                distance_l(t,m)=(q_l_SCA(t,1)-W_m(1,m))^(2)+(q_l_SCA(t,2)-W_m(2,m))^(2);
                %第l次迭代，第t个时隙UAV与第m个浮标水平距离平方
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
        q_l_SCA=q;%记录本轮优化结果，供下轮使用.潜在威胁？实在不行不输出distance2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%Q的御史大夫
        for t = 1:T %Qcvx里所必需的定量值
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
        ob(n_3) = (sum(L_buoy_real)+L_Ul)/(theta*(E_tbl+E_rbl)+(1-theta)*E_tsl);%43文Fig.2的eta
        ob_lb(n_3) = (sum(L_A_lb)+L_Ul)/(theta*(E_tbl+E_rbl)+(1-theta)*E_tsl);
        if (abs(sum(ob(n_3)-ob_lb(n_3)))/sum(ob_lb(n_3))>0.01)
            flag=0;
        end
        n_3=n_3+1;
    end %对应while
    storage_BCD(l).Q = q_l_SCA;
    q_l = storage_BCD(l).Q;
    %}
        
        
        
       
    
            %%%%%%%%%%%%%n_2%%%%%%chi2%%%%%%%%%%%%其次P^s,P^b,Q给定，求A,B(linear program)
    %
    
    for t = 1:T %给出信道增益h(t,m)和c2(t,m)
        for m = 1:M
            h(t,m) = beta0/(H^(2)+(q_l(t,1)-W_m(1,m))^(2)+(q_l(t,2)-W_m(2,m))^(2));
            c1(t,m) = h(t,m)/sigma2;
            c2(t,m) = 1/(l_s2b(1,m)^(kappa)*af^(l_s2b(1,m))*Nf*deltaf);%原
        end
    end
    %
    n_2=1; %求解A,B的dinkelbach轮次
    flag=0;
    chi2=0;    
    while (flag==0)%Dinkelbach
        flag=1;
        chi2_rec(n_2)=chi2;
        [alpha1,beta1]=ABcvx(W,deltaf,c1,c2,chi2,theta,...
                                Pr,L_buoy,L_sen,E_bub,E_bus,P_b,P_s,T,M);%调用CVX
        F_chi2=sum(sum(-W*alpha1.*rel_entr(1,1+P_b.*c1)))/0.6931...
                     +sum(sum(-deltaf*beta1.*rel_entr(1,1+P_s.*c2)))/0.6931...
                     -chi2*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1));%更新F_chi
        chi2=(sum(sum(-1*W*alpha1.*rel_entr(1,1+P_b.*c1)))/0.6931...
                     +sum(sum(-1*deltaf*beta1.*rel_entr(1,1+P_s.*c2)))/0.6931)...
                     /sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1));
         n_2=n_2+1;
         if (F_chi2>epsilon2)
             flag=0;
         end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%A&B的御史大夫
            L_buoy_real = sum(W*alpha1.*log2(1+P_b.*c1));
            L_sen_real = sum(deltaf*beta1.*log2(1+P_s.*c2));
            E_buoy_real = sum(P_b.*alpha1+Pr.*beta1);
            E_sen_real = sum(P_s.*beta1);
            storage_BCD(l).A = alpha1;
            storage_BCD(l).B = beta1;
    %}
    

    
    %%%%%%%%%%chi1%%%%%F_chi1%%%%%%%%%%%%%%%%%%%%%首先A,B,Q给定，求P^s,P^b
    %%%%%%%%%%%chi1%%%%%%%%%%%%%%%%%%%%%%%%%%首先A,B,Q给定，求P^b
    %
    for t = 1:T %给出信道增益h(t,m)和c2(t,m)
        for m = 1:M
            h(t,m) = beta0/(H^(2)+(q_l(t,1)-W_m(1,m))^(2)+(q_l(t,2)-W_m(2,m))^(2));
            c1(t,m) = h(t,m)/sigma2;
            c2(t,m) = 1/(l_s2b(1,m)^(kappa)*af^(l_s2b(1,m))*Nf*deltaf);%原
        end
    end
    %
    n_1=1;%求解功率的dinkelbach轮次
    flag=0;
    chi1=0;
    while (flag==0)%Dinkelbach
        flag=1;
        chi1_rec(n_1)=chi1;%chi1的record变量

        
        [P_b] = Pcvx_b(W,deltaf,alpha1,beta1,c1,c2,chi1,theta,...
                          Pr,L_buoy,L_sen,E_bub,E_bus,Pmax_s,Pmax_b,T,M,P_s);
%         [P_b,P_s] = Pcvx(W,deltaf,alpha1,beta1,c1,c2,chi1,theta,...
%                          Pr,L_buoy,L_sen,E_bub,E_bus,Pmax_s,Pmax_b,T,M);%调用CVX，原来联合
                      
         F_chi1=sum(sum(-W*alpha1.*rel_entr(1,1+P_b.*c1)))/0.6931...
                     +sum(sum(-deltaf*beta1.*rel_entr(1,1+P_s.*c2)))/0.6931...
                     -chi1*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1));%更新F_chi
         chi1=(sum(sum(-1*W*alpha1.*rel_entr(1,1+P_b.*c1)))/0.6931...
                     +sum(sum(-1*deltaf*beta1.*rel_entr(1,1+P_s.*c2)))/0.6931)...
                     /sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1));
        n_1=n_1+1;
         if (F_chi1>epsilon1)
             flag=0;
         end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%P的御史大夫
            for t = 1:T %Qcvx里所必需的定量值
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
      

 
%%%%%%%%%%%%%%%%%%%%%%总监察御史
    EE_total_recod(l+1) = (sum(L_buoy_real)+sum(L_sen_real))/...
        (theta*sum(E_buoy_real)+(1-theta)*sum(E_sen_real));
    if (abs(EE_total_recod(l+1)-EE_total_recod(l))/EE_total_recod(l+1)>0.01)
        flagBCD=0;
        l=l+1;
    end
    clear q q_l q_l_SCA distance_l distance_l_SCA distance2 alpha1 beta1 P_b %P_s
 end
%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%各个部分函数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%各个部分函数%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            L_A=sum(sum(W*alpha1.*(rho_l+omega_l.*(distance2-distance_l_SCA))));%凸+凹？
 %           minimize (-1*L_A)
             minimize (-1*(L_A+L_Ul)/(theta*(E_tbl+E_rbl)+(1-theta)*E_tsl));%原来的
                
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
            
           L_buoy <= sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1))/0.6931;%sum()是行向量
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
%                      +chi1*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1))); %原  
                 minimize(-sum(sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1)))/0.6931...
                     -sum(sum(-deltaf*beta1.*rel_entr(ones(T,M),1+P_s.*c2)))/0.6931...
                     +chi4*1*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1)));  
            subject to
%              L_buoy <= sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1))/0.6931;%sum()是行向量
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
%                      +chi1*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1))); %原  
                 minimize(-sum(sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1)))/0.6931...
                     -sum(sum(-deltaf*beta1.*rel_entr(ones(T,M),1+P_s.*c2)))/0.6931...
                     +chi1*1*sum(sum(theta*(P_b.*alpha1+Pr.*beta1)+(1-theta)*P_s.*beta1)));  
            subject to
             L_buoy <= sum(-W*alpha1.*rel_entr(ones(T,M),1+P_b.*c1))/0.6931;%sum()是行向量
%              L_sen <= sum(-deltaf*beta1.*rel_entr(ones(T,M),1+P_s.*c2))/0.6931;
             
             sum(P_b.*alpha1+Pr.*beta1) <= E_bub;
%              sum(P_s.*beta1) <= E_bus;
              
%               P_s >= 0;%<= Pmax_s;
%               P_s <= Pmax_s;
             P_b >= 0;%<= Pmax_b;
             P_b <= Pmax_b;           
cvx_end   