% ----------构造贝恩斯坦近似干扰约束的系数矩阵----------% 
function [Q] = rate(V1,p)
 rate(V1,p)
    I=H;%贝恩斯坦近似约束矩阵
    for i=1:5
      I(i,i)=0;%构造干扰约束中的系数矩阵
    end
    
global N
N =25;
arfa=zeros(N,5);
beta=zeros(N,5);
SINR=zeros(N,5);
Pmax=log(0.5)*ones(N,5);
B=1;
C1=linspace(1,1,5);
GF=linspace(1.5,1.5,5);%功放系数
PC=0.05;%电路中固定消耗的功率
 S=zeros(N,5);
 s=zeros(N,5);
 Gama=0.1;%原来是
 Ith=2e-6;%干扰阈值Gama要给相应的存储空间
Temp=zeros(40,5);%每次迭代收敛的最优功率暂存区 40应为Z
Z=1;

while 1
   P=zeros(N,5);
P(1,:)=[-8,-8,-8,-8,-8]; 
Lamda=zeros(N,5);
Lamda(1,:)=linspace(20,20,5);
for k=1:N
  for i=1:5
     %  w=exp(P(k,:)')H;观测量
        I1(i)=exp(P(k,:))*G(:,i)+Delta-G(i,i)*exp(P(k,i));
        SINR(k,i)=G(i,i)*exp(P(k,i))/I1(i); 
        
    arfa(k,i)=SINR(k,i)/(1+SINR(k,i));
    Total_T(k)=B*dot(log(C1+SINR(k,:))/log(2),C1);
  end
       for i=1:5
          for j=1:5
              if j~=i
                S(k,j)=B/log(2)*arfa(k,j)*G(i,j)*SINR(k,j)/((exp(P(k,j))*G(j,j)))+Lamda(k,j)*I(i,j);%%到底是G(i,j)还是G(j,i),应该是G(i,j)。当j不变是第j列，代表的是各个CM向某个簇头传输的信息。
              else
             S(k,j)=0;
              end
          end
          s(k,i)=sum(S(k,:))+Gama(Z)*GF(i);
  
  end
%  Total_T(k+1)=B*dot(log(C1+SINR(k,:))/log(2),C1); 
P(k+1,:)=min(log(B.*arfa(k,:)/log(2))-log(s(k,:)),Pmax(k,:));
   if max(abs(exp(P(k+1,:))-exp(P(k,:))))<=1e-5
  P(k+1,:)=P(k,:);
   end
    for i=1:5
   Lamda(k+1,i)=max(0,Lamda(k,i)+1000000/k*(exp(P(k,:))*I(:,i)-Ith));
  w(k,i)=exp(P(k,:))*I(:,i)-Ith;
    end 
    if  max(abs(w(k,:)))<=3e-7    %该数值可改，表示乘子在一定的程度下就认为已经稳定了
   Lamda(k+1,:)=Lamda(k,:);
    end

end
D(Z)=Total_T(N);
 Temp(Z,:)=P(N,:);
 error =abs(D(Z)-Gama(Z)*(GF*exp(Temp(Z,:)')+PC));
  Gama(Z+1)=D(Z)/(GF*exp(Temp(Z,:)')+PC);
      if  error<=1e-3
    L(F,:)=exp(Temp(Z,:));
    Q(F)=D(Z);
    A(F)=Gama(Z+1);
           break;
%       else
%        Lamda=zeros(N,5);
%     Lamda(1,:)=linspace(4,6,5);
%         P=zeros(N,5);
%         P =Temp(Z,:); 
      end 
     Z=Z+1;
end

  end