clc
close all
clear
global N
N=300;%迭代次数
f5 = 7; f4 = f5; f3 = f4; f2 = f3; f1 = f2;
f_ba=5;c_exe=2;Dmax=0.2;f_total=30;
BR=[45.8974423498470,46.6472266809215,46.2453835924562,46.7441997856849,46.4685553688578];
% BRR=[6177.64250125427,6173.29762034172,6175.56932824611,6172.76796694858,6174.29210005848];
D1=[0.0115129097838468,0.0109235695934443,0.0112317087821159,0.0108517249037377,0.0110584640913346];
for u=1:N
    fi(1,:)=[f1,f2,f3,f4,f5];
    Ksai(1,:)=linspace(3000,3000,5);
    Fai(1,:)=linspace(15,15,5);
    for i=1:5
%         RRR=(exp(P(k,:))*G(i,i))/(I1(i)+Delta);
%         Ri=(log(C1+RRR))/log(2);
%         D1=C1./(tau*Ri-varsigma);
%         BR=Ri.*c_exe/(t_imax*d_up); 
        BRR=BR+c_exe.*Ksai(u,:);
        sumfai(u,i)=sum(Fai(u,:));
        fi(u+1,:)=max(0,sqrt((BRR)/sumfai(u,i))-f_ba);  
    if  max(abs(fi(u+1,:)-fi(u,:)))<=1e-6  %该数值可改，表示乘子在一定的程度下就认为已经稳定了
        fi(u+1,:)=fi(u,:);
    end
    end    
   f1= fi(u+1,1);
   f2= fi(u+1,2);
   f3= fi(u+1,3);
   f4= fi(u+1,4);
   f5= fi(u+1,5);
   sumfi(u,1)=sum(fi(u,:));       

    for i=1:5       
        Ksai(u+1,i)=Ksai(u,i)+1/u*(Dmax-c_exe/(fi(u,i)+f_ba)+D1(1,i));
    end    
    for i=1:5       
        sumfi(u,i)=sum(fi(u,:));
        Fai(u+1,i)=max(0,Fai(u,i)+1/u*(sumfi(u,i)-f_total));
    end  
    if  max(abs(Fai(u+1,:)-Fai(u,:)))<=1e-6  %该数值可改，表示乘子在一定的程度下就认为已经稳定了
        Fai(u+1,:)=Fai(u,:);
    end
end    