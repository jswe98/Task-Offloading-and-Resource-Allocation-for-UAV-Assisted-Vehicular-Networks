clc
close all
clear
epsilon2=9;
flag=0;
Pm=0.3;
F_chi2=90;
z=1;
while (flag==0)%Dinkelbach
        flag=1;
        F_chi2=F_chi2-1;
        z=z+1;
         if (F_chi2>epsilon2)
             flag=0;
         end
         Pm=Pm+1;
end