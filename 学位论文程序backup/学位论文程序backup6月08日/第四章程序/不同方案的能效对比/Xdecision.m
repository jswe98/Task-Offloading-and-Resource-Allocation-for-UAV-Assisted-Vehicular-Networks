function [X,Y] = Xdecision(Pm,GVU,GVR,T,M)
X=zeros(M,T);
Pm=Pm';
for t=1:T
    for m=1:M
    if  Pm(m,t)*GVU(m,t)<Pm(m,t)*GVR(m,t) %车车向基站通信为1
        X(m,t)=1;
    else
        X(m,t)=0;
    end;
    end
end
Y=1-X;