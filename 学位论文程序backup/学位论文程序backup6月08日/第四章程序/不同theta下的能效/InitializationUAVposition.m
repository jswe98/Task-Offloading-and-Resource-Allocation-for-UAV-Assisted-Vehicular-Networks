function [UAVposition] = InitializationUAVposition(T,roadlength,roadwidth)
Qv=zeros(T,2);%UAV³õÊ¼»¯×ø±ê
Qv(1,:)=[0,0];
Qv(T,:)=[roadlength,roadwidth];
for n=2:T-1
    Qv(n,:)=Qv(n-1,:)+[roadlength/(T-1),roadwidth/(T-1)];
end

UAVposition=Qv;
