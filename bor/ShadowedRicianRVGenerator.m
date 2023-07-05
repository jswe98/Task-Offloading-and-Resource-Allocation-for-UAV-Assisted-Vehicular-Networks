function [shadowedricianRV] = ShadowedRicianRVGenerator(Nunit,b,m,omega)%2b是散射分量平均功率，omega是LOS分量
%Nunit是信道数量，n是仿真次数

RayleighRV=sqrt(b)*(randn(Nunit,1)+1i*randn(Nunit,1));%Rayleigh distribution
NakagamiRV=sqrt(gamrnd(m,omega/m,Nunit,1));%Nakagami distribution
zeta=0;
shadowedricianRV=RayleighRV+NakagamiRV*exp(1i*zeta);%ShadowedRician distribution

%mu=m1*ones(1:length(n));
%omega=omega1*ones(1:length(n));