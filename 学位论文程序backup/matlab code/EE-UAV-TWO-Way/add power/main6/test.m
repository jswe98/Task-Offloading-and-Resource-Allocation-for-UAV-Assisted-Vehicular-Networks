clc
close all
clear
lright=1;%�������������� right
lleft=-1;%��������������
roadlength=800;
roadwidth=30;
Vvelocity=20;
T=70;
det=1;
M=9;
while 1
for m = 1:M
random_number = randi([0, 1]); % ����0��1֮������������0����-1��1����1��
if random_number == 0
    ll(m) = lright;   % �����ܵĳ�
    random_number = randi([0, 1]);
    if random_number == 0
        Vposition_x=0;   
    else
        Vposition_x=roadlength/2;
    end   
    random_numbery = (roadwidth/2) * rand;  % ����4��8֮��������
else
    ll(m) = lleft; %�����ܵĳ�
    random_number = randi([0, 1]);
    if random_number == 0
        Vposition_x=roadlength/2;
    else
        Vposition_x=roadlength/4;
    end 
    random_numbery = roadwidth/2+(roadwidth/2) * rand; % ����0��4֮��������
end
%      Vposition(m,:,1)=[Vposition_x; roadwidth/2+roadwidth/2*ll(m); 0];
     Vposition(m,:,1)=[Vposition_x; random_numbery+roadwidth/2*ll(m); 0];
     for t=2:T
%      Vposition(m,:,t)=[Vposition(m,1,t-1)+ll(m)*det*Vvelocity; roadwidth/2+roadwidth/2*ll(m); 0];
     Vposition(m,:,t)=[Vposition(m,1,t-1)+ll(m)*det*Vvelocity; Vposition(m,2); 0];
     end
end
      if  sum(ll == 1)==0.5*M;
          break; 
          elseif sum(ll == 1)==0.5*M+0.5;
          break; 
          elseif sum(ll == 1)==0.5*M-0.5;
          break;    
      end 
end
