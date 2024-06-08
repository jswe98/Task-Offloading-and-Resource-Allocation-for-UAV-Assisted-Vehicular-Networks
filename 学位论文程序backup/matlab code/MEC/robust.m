
M=5;%BS个数
BS=zeros(2,M);%BS坐标

clc;
clear all;
A=[0 0;500 0;0 500;500 500;0 1000];
L1=sum(A.*A,2);
L2=A*A';
D=bsxfun(@plus,L1,L1')-2*L2;
D=sqrt(D);
K=zeros(1,5);
dsum=zeros(1,5);
    for i=1:5
    dsum(i)= sum(D(:,i));  
    end   
    
    
      for i=1:5
     K=D( i , : );
    end
    
   K=D( 1 , : );
   
   
BS=zeros(5,5);   
C=[1304,2312; 3639,1315; 4177,2244; 3712,1399; 3488,1535;];
for i=1:5
for j=1:5
BS(i,j)=sqrt((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2);
end
end
for i=1:5
C(i,i)=0;
end
%     C=[0,0; 500,0; 0,500; 500,500; 0,1000; ];  保存的数据
clc;
clear all;
clear
BS=zeros(5,5); 
C=[0,0; 500,0; 0,500; 500,500; 0,1000; ];
for i=1:5
for j=1:5
BS(i,j)=sqrt((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2);
end
end
for i=1:5
BS(i,i)=0;
end
   