% ���� data ����Ҫ���������  
data = rand(10,1);  
  
% ���ļ���д������  
fileID = fopen('temp.txt','w');  
  
% ������д���ļ�  
fprintf(fileID, '%f\n', data);  
  
% �ر��ļ�  
fclose(fileID);

% ���ļ��Զ�ȡ����  
fileID = fopen('temp.txt','r');  
  
% ���ļ��ж�ȡ����  
data = fscanf(fileID, '%f\n');  
  
% �ر��ļ�  
fclose(fileID);