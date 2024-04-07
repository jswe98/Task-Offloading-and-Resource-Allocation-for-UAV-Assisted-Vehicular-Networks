% 假设 data 是你要保存的数据  
data = rand(10,1);  
  
% 打开文件以写入数据  
fileID = fopen('temp.txt','w');  
  
% 将数据写入文件  
fprintf(fileID, '%f\n', data);  
  
% 关闭文件  
fclose(fileID);

% 打开文件以读取数据  
fileID = fopen('temp.txt','r');  
  
% 从文件中读取数据  
data = fscanf(fileID, '%f\n');  
  
% 关闭文件  
fclose(fileID);