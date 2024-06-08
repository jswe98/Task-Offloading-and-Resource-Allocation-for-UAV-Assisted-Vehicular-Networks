@echo off  
chcp 65001  
  
REM 删除旧的 PDF 文件  
del "..\论文.pdf"  
if exist "..\论文.pdf" (  
    echo 删除文件失败。  
    pause  
    exit /b  
)  
  
REM 复制新的 PDF 文件到相同位置  
copy ".\template.pdf" "..\论文.pdf"  
if not exist "..\论文.pdf" (  
    echo 复制文件失败。  
    pause  
    exit /b  
)  
  
echo 文件已成功替换。  
exit