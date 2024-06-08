@echo off  
chcp 65001  
setlocal EnableDelayedExpansion  
  
set "script_folder=%~dp0"  
set "parent_folder=!script_folder:~0,-1!"  
for %%i in ("!parent_folder!.") do set "parent_folder=%%~dpi"  
  
REM 获取父级目录的父级目录  
for %%i in ("!parent_folder!.") do set "grandparent_folder=%%~dpi"  
echo 当前脚本父级目录的父级目录是：%grandparent_folder%  
  
REM 删除旧的 PDF 文件  
del "%grandparent_folder%202121030188_魏建帅_学位论文.pdf"  
if exist "%grandparent_folder%202121030188_魏建帅_学位论文.pdf" (  
    echo 删除文件失败。  
    pause  
    exit /b  
)  
  
REM 复制新的 PDF 文件到相同位置  
copy "%parent_folder%202121030188_魏建帅_学位论文.pdf" "%grandparent_folder%202121030188_魏建帅_学位论文.pdf"  
if not exist "%grandparent_folder%202121030188_魏建帅_学位论文.pdf" (  
    echo 复制文件失败。  
    pause  
    exit /b  
)  
  
echo 文件已成功替换。  
exit