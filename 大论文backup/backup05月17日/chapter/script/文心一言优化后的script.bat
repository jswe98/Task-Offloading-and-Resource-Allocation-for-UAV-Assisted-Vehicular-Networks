@echo off    
chcp 65001    
setlocal EnableDelayedExpansion    
  
REM 预先定义文件名和变量  
set "filename=202121030188_魏建帅_学位论文.pdf"    
set "authorname=魏建帅"        
set "StudentNumber=202121030188"       
set "Title=高动态异构车联网络资源联合鲁棒优化分配的研究"         
set "yuanfilename=template.pdf"      
set "newfilename=!StudentNumber!_!authorname!_!Title!.pdf"     
  
REM 获取父级目录  
set "script_folder=%~dp0"  
set "parent_folder=!script_folder:~0,-1!"  
for %%i in ("!parent_folder!.") do set "parent_folder=%%~dpi"  
  
REM 获取爷级目录  
for %%i in ("!parent_folder!.") do set "grandparent_folder=%%~dpi"  
echo 当前脚本的爷级目录是：%grandparent_folder%  
  
REM 删除旧的 PDF 文件  
if exist "%grandparent_folder%!filename!" (  
    del "%grandparent_folder%!filename!"  
    if exist "%grandparent_folder%!filename!" (  
        echo 删除文件失败。  
        pause  
        exit /b  
    )  
)  
  
REM 复制新的 PDF 文件到爷级目录  
if exist "%grandparent_folder%!yuanfilename!" (  
    copy "%grandparent_folder%!yuanfilename!" "%grandparent_folder%!newfilename!"  
    if not exist "%grandparent_folder%!newfilename!" (  
        echo 复制文件失败。  
        pause  
        exit /b  
    )  
) else (  
    echo 找不到模板文件：%grandparent_folder%!yuanfilename!  
    pause  
    exit /b  
)  
  
echo 文件已成功替换。  
exit /b