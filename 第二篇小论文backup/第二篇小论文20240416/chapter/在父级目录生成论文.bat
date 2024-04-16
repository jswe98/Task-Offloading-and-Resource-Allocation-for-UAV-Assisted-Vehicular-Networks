@echo off    
chcp 65001    
setlocal EnableDelayedExpansion    
set "script_folder=%~dp0"    
set "parent_folder=!script_folder:~0,-1!"    
for %%i in ("!parent_folder!.") do set "parent_folder=%%~dpi"    
echo 当前脚本父级目录是：%parent_folder%    
    
REM 预先定义文件名和变量    
set "filename=论文.pdf"    
set "authorname=张三"    
set "StudentNumber=1"    
set "Title=Energy-Efficient for UAV-Assisted Vehicular Networks Under the Two-ways Street: Joint UAV Trajectory Optimization and Robust Power Control Approach"    
set "yuanfilename=Energy_Efficient_for_UAV_Assisted_Vehicular.pdf"  
set "Title=!Title::=_!"    
set "newfilename=!Title: =_!.pdf"    
  
REM 删除旧的 PDF 文件    
if exist "%parent_folder%!filename!" (    
    del "%parent_folder%!filename!"    
    if exist "%parent_folder%!filename!" (    
        echo 删除文件失败。    
        pause    
        exit /b    
    )    
)    
    
REM 复制新的 PDF 文件到相同位置    
if exist "%script_folder%!yuanfilename!" (    
    REM 注意这里，我们使用双引号将整个文件名括起来，以确保冒号不会导致问题  
    copy "%script_folder%!yuanfilename!" "!parent_folder!\!newfilename!"  
    if not exist "!parent_folder!\!newfilename!" (    
        echo 复制文件失败。    
        pause    
        exit /b    
    )    
) else (    
    echo 找不到模板文件：%parent_folder%!yuanfilename!    
    pause    
    exit /b    
)    
    
echo 文件已成功替换。    
exit /b