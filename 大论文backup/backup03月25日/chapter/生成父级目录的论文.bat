@echo off    
chcp 65001    
setlocal EnableDelayedExpansion    
set "script_folder=%~dp0"    
set "parent_folder=!script_folder:~0,-1!"    
for %%i in ("!parent_folder!.") do set "parent_folder=%%~dpi"    
echo 当前脚本父级目录是：%parent_folder%    
  
REM 预先定义文件名    
set "filename=学位论文.pdf"    
set "authorname=张三"    
set "StudentNumber=202121030188"   
set "Title=无人机辅助的车联网络任务卸载与资源分配的研究"     
set "yuanfilename=template.pdf"   
  
REM 删除旧的 PDF 文件    
del "%parent_folder%!filename!"    
if exist "%parent_folder%!filename!" (    
    echo 删除文件失败。    
    pause    
    exit /b    
)    
  
REM 复制新的 PDF 文件到相同位置    
copy "%parent_folder%!yuanfilename!" "%parent_folder%!StudentNumber!_!Title!.pdf"    
if not exist "%parent_folder%!StudentNumber!_!Title!.pdf" (    
    echo 复制文件失败。    
    pause    
    exit /b    
)    
  
echo 文件已成功替换。    
exit