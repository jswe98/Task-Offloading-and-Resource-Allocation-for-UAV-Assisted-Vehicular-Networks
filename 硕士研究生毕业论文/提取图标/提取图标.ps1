Get-Icon -folder "F:\OneDrive\组内文献\硕士研究生毕业论文\提取图标"

Function Get-Icon {
    [CmdletBinding()]
    Param (
        [Parameter(Mandatory=$True, HelpMessage="输入.EXE文件的位置")]
        [string]$folder
    )

    [System.Reflection.Assembly]::LoadWithPartialName('System.Drawing') | Out-Null

    md $folder -ea 0 | Out-Null

    dir $folder *.exe -ea 0 -rec |
    ForEach-Object {
        $baseName = [System.IO.Path]::GetFileNameWithoutExtension($_.FullName)
        Write-Progress "正在提取图标" $basename
        [System.Drawing.Icon]::ExtractAssociatedIcon($_.FullName).ToBitmap().Save("$folder\$basename.ico")
    }
}

Get-Help Get-Icon -Full