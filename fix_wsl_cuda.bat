@echo off
chcp 65001 >nul
echo ============================================
echo   WSL2 CUDA 修复脚本
echo ============================================
echo.

echo [1/4] 关闭所有 WSL 实例...
wsl --shutdown
timeout /t 3 /nobreak >nul
echo 完成.

echo [2/4] 将 Ubuntu 设为默认 WSL 发行版...
wsl --set-default Ubuntu
if %errorlevel% neq 0 (
    echo 错误: 无法设置默认发行版
    pause
    exit /b 1
)
echo 完成.

echo [3/4] 验证 Ubuntu 中的 CUDA...
echo.
wsl -d Ubuntu bash -c "nvidia-smi"
echo.

echo [4/4] 验证 PyTorch CUDA...
wsl -d Ubuntu bash -c "python3 -c 'import torch; print(\"CUDA available:\", torch.cuda.is_available()); print(\"PyTorch:\", torch.__version__); print(\"CUDA version:\", torch.version.cuda)'"
echo.

echo ============================================
echo   修复完成！
echo.
echo   现在可以使用以下方式启动项目：
echo     wsl -d Ubuntu
echo     cd /mnt/e/AI_Agents/protein_designer
echo     python run.py
echo.
echo   或者在 Windows 终端直接运行：
echo     wsl -d Ubuntu bash -c "cd /mnt/e/AI_Agents/protein_designer && python run.py"
echo ============================================
pause
