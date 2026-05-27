@echo off
chcp 65001 >nul
echo ============================================
echo   Protein Design Studio v0.1.0
echo ============================================
echo.
echo 在 WSL2 Ubuntu 中启动服务...
echo 浏览器打开: http://localhost:8899
echo.

wsl -d Ubuntu bash -c "cd /mnt/e/AI_Agents/protein_designer && python3 run.py"

pause
