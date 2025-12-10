@echo off
REM Batch script to build standalone executable

echo ============================================================
echo CFD Solver GUI - Standalone Executable Builder
echo ============================================================
echo.

REM Check if PyInstaller is installed
python -c "import PyInstaller" 2>nul
if errorlevel 1 (
    echo PyInstaller not found. Installing...
    pip install pyinstaller
    if errorlevel 1 (
        echo ERROR: Failed to install PyInstaller
        pause
        exit /b 1
    )
)

echo PyInstaller is ready!
echo.
echo Building executable...
echo This will take a few minutes...
echo.

REM Run PyInstaller directly
pyinstaller cfd_gui.py ^
    --name=CFDSolverGUI ^
    --onefile ^
    --windowed ^
    --clean ^
    --noconfirm ^
    --icon=app_icon.ico ^
    --hidden-import=numpy ^
    --hidden-import=matplotlib ^
    --hidden-import=matplotlib.backends.backend_tkagg ^
    --hidden-import=PIL ^
    --add-data="analysis;analysis" ^
    --optimize=2

if errorlevel 1 (
    echo.
    echo ERROR: Build failed!
    pause
    exit /b 1
)

echo.
echo ============================================================
echo BUILD SUCCESSFUL!
echo ============================================================
echo.
echo Executable location: dist\CFDSolverGUI.exe
echo.
echo You can now:
echo   1. Run dist\CFDSolverGUI.exe
echo   2. Share it with others (no Python installation needed)
echo   3. The .exe is standalone but still needs CMake/C++ solver
echo.
echo ============================================================
pause
