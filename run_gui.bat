@echo off
REM Simple launcher for CFD Solver GUI
REM Uses Python directly (no exe needed)

echo ============================================================
echo CFD Solver GUI - Launcher
echo ============================================================
echo.

REM Check if Python is available
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python not found!
    echo Please install Python or Anaconda
    pause
    exit /b 1
)

echo Starting CFD Solver GUI...
echo.

REM Launch the GUI
python cfd_gui.py

if errorlevel 1 (
    echo.
    echo ERROR: Failed to start GUI
    echo Check that required packages are installed:
    echo   pip install numpy matplotlib pillow
    pause
)
