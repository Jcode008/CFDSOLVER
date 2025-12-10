@echo off
REM CFD Solver GUI - Simple Python Launcher
REM This runs the GUI directly with Python (no exe building needed!)

echo ============================================================
echo          CFD Solver GUI - Starting...
echo ============================================================
echo.

REM Use Anaconda Python
set PYTHON=C:\Users\graha\anaconda3\python.exe

REM Check if Python exists
if not exist "%PYTHON%" (
    echo ERROR: Anaconda Python not found at %PYTHON%
    echo Please update the path in this batch file
    pause
    exit /b 1
)

echo Using: %PYTHON%
echo.
echo Starting GUI with live plotting capabilities...
echo.

REM Launch the GUI
"%PYTHON%" cfd_gui.py

if errorlevel 1 (
    echo.
    echo ERROR: Failed to start GUI
    echo Check the error message above
    pause
)
