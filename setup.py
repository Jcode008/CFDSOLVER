"""
Setup script to create standalone executable for CFD Solver GUI
Uses PyInstaller to bundle everything into a single .exe file

Usage: python setup.py
"""

from PyInstaller.__main__ import run
import os
import sys

# Application info
app_name = "CFDSolverGUI"
script_file = "cfd_gui.py"
icon_file = "app_icon.ico"  # Optional, will create if needed

# PyInstaller configuration
args = [
    script_file,
    '--name=%s' % app_name,
    '--onefile',  # Single executable file
    '--windowed',  # No console window (GUI only)
    '--clean',
    '--noconfirm',
    
    # Include hidden imports
    '--hidden-import=numpy',
    '--hidden-import=matplotlib',
    '--hidden-import=matplotlib.backends.backend_tkagg',
    '--hidden-import=PIL',
    '--hidden-import=tkinter',
    
    # Include data files (analysis scripts)
    '--add-data=analysis;analysis',
    
    # Optimization
    '--optimize=2',
    
    # Add icon if it exists
]

# Add icon if available
if os.path.exists(icon_file):
    args.append('--icon=%s' % icon_file)

print("=" * 60)
print("Building standalone executable for CFD Solver GUI")
print("=" * 60)
print(f"\nApplication: {app_name}")
print(f"Script: {script_file}")
print("\nThis will create:")
print(f"  - dist/{app_name}.exe (standalone executable)")
print("\nBuilding... (this may take a few minutes)")
print("=" * 60 + "\n")

# Run PyInstaller
run(args)

print("\n" + "=" * 60)
print("BUILD COMPLETE!")
print("=" * 60)
print(f"\n✓ Executable created: dist/{app_name}.exe")
print("\nYou can now:")
print(f"  1. Run dist/{app_name}.exe directly (no Python needed)")
print("  2. Share the .exe file with others")
print("  3. The executable is ~100-200 MB (includes Python + all libraries)")
print("\nNote: The .exe still needs the C++ solver (CMake build) to function.")
print("=" * 60)
