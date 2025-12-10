"""
CFD Solver GUI Application
Interactive interface for running airfoil flow simulations
"""
import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
import subprocess
import threading
import os
import sys
from pathlib import Path
import webbrowser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.animation as animation
import time

class CFDSolverGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("CFD Airfoil Solver - Interactive GUI")
        self.root.geometry("900x800")
        
        # Set up paths
        self.base_dir = Path(__file__).parent
        self.build_dir = self.base_dir / "build" / "Release"
        self.src_dir = self.base_dir / "src"
        self.analysis_dir = self.base_dir / "analysis"
        
        self.running = False
        self.process = None
        
        # Live plotting
        self.live_plot_enabled = tk.BooleanVar(value=True)
        self.plot_window = None
        self.latest_snapshot = -1
        
        # Create GUI
        self.create_widgets()
        
    def create_widgets(self):
        # Main container with padding
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Title
        title = ttk.Label(main_frame, text="CFD Airfoil Flow Simulator", 
                         font=('Arial', 16, 'bold'))
        title.grid(row=0, column=0, columnspan=2, pady=10)
        
        # Create notebook for tabs
        notebook = ttk.Notebook(main_frame)
        notebook.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Tab 1: Simulation Parameters
        param_frame = ttk.Frame(notebook, padding="10")
        notebook.add(param_frame, text="Simulation Setup")
        
        # Tab 2: Run & Monitor
        run_frame = ttk.Frame(notebook, padding="10")
        notebook.add(run_frame, text="Run & Monitor")
        
        # Tab 3: Results
        results_frame = ttk.Frame(notebook, padding="10")
        notebook.add(results_frame, text="Results & Analysis")
        
        # === TAB 1: PARAMETERS ===
        self.create_parameter_inputs(param_frame)
        
        # === TAB 2: RUN & MONITOR ===
        self.create_run_interface(run_frame)
        
        # === TAB 3: RESULTS ===
        self.create_results_interface(results_frame)
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(1, weight=1)
        
    def create_parameter_inputs(self, parent):
        """Create input fields for simulation parameters"""
        
        # Flow Conditions Section
        flow_frame = ttk.LabelFrame(parent, text="Flow Conditions", padding="10")
        flow_frame.grid(row=0, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        
        row = 0
        
        # Velocity
        ttk.Label(flow_frame, text="Freestream Velocity (m/s):").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.velocity_var = tk.StringVar(value="65.0")
        ttk.Entry(flow_frame, textvariable=self.velocity_var, width=15).grid(row=row, column=1, sticky=tk.W, padx=5)
        ttk.Label(flow_frame, text="(typical: 50-100 m/s)", font=('Arial', 8)).grid(row=row, column=2, sticky=tk.W)
        row += 1
        
        # Angle of attack
        ttk.Label(flow_frame, text="Angle of Attack (degrees):").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.alpha_var = tk.StringVar(value="0.0")
        ttk.Entry(flow_frame, textvariable=self.alpha_var, width=15).grid(row=row, column=1, sticky=tk.W, padx=5)
        ttk.Label(flow_frame, text="(positive = nose up)", font=('Arial', 8)).grid(row=row, column=2, sticky=tk.W)
        row += 1
        
        # Density
        ttk.Label(flow_frame, text="Density (kg/m³):").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.density_var = tk.StringVar(value="1.225")
        ttk.Entry(flow_frame, textvariable=self.density_var, width=15).grid(row=row, column=1, sticky=tk.W, padx=5)
        ttk.Label(flow_frame, text="(1.225 = sea level air)", font=('Arial', 8)).grid(row=row, column=2, sticky=tk.W)
        row += 1
        
        # Viscosity
        ttk.Label(flow_frame, text="Kinematic Viscosity (m²/s):").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.viscosity_var = tk.StringVar(value="1.0e-5")
        ttk.Entry(flow_frame, textvariable=self.viscosity_var, width=15).grid(row=row, column=1, sticky=tk.W, padx=5)
        ttk.Label(flow_frame, text="(1.48e-5 = air at 15°C)", font=('Arial', 8)).grid(row=row, column=2, sticky=tk.W)
        row += 1
        
        # Airfoil Section
        airfoil_frame = ttk.LabelFrame(parent, text="Airfoil Geometry (NACA 4-digit)", padding="10")
        airfoil_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        
        row = 0
        
        # Chord
        ttk.Label(airfoil_frame, text="Chord Length (m):").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.chord_var = tk.StringVar(value="0.6")
        ttk.Entry(airfoil_frame, textvariable=self.chord_var, width=15).grid(row=row, column=1, sticky=tk.W, padx=5)
        row += 1
        
        # NACA digits
        ttk.Label(airfoil_frame, text="Max Camber (fraction):").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.camber_var = tk.StringVar(value="0.02")
        ttk.Entry(airfoil_frame, textvariable=self.camber_var, width=15).grid(row=row, column=1, sticky=tk.W, padx=5)
        ttk.Label(airfoil_frame, text="(0.02 = 2% for NACA 2xxx)", font=('Arial', 8)).grid(row=row, column=2, sticky=tk.W)
        row += 1
        
        ttk.Label(airfoil_frame, text="Camber Position (fraction):").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.camber_pos_var = tk.StringVar(value="0.4")
        ttk.Entry(airfoil_frame, textvariable=self.camber_pos_var, width=15).grid(row=row, column=1, sticky=tk.W, padx=5)
        ttk.Label(airfoil_frame, text="(0.4 = 40% for NACA x4xx)", font=('Arial', 8)).grid(row=row, column=2, sticky=tk.W)
        row += 1
        
        ttk.Label(airfoil_frame, text="Max Thickness (fraction):").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.thickness_var = tk.StringVar(value="0.12")
        ttk.Entry(airfoil_frame, textvariable=self.thickness_var, width=15).grid(row=row, column=1, sticky=tk.W, padx=5)
        ttk.Label(airfoil_frame, text="(0.12 = 12% for NACA xx12)", font=('Arial', 8)).grid(row=row, column=2, sticky=tk.W)
        row += 1
        
        # Simulation Section
        sim_frame = ttk.LabelFrame(parent, text="Simulation Settings", padding="10")
        sim_frame.grid(row=2, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        
        row = 0
        
        # Grid size
        ttk.Label(sim_frame, text="Grid Points (nx × ny):").grid(row=row, column=0, sticky=tk.W, pady=2)
        grid_subframe = ttk.Frame(sim_frame)
        grid_subframe.grid(row=row, column=1, sticky=tk.W, padx=5)
        self.nx_var = tk.StringVar(value="800")
        ttk.Entry(grid_subframe, textvariable=self.nx_var, width=7).pack(side=tk.LEFT)
        ttk.Label(grid_subframe, text=" × ").pack(side=tk.LEFT)
        self.ny_var = tk.StringVar(value="400")
        ttk.Entry(grid_subframe, textvariable=self.ny_var, width=7).pack(side=tk.LEFT)
        ttk.Label(sim_frame, text="(higher = better resolution, slower)", font=('Arial', 8)).grid(row=row, column=2, sticky=tk.W)
        row += 1
        
        # Domain size
        ttk.Label(sim_frame, text="Domain Size (Lx × Ly, m):").grid(row=row, column=0, sticky=tk.W, pady=2)
        domain_subframe = ttk.Frame(sim_frame)
        domain_subframe.grid(row=row, column=1, sticky=tk.W, padx=5)
        self.lx_var = tk.StringVar(value="4.0")
        ttk.Entry(domain_subframe, textvariable=self.lx_var, width=7).pack(side=tk.LEFT)
        ttk.Label(domain_subframe, text=" × ").pack(side=tk.LEFT)
        self.ly_var = tk.StringVar(value="2.0")
        ttk.Entry(domain_subframe, textvariable=self.ly_var, width=7).pack(side=tk.LEFT)
        row += 1
        
        # Timesteps
        ttk.Label(sim_frame, text="Number of Timesteps:").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.nt_var = tk.StringVar(value="5000")
        ttk.Entry(sim_frame, textvariable=self.nt_var, width=15).grid(row=row, column=1, sticky=tk.W, padx=5)
        ttk.Label(sim_frame, text="(5000 = ~1 hour)", font=('Arial', 8)).grid(row=row, column=2, sticky=tk.W)
        row += 1
        
        # Snapshot interval
        ttk.Label(sim_frame, text="Snapshot Interval:").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.snapshot_var = tk.StringVar(value="50")
        ttk.Entry(sim_frame, textvariable=self.snapshot_var, width=15).grid(row=row, column=1, sticky=tk.W, padx=5)
        ttk.Label(sim_frame, text="(save every N steps)", font=('Arial', 8)).grid(row=row, column=2, sticky=tk.W)
        row += 1
        
        # Reynolds number display
        self.reynolds_label = ttk.Label(sim_frame, text="Reynolds Number: (calculate)", 
                                       font=('Arial', 10, 'bold'), foreground='blue')
        self.reynolds_label.grid(row=row, column=0, columnspan=3, pady=10)
        
        # Calculate button
        calc_btn = ttk.Button(sim_frame, text="Calculate Reynolds Number", 
                             command=self.calculate_reynolds)
        calc_btn.grid(row=row+1, column=0, columnspan=3, pady=5)
        
        # Apply parameters button
        apply_btn = ttk.Button(parent, text="✓ Apply Parameters & Update Code", 
                              command=self.apply_parameters, style='Accent.TButton')
        apply_btn.grid(row=3, column=0, columnspan=2, pady=20)
        
    def create_run_interface(self, parent):
        """Create run and monitoring interface"""
        
        # Status frame
        status_frame = ttk.LabelFrame(parent, text="Simulation Status", padding="10")
        status_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        status_frame.columnconfigure(0, weight=1)
        status_frame.rowconfigure(1, weight=1)
        
        # Status label
        self.status_label = ttk.Label(status_frame, text="Ready to run", 
                                     font=('Arial', 12, 'bold'), foreground='green')
        self.status_label.grid(row=0, column=0, pady=5)
        
        # Progress bar
        self.progress = ttk.Progressbar(status_frame, mode='indeterminate')
        self.progress.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=10)
        
        # Output console
        console_frame = ttk.LabelFrame(parent, text="Console Output", padding="10")
        console_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        console_frame.columnconfigure(0, weight=1)
        console_frame.rowconfigure(0, weight=1)
        
        self.console = scrolledtext.ScrolledText(console_frame, height=20, width=80, 
                                                 font=('Consolas', 9))
        self.console.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Control buttons
        btn_frame = ttk.Frame(parent)
        btn_frame.grid(row=2, column=0, pady=10)
        
        self.run_btn = ttk.Button(btn_frame, text="▶ Run Simulation", 
                                 command=self.run_simulation, style='Accent.TButton')
        self.run_btn.pack(side=tk.LEFT, padx=5)
        
        self.stop_btn = ttk.Button(btn_frame, text="⬛ Stop", 
                                  command=self.stop_simulation, state=tk.DISABLED)
        self.stop_btn.pack(side=tk.LEFT, padx=5)
        
        self.build_btn = ttk.Button(btn_frame, text="🔨 Rebuild Only", 
                                   command=self.rebuild_only)
        self.build_btn.pack(side=tk.LEFT, padx=5)
        
        # Live plot checkbox
        live_check = ttk.Checkbutton(btn_frame, text="📈 Enable Live Plotting", 
                                    variable=self.live_plot_enabled)
        live_check.pack(side=tk.LEFT, padx=15)
        
        # Configure grid weights
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(1, weight=1)
        
    def create_results_interface(self, parent):
        """Create results viewing interface"""
        
        info_label = ttk.Label(parent, text="After simulation completes, results will be generated automatically",
                              font=('Arial', 10))
        info_label.grid(row=0, column=0, pady=10)
        
        # Results buttons
        btn_frame = ttk.Frame(parent)
        btn_frame.grid(row=1, column=0, pady=10)
        
        ttk.Button(btn_frame, text="📊 Generate Analysis", 
                  command=self.generate_analysis).pack(pady=5, fill=tk.X)
        
        ttk.Button(btn_frame, text="🎬 Generate Animation", 
                  command=self.generate_animation).pack(pady=5, fill=tk.X)
        
        ttk.Button(btn_frame, text="📂 Open Results Folder", 
                  command=self.open_results_folder).pack(pady=5, fill=tk.X)
        
        ttk.Button(btn_frame, text="🖼️ View Analysis Image", 
                  command=self.view_analysis).pack(pady=5, fill=tk.X)
        
        ttk.Button(btn_frame, text="🎥 View Animation", 
                  command=self.view_animation).pack(pady=5, fill=tk.X)
        
        # Results info
        self.results_info = scrolledtext.ScrolledText(parent, height=15, width=80, 
                                                      font=('Consolas', 9))
        self.results_info.grid(row=2, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=10)
        
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(2, weight=1)
        
    def calculate_reynolds(self):
        """Calculate and display Reynolds number"""
        try:
            U = float(self.velocity_var.get())
            chord = float(self.chord_var.get())
            nu = float(self.viscosity_var.get())
            
            Re = U * chord / nu
            
            self.reynolds_label.config(text=f"Reynolds Number: {Re:.2e}")
            
        except ValueError:
            messagebox.showerror("Error", "Invalid numerical values for Reynolds calculation")
    
    def apply_parameters(self):
        """Update main.cpp with current parameter values"""
        try:
            # Read current parameters
            params = {
                'Lx': float(self.lx_var.get()),
                'Ly': float(self.ly_var.get()),
                'nx': int(self.nx_var.get()),
                'ny': int(self.ny_var.get()),
                'rho': float(self.density_var.get()),
                'nu': float(self.viscosity_var.get()),
                'U_infty': float(self.velocity_var.get()),
                'm': float(self.camber_var.get()),
                'p': float(self.camber_pos_var.get()),
                't': float(self.thickness_var.get()),
                'chord': float(self.chord_var.get()),
                'alpha': float(self.alpha_var.get()),
                'nt': int(self.nt_var.get()),
                'snapshot': int(self.snapshot_var.get())
            }
            
            # Update main.cpp
            main_cpp = self.src_dir / "main.cpp"
            
            with open(main_cpp, 'r') as f:
                content = f.read()
            
            # Replace values (this is simplified - you might want more robust parsing)
            # For now, we'll create a new version
            self.log_console("Updating main.cpp with new parameters...")
            
            # Simple replacement approach
            replacements = [
                (r'double Lx = [\d.]+', f'double Lx = {params["Lx"]}'),
                (r'Ly = [\d.]+', f'Ly = {params["Ly"]}'),
                (r'int nx = \d+', f'int nx = {params["nx"]}'),
                (r'ny = \d+', f'ny = {params["ny"]}'),
                (r'double rho = [\d.]+', f'double rho = {params["rho"]}'),
                (r'double nu\s*=\s*[\deE.+-]+', f'double nu  = {params["nu"]}'),
                (r'double U_infty = [\d.]+', f'double U_infty = {params["U_infty"]}'),
                (r'double m = [\d.]+', f'double m = {params["m"]}'),
                (r'double p = [\d.]+', f'double p = {params["p"]}'),
                (r'double t = [\d.]+', f'double t = {params["t"]}'),
                (r'double chord = [\d.]+', f'double chord = {params["chord"]}'),
                (r'double alpha = [-\d.]+', f'double alpha = {params["alpha"]}'),
                (r'int nt = \d+', f'int nt = {params["nt"]}'),
                (r'int snapshotInterval = \d+', f'int snapshotInterval = {params["snapshot"]}'),
            ]
            
            import re
            for pattern, replacement in replacements:
                content = re.sub(pattern, replacement, content)
            
            with open(main_cpp, 'w') as f:
                f.write(content)
            
            self.log_console("✓ Parameters applied successfully!")
            self.calculate_reynolds()
            
            messagebox.showinfo("Success", "Parameters updated in main.cpp!\n\nClick 'Run Simulation' to build and execute.")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to apply parameters:\n{str(e)}")
            self.log_console(f"ERROR: {str(e)}")
    
    def rebuild_only(self):
        """Rebuild the solver without running"""
        self.log_console("\n" + "="*60)
        self.log_console("Rebuilding CFD Solver...")
        self.log_console("="*60 + "\n")
        
        def build_thread():
            try:
                cmd = ["cmake", "--build", "build", "--config", "Release", "--target", "CFDSolver"]
                
                process = subprocess.Popen(
                    cmd,
                    cwd=str(self.base_dir),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )
                
                for line in process.stdout:
                    self.log_console(line.rstrip())
                
                process.wait()
                
                if process.returncode == 0:
                    self.log_console("\n✓ Build successful!")
                    messagebox.showinfo("Build Complete", "Solver rebuilt successfully!")
                else:
                    self.log_console(f"\n✗ Build failed with code {process.returncode}")
                    messagebox.showerror("Build Failed", "Check console for errors")
                    
            except Exception as e:
                self.log_console(f"\nERROR: {str(e)}")
                messagebox.showerror("Build Error", str(e))
        
        threading.Thread(target=build_thread, daemon=True).start()
    
    def run_simulation(self):
        """Build and run the simulation"""
        if self.running:
            messagebox.showwarning("Already Running", "Simulation is already in progress")
            return
        
        self.running = True
        self.run_btn.config(state=tk.DISABLED)
        self.stop_btn.config(state=tk.NORMAL)
        self.status_label.config(text="Building...", foreground='orange')
        self.progress.start()
        
        self.log_console("\n" + "="*60)
        self.log_console("STARTING CFD SIMULATION")
        self.log_console("="*60 + "\n")
        
        def run_thread():
            try:
                # Step 1: Build
                self.log_console("Step 1: Building solver...\n")
                build_cmd = ["cmake", "--build", "build", "--config", "Release", "--target", "CFDSolver"]
                
                build_process = subprocess.Popen(
                    build_cmd,
                    cwd=str(self.base_dir),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )
                
                for line in build_process.stdout:
                    self.log_console(line.rstrip())
                
                build_process.wait()
                
                if build_process.returncode != 0:
                    self.log_console("\n✗ Build failed!")
                    self.status_label.config(text="Build Failed", foreground='red')
                    self.running = False
                    self.run_btn.config(state=tk.NORMAL)
                    self.stop_btn.config(state=tk.DISABLED)
                    self.progress.stop()
                    return
                
                self.log_console("\n✓ Build successful!\n")
                
                # Step 2: Run simulation
                self.log_console("Step 2: Running simulation...\n")
                self.status_label.config(text="Running Simulation...", foreground='blue')
                
                exe_path = self.build_dir / "CFDSolver.exe"
                
                self.process = subprocess.Popen(
                    [str(exe_path)],
                    cwd=str(self.build_dir),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )
                
                # Start live plotting if enabled
                if self.live_plot_enabled.get():
                    self.root.after(0, self.start_live_plot)
                
                for line in self.process.stdout:
                    if not self.running:  # Check if stopped
                        break
                    self.log_console(line.rstrip())
                
                self.process.wait()
                
                if self.process.returncode == 0:
                    self.log_console("\n✓ Simulation complete!")
                    self.status_label.config(text="Simulation Complete!", foreground='green')
                    
                    # Stop live plotting
                    if self.plot_window:
                        self.root.after(0, self.close_live_plot)
                    
                    # Automatically generate results
                    self.log_console("\nGenerating results...")
                    self.auto_generate_results()
                    
                else:
                    self.log_console(f"\n✗ Simulation stopped (code {self.process.returncode})")
                    self.status_label.config(text="Simulation Stopped", foreground='orange')
                    
                    # Stop live plotting
                    if self.plot_window:
                        self.root.after(0, self.close_live_plot)
                
            except Exception as e:
                self.log_console(f"\nERROR: {str(e)}")
                self.status_label.config(text="Error", foreground='red')
                messagebox.showerror("Simulation Error", str(e))
            
            finally:
                self.running = False
                self.run_btn.config(state=tk.NORMAL)
                self.stop_btn.config(state=tk.DISABLED)
                self.progress.stop()
        
        threading.Thread(target=run_thread, daemon=True).start()
    
    def stop_simulation(self):
        """Stop the running simulation"""
        if self.process:
            self.log_console("\n⚠ Stopping simulation...")
            self.process.terminate()
            self.running = False
            self.status_label.config(text="Stopped by User", foreground='orange')
    
    def auto_generate_results(self):
        """Automatically generate animation and analysis"""
        self.generate_animation()
        self.generate_analysis()
    
    def generate_animation(self):
        """Generate animation from results"""
        self.log_console("\n📹 Generating animation...")
        
        animate_script = self.analysis_dir / "animate_current.py"
        
        if not animate_script.exists():
            self.log_console("✗ Animation script not found!")
            return
        
        try:
            result = subprocess.run(
                [sys.executable, str(animate_script)],
                cwd=str(self.build_dir),
                capture_output=True,
                text=True
            )
            
            self.log_console(result.stdout)
            if result.returncode == 0:
                self.log_console("✓ Animation generated!")
                self.results_info.insert(tk.END, "✓ Animation: flow_animation_100kts.gif\n")
            else:
                self.log_console(f"✗ Animation failed: {result.stderr}")
                
        except Exception as e:
            self.log_console(f"✗ Animation error: {str(e)}")
    
    def generate_analysis(self):
        """Generate analysis plots from results"""
        self.log_console("\n📊 Generating analysis...")
        
        analyze_script = self.analysis_dir / "analyze_cartesian_results.py"
        
        if not analyze_script.exists():
            self.log_console("✗ Analysis script not found!")
            return
        
        try:
            result = subprocess.run(
                [sys.executable, str(analyze_script)],
                cwd=str(self.build_dir),
                capture_output=True,
                text=True
            )
            
            self.log_console(result.stdout)
            if result.returncode == 0:
                self.log_console("✓ Analysis generated!")
                self.results_info.insert(tk.END, "✓ Analysis: complete_analysis_100kts.png\n")
            else:
                self.log_console(f"✗ Analysis failed: {result.stderr}")
                
        except Exception as e:
            self.log_console(f"✗ Analysis error: {str(e)}")
    
    def open_results_folder(self):
        """Open the results folder in file explorer"""
        if sys.platform == 'win32':
            os.startfile(str(self.build_dir))
        elif sys.platform == 'darwin':
            subprocess.run(['open', str(self.build_dir)])
        else:
            subprocess.run(['xdg-open', str(self.build_dir)])
    
    def view_analysis(self):
        """Open analysis image"""
        img_path = self.build_dir / "complete_analysis_100kts.png"
        if img_path.exists():
            if sys.platform == 'win32':
                os.startfile(str(img_path))
            else:
                webbrowser.open(str(img_path))
        else:
            messagebox.showwarning("Not Found", "Analysis image not generated yet")
    
    def view_animation(self):
        """Open animation GIF"""
        gif_path = self.build_dir / "flow_animation_100kts.gif"
        if gif_path.exists():
            if sys.platform == 'win32':
                os.startfile(str(gif_path))
            else:
                webbrowser.open(str(gif_path))
        else:
            messagebox.showwarning("Not Found", "Animation not generated yet")
    
    def log_console(self, message):
        """Thread-safe console logging"""
        def append():
            self.console.insert(tk.END, message + '\n')
            self.console.see(tk.END)
        
        self.root.after(0, append)
    
    def start_live_plot(self):
        """Open live plotting window"""
        self.plot_window = tk.Toplevel(self.root)
        self.plot_window.title("Live Flow Visualization")
        self.plot_window.geometry("1200x700")
        
        # Create matplotlib figure
        self.fig = Figure(figsize=(12, 6), dpi=100)
        
        # Create subplots
        self.ax1 = self.fig.add_subplot(121)
        self.ax2 = self.fig.add_subplot(122)
        
        # Create canvas
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_window)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Status label
        self.plot_status = ttk.Label(self.plot_window, text="Waiting for data...", 
                                     font=('Arial', 10, 'bold'))
        self.plot_status.pack(pady=5)
        
        # Start update loop
        self.latest_snapshot = -1
        self.update_live_plot()
    
    def update_live_plot(self):
        """Update live plot with latest snapshot"""
        if not self.running or not self.plot_window:
            return
        
        try:
            # Find latest snapshot
            snapshot_files = sorted(self.build_dir.glob("u_*.csv"))
            
            if snapshot_files:
                latest_file = snapshot_files[-1]
                timestep = int(latest_file.stem.split('_')[1])
                
                # Only update if new snapshot available
                if timestep > self.latest_snapshot:
                    self.latest_snapshot = timestep
                    
                    # Load data
                    u = np.genfromtxt(self.build_dir / f"u_{timestep}.csv", delimiter=',')
                    v = np.genfromtxt(self.build_dir / f"v_{timestep}.csv", delimiter=',')
                    p = np.genfromtxt(self.build_dir / f"p_{timestep}.csv", delimiter=',')
                    
                    # Get grid
                    Lx = float(self.lx_var.get())
                    Ly = float(self.ly_var.get())
                    nx = int(self.nx_var.get())
                    ny = int(self.ny_var.get())
                    
                    x = np.linspace(0, Lx, nx)
                    y = np.linspace(0, Ly, ny)
                    X, Y = np.meshgrid(x, y)
                    
                    # Compute velocity magnitude
                    vel_mag = np.sqrt(u**2 + v**2)
                    
                    # Clear and redraw
                    self.ax1.clear()
                    self.ax2.clear()
                    
                    # Plot 1: Velocity magnitude
                    levels_v = np.linspace(0, np.nanmax(vel_mag), 25)
                    cf1 = self.ax1.contourf(X, Y, vel_mag, levels=levels_v, cmap='jet')
                    self.ax1.set_xlabel('x (m)')
                    self.ax1.set_ylabel('y (m)')
                    self.ax1.set_title(f'Velocity Magnitude (m/s) - Step {timestep}')
                    self.ax1.set_aspect('equal')
                    self.ax1.set_xlim(0.3, 1.5)
                    self.ax1.set_ylim(0.6, 1.4)
                    
                    # Plot 2: Pressure
                    levels_p = np.linspace(np.nanmin(p), np.nanmax(p), 25)
                    cf2 = self.ax2.contourf(X, Y, p, levels=levels_p, cmap='RdBu_r')
                    self.ax2.set_xlabel('x (m)')
                    self.ax2.set_ylabel('y (m)')
                    self.ax2.set_title(f'Pressure (Pa) - Step {timestep}')
                    self.ax2.set_aspect('equal')
                    self.ax2.set_xlim(0.3, 1.5)
                    self.ax2.set_ylim(0.6, 1.4)
                    
                    self.fig.tight_layout()
                    self.canvas.draw()
                    
                    # Update status
                    U = float(self.velocity_var.get())
                    dt = 1.25e-5  # From main.cpp
                    physical_time = timestep * dt
                    self.plot_status.config(
                        text=f"Timestep: {timestep} | Time: {physical_time:.4f}s | "
                             f"Vel: {np.nanmin(vel_mag):.1f}-{np.nanmax(vel_mag):.1f} m/s | "
                             f"P: {np.nanmin(p):.1f}-{np.nanmax(p):.1f} Pa"
                    )
        
        except Exception as e:
            # Silently pass errors during plotting (file might be writing)
            pass
        
        # Schedule next update (every 2 seconds)
        if self.running and self.plot_window:
            self.plot_window.after(2000, self.update_live_plot)
    
    def close_live_plot(self):
        """Close live plotting window"""
        if self.plot_window:
            self.plot_window.destroy()
            self.plot_window = None


def main():
    root = tk.Tk()
    
    # Try to use modern theme
    try:
        style = ttk.Style()
        style.theme_use('vista')  # Windows modern theme
    except:
        pass
    
    app = CFDSolverGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
