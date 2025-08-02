import customtkinter as ctk
import subprocess, threading, os, zipfile, sys, subprocess
from cfd_visuals import *
import numpy as np
from PIL import Image
from customtkinter import CTkImage
from tkinter import filedialog
import pandas as pd

#sets the appearance of the oveall window
ctk.set_default_color_theme("theme.json")
ctk.set_appearance_mode("system")

#class that contains the window and all widgets
class CFDApp(ctk.CTk):
    #creates and configures the root
    def __init__(self):
        super().__init__()
        self.title("2D CFD Simulator")
        self.geometry("1000x700")
        self.resizable(False, False)
        self.build_gui()

    #function that builds all the gui components of the window
    def build_gui(self):
        #configures the grid layout of the window
        for i in range(50): self.grid_rowconfigure(i, weight=1)
        for j in range(4): self.grid_columnconfigure(j, weight=1)

        #defines the inputs and the kind of widget needed
        inputs = [("Domain Length (m):", "entry"), ("Grid Size:", "entry"), ("Inlet Velocity (m/s):", "entry"), ("Kinematic Viscosity (m²/s):", "entry"),
            ("Time Steps:", "entry"), ("Total Time (s):", "entry"), ("Time Solver:", "dropdown", ["Euler", "Heun", "RK4", "RK6", "RK8"]),
            ("Poisson Solver:", "dropdown", ["Jacobi", "Gauss-Seidel", "SOR"]), ("Cylinder Position \"x,y\":", "entry"), ("Cylinder Radius (m):", "entry")]
        
        #creates and places the main label
        self.main_label = ctk.CTkLabel(self, text="2D CFD Simulation With Optional Cylinder - Please Enter Conditions", font=("Times New Roman", 32))
        self.main_label.grid(row = 0, rowspan=10, column=0, columnspan=4)

        #creates a dictionary of inputs to access
        self.widgets = {}

        #places the labels and widgets on the grid
        for i, (label_text, widget_type, *options) in enumerate(inputs):
            #defines the correct row and column
            col_offset = 0 if i < 5 else 2
            row_base = 10 + (i % 5) * 7

            #creates and places the column
            label = ctk.CTkLabel(self, text=label_text)
            label.grid(row=row_base, column=col_offset, padx=10, pady=3, sticky="e")

            #configures and places the widget according to the type
            if widget_type == "entry":
                widget = ctk.CTkEntry(self, justify="left", placeholder_text="Enter...")
            elif widget_type == "dropdown":
                values = options[0]
                widget = ctk.CTkOptionMenu(self, values=values, anchor="w")
                widget.set(values[0])
            widget.grid(row=row_base, column=col_offset+1, padx=10, pady=5, sticky="w")
            self.widgets[label_text] = widget

        #defines the preset viscosities
        presets = [("Water", 1e-6), ("Oil", 1e-4), ("Honey", 1e-2), ("Air", 1.5e-5), ("Glycerin", 1.2e-3), ("Mercury", 1.5e-7),
            ("Blood", 3.5e-6), ("Molasses", 5), ("Alcohol", 1.2e-6), ("Gasoline", 5e-7)]

        #places the preset buttons on the grid
        for i, (name, value) in enumerate(presets):
            button = ctk.CTkButton(self, text=name, command=lambda v=value: self.set_viscosity(v), font=("Times New Roman", 27))
            button.grid(row=45 + i//4, column=i%4, padx=2, pady=6)

        #creates the button to begin the calculations
        self.run_button = ctk.CTkButton(self, text="Begin", command=self.run_simulation, font=("Times New Roman", 27))
        self.run_button.grid(row=47, column=3, padx = 2, pady=6)

        #creates the button to toggle the theme of the window
        self.theme_toggle = ctk.CTkButton(self, text="Theme", command=self.toggle_theme, font=("Times New Roman", 27))
        self.theme_toggle.grid(row=47, column=2, padx=2, pady=6)

    #function that begins the calculations
    def run_simulation(self):
        try:
            #collects the inputs from the user and makes sure that they are all valid numbers
            length = float(self.widgets["Domain Length (m):"].get())
            grid_size = int(self.widgets["Grid Size:"].get())
            U = float(self.widgets["Inlet Velocity (m/s):"].get())
            v = float(self.widgets["Kinematic Viscosity (m²/s):"].get())
            pos_str = self.widgets["Cylinder Position \"x,y\":"].get()
            x_str, y_str = pos_str.split(",")
            x = float(x_str.strip())
            y = float(y_str.strip())
            radius = float(self.widgets["Cylinder Radius (m):"].get())
            total_time = float(self.widgets["Total Time (s):"].get())
            time_steps = int(self.widgets["Time Steps:"].get())

            #checks for any logical errors in the code
            if(length <= 0 or U <= 0 or v <= 0 or total_time <= 0): 
                self.main_label.configure(text="Length, Velocity, Time, and Viscosity Must be Positive")
                return

            #checks for grid or time steps below two
            if grid_size < 2 or time_steps < 2: 
                self.main_label.configure(text="Grid Size and Time Divisions Must be at Least 2")
                return

            #checks for negative radius for the cylinder
            if radius < 0: 
                self.main_label.configure(text="Cylinder Radius Must be Positive or 0 for Off")
                return
            
            #checks for cylinder center in grid
            if not (0 < x < length) or not (0 < y < length): 
                self.main_label.configure(text="Cylinder Center Must be Within Grid")
                return
        
            #checks for cylinder fit in grid
            if (x + radius >= length or x - radius <= 0 or y + radius >= length or y - radius <= 0): 
                self.main_label.configure(text="Cylinder Must be Completly Within Grid")
                return        
            
            #saves all the files to the parameter file
            with open("parameters.txt", "w") as file:
                file.write(f"{length}\n")
                file.write(f"{grid_size}\n")
                file.write(f"{U}\n")
                file.write(f"{v}\n")
                file.write(f"{total_time}\n")
                file.write(f"{time_steps}\n")
                file.write(f"{self.widgets['Poisson Solver:'].get()}\n")
                file.write(f"{self.widgets['Time Solver:'].get()}\n")
                file.write(f"{x}\n")
                file.write(f"{y}\n")
                file.write(f"{radius}\n")

            #saves the results as class variables
            self.inlet = U
            self.visc = v
            self.domain_length = length
            self.step_length = grid_size
            self.time = total_time
            self.iterations = time_steps
            self.cyl_x = x
            self.cyl_y = y
            self.cyl_r = radius

            #calculates values based on inputs
            self.dx = length/grid_size
            self.dt = total_time/time_steps
            
        #displays an error message if type is invalid
        except Exception as e:
            self.main_label.configure(text="Please Make Sure Inputs are Valid Numbers")
            print(e)
            return
        
        #clears the results directory
        results_dir = "results"
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
        for file in os.listdir(results_dir):
            file_path = os.path.join(results_dir, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
        
        #clears the graphs directory
        results_dir = "graphs"
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
        for file in os.listdir(results_dir):
            file_path = os.path.join(results_dir, file)
            if os.path.isfile(file_path):
                os.remove(file_path)
                
        #removes all widgets from the screen
        for widget in self.winfo_children():
            if widget != self.main_label:
                widget.destroy()
        
        #restes the progress file
        with open("progress.txt", "w") as f:
            f.write("0")
            
        #creates a progress bar to track the backend progress
        self.progress_bar = ctk.CTkProgressBar(self, orientation="horizontal", width = 500, height=20, corner_radius=10, fg_color=("#d9e6f5", "#1f3a5f"), progress_color=("#173a5e", "#d8ecff"))
        self.progress_bar.set(0)
        self.progress_bar.grid(row=23, column=0, columnspan=4)
        self.simulation_running = True
        self.update_progress()

        #configures the main label
        self.main_label.configure(text=f"Running Simulation\nProgress: 0/{self.iterations}", font=("Times New Roman", 54))
        self.main_label.grid(row=21, rowspan=1)

        #calls the function to run the backend solver in a parralel thread
        thread = threading.Thread(target=self.call_fortran, daemon=True)
        thread.start()
    
    #function that sets the viscosity according to the presets
    def set_viscosity(self, visc):
        visc_entry = self.widgets["Kinematic Viscosity (m²/s):"]
        visc_entry.delete(0, "end")
        visc_entry.insert(0, str(visc))

    #function that toggles the theme of the window
    def toggle_theme(self):
        current = ctk.get_appearance_mode().lower()
        ctk.set_appearance_mode("light") if current == "dark" else ctk.set_appearance_mode("dark")

    #function that calls the fortran backend solver
    def call_fortran(self):
        try:
            result = subprocess.run(["cfd.exe"], capture_output=True, text=True, check=True)
        #displays an error message if simulation failed
        except subprocess.CalledProcessError as e:
            self.main_label.configure(text="Simulation failed. Please try again")
            print("Simulation failed:", e)
            print("stdout:\n", e.stdout)
            print("stderr:\n", e.stderr)
            return
        except Exception as e:
            self.main_label.configure(text="Unexpected error. Check console.")
            print("Error:", e)
            return
        self.after(0, self.after_simulation)

    #function that updates the progress bar
    def update_progress(self):
        try:
            with open("progress.txt", "r") as f:
                step = int(f.readline().strip())
                total = int(self.iterations)
                self.progress_bar.set(step / total)
                self.main_label.configure(text=f"Running Simulation\nProgress: {step}/{total}")
        except:
            pass

        if self.simulation_running:
            self.after(100, self.update_progress)

    #function that processes the results
    def after_simulation(self):
        #configures gui widgets
        self.main_label.configure(text="Simulation Completed:\nProcessing Results")
        self.progress_bar.destroy()
        self.update()

        #reads the scalar values
        self.time = np.loadtxt("results/time.csv")
        self.cd = np.loadtxt("results/cd.csv")
        self.cl = np.loadtxt("results/cl.csv")
        self.lift = np.loadtxt("results/lift.csv")
        self.drag = np.loadtxt("results/drag.csv")

        #sets the first unstable value to the stable second value
        self.cd[0] = self.cd[1]
        self.cl[0] = self.cl[1]
        self.drag[0] = self.drag[1]
        self.lift[0] = self.lift[1]

        #function to load the value fields
        def load_field(path):
            frames = []
            block = []
            with open(path) as f:
                for line in f:
                    if line.startswith("#"):
                        if block:
                            frames.append(np.array(block))
                            block = []
                    else:
                        row = [float(val) for val in line.strip().split()]
                        block.append(row)
            data = np.array(frames)
            data = np.flip(np.transpose(data, axes=(0, 2, 1)), axis=2)
            data = np.flip(data, axis=2)
            return data
        
        #calls the function to read all fields
        self.psi = load_field("results/psi.csv")
        self.omega = load_field("results/omega.csv")
        self.u = load_field("results/x_vel.csv")
        self.v = load_field("results/y_vel.csv")
        self.pressure = load_field("results/pressure.csv")

        #gets derived quantities
        dx = self.dx
        dy = self.dx
        U = self.inlet

        #defines density as constant (incompressible)
        rho = 1.0

        #calculates the velocity and vorticity magnitudes
        self.vel_mag = np.sqrt(self.u**2 + self.v**2)
        self.vort_mag = np.abs(self.omega)

        #calculates the kinetic energy
        self.kinetic_energy = 0.5 * np.sum(self.u**2 + self.v**2, axis=(1, 2)) * dx * dy

        #calculates enstrophy
        self.enstrophy = 0.5 * np.sum(self.omega**2, axis=(1, 2)) * dx * dy

        #calculates total circulation
        self.total_vorticity = np.sum(self.omega, axis=(1, 2)) * dx * dy

        #calculates average pressure drop
        inlet_p = self.pressure[:, :, 0]
        outlet_p = self.pressure[:, :, -1]
        self.delta_p = np.mean(inlet_p, axis=1) - np.mean(outlet_p, axis=1)
    
        #calculates mass flow rate
        self.mass_flow_in = np.sum(self.u[:, :, 0], axis=1) * dy
        self.mass_flow_out = np.sum(self.u[:, :, -1], axis=1) * dy

        #sets first value of scalars to stable time
        self.kinetic_energy = np.insert(self.kinetic_energy, 0, self.kinetic_energy[0])
        self.enstrophy = np.insert(self.enstrophy, 0, self.enstrophy[0])
        self.total_vorticity = np.insert(self.total_vorticity, 0, self.total_vorticity[0])
        self.delta_p = np.insert(self.delta_p, 0, self.delta_p[0])
        self.mass_flow_in = np.insert(self.mass_flow_in, 0, self.mass_flow_in[0])
        self.mass_flow_out = np.insert(self.mass_flow_out, 0, self.mass_flow_out[0])

        #smooths all the graphs for instability
        def smooth_initial(data, count=5):
            if len(data) > count: data[:count] = data[count]
            return data
        self.kinetic_energy = smooth_initial(self.kinetic_energy)
        self.enstrophy = smooth_initial(self.enstrophy)
        self.total_vorticity = smooth_initial(self.total_vorticity)
        self.delta_p = smooth_initial(self.delta_p)
        self.mass_flow_in = smooth_initial(self.mass_flow_in)
        self.mass_flow_out = smooth_initial(self.mass_flow_out)

        #configures main label for results
        self.main_label.configure(text="Simulation Results", font=("Times New Roman", 42))
        self.main_label.grid(row=0, column=0, columnspan=4, pady=10)

        #defines a dictionary of all the plots
        self.plot_paths = {
            "Drag Coefficient": "drag_coefficient.png",
            "Lift Coefficient": "lift_coefficient.png",
            "Drag Force": "drag.png",
            "Lift Force": "lift.png",
            "Kinetic Energy": "kinetic_energy.png",
            "Enstrophy": "enstrophy.png",
            "Total Circulation": "circulation.png",
            "Pressure Drop": "pressure_drop.png",
            "Mass Flow In": "mass_in_flow.png",
            "Mass Flow Out": "mass_out_flow.png"
        }

        #creates buttons for each calar plot
        for i, (label, filename) in enumerate(self.plot_paths.items()):
            btn = ctk.CTkButton(self, text=label, font=("Times New Roman", 21), command=lambda f=filename, l=label: self.open_plot_window(f, l))
            btn.grid(row= 11 + (i // 2) * 3, column=2*(i % 2), sticky="nsew", padx=10)

        #adds buttons for the field viewers
        self.interactive_fields = {"Velocity Magnitude": self.vel_mag, "Pressure": self.pressure, "Vorticity": self.omega, "Streamfunction": self.psi}
        plot_types = ["Heatmap", "Streamlines", "Quiver"]
        col_map = {"Heatmap": 0, "Streamlines": 1, "Quiver": 2}
        for i, field_name in enumerate(self.interactive_fields):
            for plot_type in plot_types:
                btn = ctk.CTkButton(self, text=f"{field_name} - {plot_type}", font=("Times New Roman", 21), command=lambda f=field_name, p=plot_type: self.open_field_viewer(f, p))
                row = 26 + i * 3
                col = col_map[plot_type]
                btn.grid(row=row, column=col, sticky="nsew", padx=10)

        #creates the buttons to export graphs
        export_graphs_btn = ctk.CTkButton(self, text="Export Opened Graphs", font=("Times New Roman", 21), command=self.export_all_graphs)
        export_graphs_btn.grid(row=38, column=0, sticky="nsew", padx=10)

        #creates the button to export results
        export_data_btn = ctk.CTkButton(self, text="Export All Data", font=("Times New Roman", 21), command=self.export_all_data)
        export_data_btn.grid(row=38, column=2, sticky="nsew", padx=10)

        #creates the button to restart the program
        restart_btn = ctk.CTkButton(self, text="Restart", font=("Times New Roman", 21), command=self.restart_app)
        restart_btn.grid(row=41, column=0, sticky="nsew", padx=10)

        #creates the button to exit the program
        exit_btn = ctk.CTkButton(self, text="Exit", font=("Times New Roman", 21), command=self.exit_app)
        exit_btn.grid(row=41, column=2, sticky="nsew", padx=10)
        
    #function that plots the static plots
    def open_plot_window(self, filename, title):
        #creates the graph that the user specifcies if not already craeted
        img_path = os.path.join("graphs", filename)
        theme = ctk.get_appearance_mode().lower()
        if not os.path.exists(img_path):
            #
            data_map = {
                "Drag Coefficient": self.cd,
                "Lift Coefficient": self.cl,
                "Drag Force": self.drag,
                "Lift Force": self.lift,
                "Kinetic Energy": self.kinetic_energy,
                "Enstrophy": self.enstrophy,
                "Total Circulation": self.total_vorticity,
                "Pressure Drop": self.delta_p,
                "Mass Flow In": self.mass_flow_in,
                "Mass Flow Out": self.mass_flow_out
            }
        if title in data_map: plot_scalar_vs_time(self.time, data_map[title], filename, title, theme)

        #creates and configures a top level to the root
        win = ctk.CTkToplevel(self)
        win.title(title)
        win.geometry("720x520")
        win.resizable(False, False)

        #loads and resizes the image
        img_path = os.path.join("graphs", filename)
        pil_img = Image.open(img_path)
        img = CTkImage(light_image=pil_img, dark_image=pil_img, size=(700, 500))
        panel = ctk.CTkLabel(win, image=img, text="")
        panel.image = img
        panel.pack(pady=5)

    #function that visualzes the heatmaps
    def open_field_viewer(self, field_name, plot_type):
        #collects the data and makes the slider max
        data = self.interactive_fields[field_name]
        slider_max = data.shape[0] - 1

        #creates the top level from the root
        viewer = ctk.CTkToplevel(self)
        viewer.title(f"{field_name} Viewer")
        viewer.geometry("900x600")
        viewer.resizable(False, False)

        #configures the layout
        for i in range(11):
            viewer.grid_rowconfigure(i, weight=1 if i < 10 else 0)
        for j in range(10):
            viewer.grid_columnconfigure(j, weight=1)

        #creates the image
        img_label = ctk.CTkLabel(viewer, text="")
        img_label.grid(row=0, column=0, rowspan=10, columnspan=10, sticky="nsew")

        #function that updates the plot according to slider position
        def update_plot(t_index):
            #gets the specified data and calls the creation function
            time_val = self.time[t_index]
            theme = ctk.get_appearance_mode().lower()
            filename = f"__temp_{field_name.lower().replace(' ', '_')}_{t_index}.png"
            filepath = os.path.join("graphs", filename)
            if not os.path.exists(filepath):
                if plot_type == "Heatmap": plot_heatmap_field(data[t_index], filename, field_name, theme, time_val, cyl_x=self.cyl_x, cyl_y=self.cyl_y, cyl_r=self.cyl_r, L=self.domain_length, N=self.step_length)
                elif plot_type == "Streamlines": plot_streamline_field(self.u[t_index], self.v[t_index], data[t_index], filename, f"{field_name}", theme, time_val, cyl_x=self.cyl_x, cyl_y=self.cyl_y, cyl_r=self.cyl_r, L=self.domain_length, N=self.step_length)
                elif plot_type == "Quiver": plot_quiver_field(self.u[t_index], self.v[t_index], data[t_index], filename, f"{field_name}", theme, time_val, cyl_x=self.cyl_x, cyl_y=self.cyl_y, cyl_r=self.cyl_r, L=self.domain_length, N=self.step_length)

            #loads the image to the gui
            pil_img = Image.open(os.path.join("graphs", filename))
            img = CTkImage(light_image=pil_img, dark_image=pil_img, size=(760, 550))
            img_label.configure(image=img, text="")
            img_label.image = img

        #creates the slider and bind it to the updating function
        update_plot(0)
        slider = ctk.CTkSlider(viewer, from_=0, to=slider_max, number_of_steps=slider_max, orientation="horizontal", width=600)
        slider.grid(row=10, column=0, columnspan=10, sticky="swe", pady=(10, 5))
        def on_slider_change(val): update_plot(int(round(val)))
        slider.set(0)
        slider.configure(command=on_slider_change)
    
    #function that exports the opened graphs graphs
    def export_all_graphs(self):
        #checks which graphs have been opened
        output_dir = "graphs"
        if not os.listdir(output_dir):
            self.main_label.configure(text="No graphs to export\nPlease open a graph to save it for export")
            self.after(3000, self.main_label.configure(text="Simulation Results"))
            return

        #asks user to select the destination
        zip_path = filedialog.asksaveasfilename(defaultextension=".zip", initialfile="CFD_Graphs.zip", title="Save All Graphs", filetypes=[("ZIP Archive", "*.zip")])
        if not zip_path: return

        #exports all the current graphs to the selected directory
        with zipfile.ZipFile(zip_path, 'w') as z:
            for file in os.listdir(output_dir):
                if file.endswith(".png"): z.write(os.path.join(output_dir, file), arcname=file)

    #function that exports all the data
    def export_all_data(self):
        #asks user to select the destination
        zip_path = filedialog.asksaveasfilename(defaultextension=".zip", initialfile="CFD_Data.zip", title="Save All Data", filetypes=[("ZIP Archive", "*.zip")])
        if not zip_path: return

        #makes a temporary data directory to export all the data
        temp_dir = "temp_data_export"
        os.makedirs(temp_dir, exist_ok=True)

        #loads the field data
        full_fields = {
            "psi": self.psi, "omega": self.omega, "pressure": self.pressure,
            "u": self.u, "v": self.v, "velocity_magnitude": self.vel_mag,
            "vorticity_magnitude": self.vort_mag
        }

        #saves the field data to the temp directoy
        for name, array in full_fields.items():
            flat = array.reshape(array.shape[0], -1)
            df = pd.DataFrame(flat)
            df.insert(0, "time", self.time[:array.shape[0]])
            df.to_csv(os.path.join(temp_dir, f"{name}.csv"), index=False)

        #saves the scalar data
        scalar_fields = {
            "cd": self.cd, "cl": self.cl, "drag": self.drag, "lift": self.lift,
            "kinetic_energy": self.kinetic_energy, "enstrophy": self.enstrophy,
            "total_circulation": self.total_vorticity, "pressure_drop": self.delta_p,
            "mass_flow_in": self.mass_flow_in, "mass_flow_out": self.mass_flow_out
        }

        #saves the scalar data to the temp directoy
        for name, data in scalar_fields.items():
            df = pd.DataFrame({"time": self.time, name: data})
            df.to_csv(os.path.join(temp_dir, f"{name}.csv"), index=False)

        #saves all the files to the selected path
        with zipfile.ZipFile(zip_path, "w") as zipf:
            for file in os.listdir(temp_dir):
                zipf.write(os.path.join(temp_dir, file), arcname=file)

        #removes the temporary directory
        for file in os.listdir(temp_dir):
            os.remove(os.path.join(temp_dir, file))
        os.rmdir(temp_dir)

    #function that restarts the program
    def restart_app(self):
        self.destroy()
        subprocess.call([sys.executable, sys.argv[0]])

    #function that exits the program
    def exit_app(self):
        self.destroy()

#begins the program
if __name__ == "__main__":
    app = CFDApp()
    app.mainloop()
