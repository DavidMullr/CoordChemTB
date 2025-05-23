import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from PIL import Image, ImageTk
from complexdrawerfinal import create_complex_from_ligand_dict
import complexorbitalssplitting as orbital
from metals_db import METALS

class CoordinationGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Coordination Chemistry Toolkit")
        self.geometry("1000x600")
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill='both', expand=True)
        self.init_complex_tab()
        self.init_orbital_tab()

    def init_complex_tab(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Complex Builder")

        ttk.Label(frame, text="Metal:").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        metals = [m['symbol'] for m in METALS]
        self.metal_var = tk.StringVar()
        ttk.Combobox(frame, textvariable=self.metal_var, values=metals, state="readonly").grid(row=0, column=1, padx=5, pady=5)

        ttk.Label(frame, text="Ligands (Name:Count, e.g. py:2, bipy:1):").grid(row=1, column=0, columnspan=2, sticky='w', padx=5)
        self.ligand_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.ligand_var, width=60).grid(row=2, column=0, columnspan=2, padx=5)

        ttk.Label(frame, text="Geometry:").grid(row=3, column=0, sticky='w', padx=5, pady=5)
        self.geom_var = tk.StringVar()
        ttk.Combobox(
            frame,
            textvariable=self.geom_var,
            values=["octahedral", "tetrahedral", "square planar"],
            state="readonly"
        ).grid(row=3, column=1, padx=5, pady=5)

        ttk.Button(frame, text="Draw Complex", command=self.draw_complex).grid(row=4, column=0, columnspan=2, pady=10)

        self.canvas = tk.Canvas(frame, bg='white')
        self.canvas.grid(row=5, column=0, columnspan=2, sticky='nsew', padx=5, pady=5)
        frame.rowconfigure(5, weight=1)
        frame.columnconfigure(1, weight=1)

    def init_orbital_tab(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Orbital Splitting")

        self.sdf_file_var = tk.StringVar()
        self.metal_orb_var = tk.StringVar()
        self.charge_var = tk.StringVar()
        self.lig_list_var = tk.StringVar()
        self.geom_orb_var = tk.StringVar()

        ttk.Label(frame, text="Ligand Data SDF:").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.sdf_file_var, width=50).grid(row=0, column=1, padx=5)
        ttk.Button(frame, text="Browse", command=self.browse_sdf_file).grid(row=0, column=2, padx=5)

        ttk.Label(frame, text="Metal:").grid(row=1, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.metal_orb_var).grid(row=1, column=1, padx=5)

        ttk.Label(frame, text="Total Complex Charge:").grid(row=2, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.charge_var).grid(row=2, column=1, padx=5)

        ttk.Label(frame, text="Ligands (comma-separated):").grid(row=3, column=0, columnspan=2, sticky='w', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.lig_list_var, width=60).grid(row=4, column=0, columnspan=3, padx=5)

        ttk.Label(frame, text="Geometry:").grid(row=5, column=0, sticky='w', padx=5, pady=5)
        ttk.Combobox(
            frame,
            textvariable=self.geom_orb_var,
            values=["octahedral", "tetrahedral", "square planar"],
            state="readonly"
        ).grid(row=5, column=1, padx=5)

        ttk.Button(frame, text="Compute Orbital Splitting", command=self.compute_orbital).grid(row=6, column=0, columnspan=3, pady=10)
        ttk.Button(frame, text="Manual d-Electron Tool", command=self.open_manual_visualizer).grid(row=7, column=0, columnspan=3, pady=5)

    def browse_sdf_file(self):
        path = filedialog.askopenfilename(title="Select all_ligands.sdf", filetypes=[("SDF Files", "*.sdf")])
        if path:
            self.sdf_file_var.set(path)

    def draw_complex(self):
        try:
            metal = self.metal_var.get()
            geometry = self.geom_var.get().replace(' ', '_')
            ligand_counts = {
                name.strip(): int(cnt)
                for part in self.ligand_var.get().split(',')
                if ':' in part for name, cnt in [part.split(':')]
            }
            img = create_complex_from_ligand_dict(metal, ligand_counts, geometry)
            self.photo = ImageTk.PhotoImage(img)
            self.canvas.update_idletasks()
            w, h = self.canvas.winfo_width(), self.canvas.winfo_height()
            self.canvas.delete("all")
            self.canvas.create_image(w/2, h/2, image=self.photo, anchor='center')
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def compute_orbital(self):
        try:
            sdf = self.sdf_file_var.get()
            metal_sym = self.metal_orb_var.get()
            total_charge = int(self.charge_var.get())
            lig_names = [l.strip() for l in self.lig_list_var.get().split(',') if l.strip()]
            geom = self.geom_orb_var.get()
            data = orbital.load_ligand_data(sdf)
            metal = orbital.get_metal_data(metal_sym)
            deltas, charges = [], []
            for name in lig_names:
                info = data.get(name) or data.get(orbital.find_closest_ligand(name, data)) or {"Delta": 20000, "Charge": 0}
                deltas.append(info["Delta"])
                charges.append(info["Charge"])
            ligand_charge = sum(charges)
            ox_state = total_charge - ligand_charge
            d_elec = orbital.get_d_electron_count(metal["atomic_number"], ox_state)
            avg_delta = sum(deltas) / len(deltas)
            spin = orbital.predict_spin_state(avg_delta, orbital.get_pairing_energy(metal["block"]))
            orbital.plot_cf(d_elec, avg_delta, spin, geom, metal["block"], label=f"{metal_sym} {spin}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def open_manual_visualizer(self):
        ManualFieldVisualizer(self)

class ManualFieldVisualizer(tk.Toplevel):
    def __init__(self, master=None):
        super().__init__(master)
        self.title("Manual d-Electron Field Visualization")
        self.geometry("400x300")

        self.geom_var = tk.StringVar(value="octahedral")
        self.delec_var = tk.StringVar(value="6")
        self.delta_var = tk.StringVar(value="20000")
        self.spin_var = tk.StringVar(value="High-spin")
        self.block_var = tk.StringVar(value="3d")

        ttk.Label(self, text="Geometry:").grid(row=0, column=0, sticky="w", padx=10, pady=5)
        ttk.Combobox(self, textvariable=self.geom_var, values=["octahedral", "tetrahedral", "square planar"], state="readonly").grid(row=0, column=1, padx=10, pady=5)

        ttk.Label(self, text="d-Electron Count:").grid(row=1, column=0, sticky="w", padx=10, pady=5)
        ttk.Entry(self, textvariable=self.delec_var).grid(row=1, column=1, padx=10, pady=5)

        ttk.Label(self, text="Δ (cm⁻¹):").grid(row=2, column=0, sticky="w", padx=10, pady=5)
        ttk.Entry(self, textvariable=self.delta_var).grid(row=2, column=1, padx=10, pady=5)

        ttk.Label(self, text="Spin State:").grid(row=3, column=0, sticky="w", padx=10, pady=5)
        ttk.Combobox(self, textvariable=self.spin_var, values=["High-spin", "Low-spin"], state="readonly").grid(row=3, column=1, padx=10, pady=5)

        ttk.Label(self, text="Block:").grid(row=4, column=0, sticky="w", padx=10, pady=5)
        ttk.Combobox(self, textvariable=self.block_var, values=["3d", "4d", "5d", "f"], state="readonly").grid(row=4, column=1, padx=10, pady=5)

        ttk.Button(self, text="Plot", command=self.plot_diagram).grid(row=5, column=0, columnspan=2, pady=15)

    def plot_diagram(self):
        try:
            orbital.plot_cf(
                d_elec=int(self.delec_var.get()),
                delta=float(self.delta_var.get()),
                spin=self.spin_var.get(),
                geom=self.geom_var.get(),
                geom_block=self.block_var.get(),
                label=f"{self.spin_var.get()} {self.geom_var.get()}"
            )
        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    app = CoordinationGUI()
    app.mainloop()
