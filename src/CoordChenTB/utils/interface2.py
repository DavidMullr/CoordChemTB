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
        self.sdf_file_list = []
        self.ligand_names = []
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
        self.oxstate_var = tk.StringVar()
        self.selected_ligand = tk.StringVar()
        self.geom_orb_var = tk.StringVar()
        self.lig_list_accum = []

        ttk.Label(frame, text="Ligand Data SDFs:").grid(row=0, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.sdf_file_var, width=50).grid(row=0, column=1, padx=5)
        ttk.Button(frame, text="Browse", command=self.browse_sdf_file).grid(row=0, column=2, padx=5)

        ttk.Label(frame, text="Metal:").grid(row=1, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.metal_orb_var).grid(row=1, column=1, padx=5)

        ttk.Label(frame, text="Oxidation State of Metal:").grid(row=2, column=0, sticky='w', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.oxstate_var).grid(row=2, column=1, padx=5)

        ttk.Label(frame, text="Select Ligand:").grid(row=3, column=0, sticky='w', padx=5, pady=5)
        self.ligand_dropdown = ttk.Combobox(frame, textvariable=self.selected_ligand, values=[], state="readonly")
        self.ligand_dropdown.grid(row=3, column=1, padx=5)
        ttk.Button(frame, text="Add Ligand", command=self.add_ligand).grid(row=3, column=2, padx=5)

        ttk.Label(frame, text="Ligands Selected:").grid(row=4, column=0, sticky='w', padx=5)
        self.lig_list_display = tk.StringVar()
        ttk.Entry(frame, textvariable=self.lig_list_display, state="readonly", width=60).grid(row=5, column=0, columnspan=3, padx=5)

        ttk.Label(frame, text="Geometry:").grid(row=6, column=0, sticky='w', padx=5, pady=5)
        ttk.Combobox(
            frame,
            textvariable=self.geom_orb_var,
            values=["octahedral", "tetrahedral", "square planar"],
            state="readonly"
        ).grid(row=6, column=1, padx=5)

        ttk.Button(frame, text="Compute Orbital Splitting", command=self.compute_orbital).grid(row=7, column=0, columnspan=3, pady=10)

    def browse_sdf_file(self):
        import os
        paths = filedialog.askopenfilenames(title="Select Ligand SDF Files", filetypes=[("SDF Files", "*.sdf")])
        if paths:
            self.sdf_file_list = list(paths)
            self.sdf_file_var.set("; ".join(paths))
            ligand_set = set()
            for path in paths:
                d = orbital.load_ligand_data(path)
                if d:
                    ligand_set.update(d.keys())
            self.ligand_dropdown['values'] = sorted(ligand_set)

    def add_ligand(self):
        ligand = self.selected_ligand.get()
        if ligand and ligand not in self.lig_list_accum:
            self.lig_list_accum.append(ligand)
            self.lig_list_display.set(", ".join(self.lig_list_accum))

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
            import io, sys

            metal_sym = self.metal_orb_var.get()
            ox_state = int(self.oxstate_var.get())
            lig_names = self.lig_list_accum
            geom = self.geom_orb_var.get()

            data = {}
            for path in self.sdf_file_list:
                d = orbital.load_ligand_data(path)
                if d:
                    data.update(d)

            metal = orbital.get_metal_data(metal_sym)
            deltas, charges = [], []
            for name in lig_names:
                info = data.get(name) or data.get(orbital.find_closest_ligand(name, data)) or {"Delta": 20000, "Charge": 0}
                deltas.append(info["Delta"])
                charges.append(info["Charge"])

            d_elec = orbital.get_d_electron_count(metal["atomic_number"], ox_state)
            avg_delta = sum(deltas) / len(deltas)
            spin = orbital.predict_spin_state(avg_delta, orbital.get_pairing_energy(metal["block"]))

            # Capture CFSE output
            buf = io.StringIO()
            old_stdout = sys.stdout
            sys.stdout = buf
            orbital.plot_cf(d_elec, avg_delta, spin, geom, metal["block"], label=f"{metal_sym} {spin}")
            sys.stdout = old_stdout
            output = buf.getvalue()
            lines = output.splitlines()

            pre = next((ln for ln in lines if "📉 CFSE" in ln and "Total" not in ln and "pairing" not in ln), None)
            pen = next((ln for ln in lines if "pairing penalty" in ln), None)
            tot = next((ln for ln in lines if "Total CFSE" in ln), None)

            msg = "\n".join(filter(None, [pre, pen, tot]))
            if msg:
                messagebox.showinfo("CFSE Energies", msg)
            else:
                messagebox.showinfo("CFSE Energies", "No energy details found.")
        except Exception as e:
            sys.stdout = old_stdout
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    app = CoordinationGUI()
    app.mainloop()
