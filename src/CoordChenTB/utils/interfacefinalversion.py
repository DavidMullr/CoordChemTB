import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from PIL import Image, ImageTk
from complexdrawerfinal import create_complex_from_ligand_dict, load_ligands_from_folder, LIGAND_FOLDER
import complexorbitalssplitting as orbital
from metals_db import METALS

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


class CoordinationGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Coordination Chemistry Toolkit")
        self.geometry("1000x600")
        self.configure(bg='#2e2e2e')  # Set dark background
        title_label = ttk.Label(self, text="CoordChemTB", font=("Helvetica", 28, "bold"))
        title_label.pack(pady=15)

        # Global dark theme styling
        style = ttk.Style(self)
        style.theme_use('clam')
        style.configure('.', background='#2e2e2e', foreground='white', fieldbackground='#3a3a3a')
        style.configure('TLabel', background='#2e2e2e', foreground='white', anchor='center')
        style.configure('TEntry', fieldbackground='#3a3a3a', foreground='white', justify='center')
        style.configure('TCombobox', fieldbackground='#3a3a3a', foreground='white', justify='center')
        style.map('TCombobox',
                fieldbackground=[('readonly', '#3a3a3a')],
                foreground=[('readonly', 'white')])
        style.configure('TButton', background='#444', foreground='white')

        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill='both', expand=True)

        self.sdf_file_list = []
        self.lig_map = load_ligands_from_folder(LIGAND_FOLDER)
        self.ligand_names = sorted(self.lig_map.keys())

        self.init_complex_tab()
        self.init_orbital_tab()
        self.init_info_tab()


    def init_complex_tab(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Complex Builder")

        # Configure columns to expand for centering
        for i in range(2):
            frame.columnconfigure(i, weight=1)

        ttk.Label(frame, text="Metal:").grid(row=0, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        metals = [m['symbol'] for m in METALS]
        self.metal_var = tk.StringVar()
        ttk.Combobox(frame, textvariable=self.metal_var, values=metals, state="readonly").grid(row=1, column=0, columnspan=2, sticky='ew', padx=5)

        ttk.Label(frame, text="Ligands (Name:Count, e.g. py:2, bipy:1):").grid(row=2, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        self.ligand_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.ligand_var, width=60).grid(row=3, column=0, columnspan=2, sticky='ew', padx=5)

        ttk.Label(frame, text="Geometry:").grid(row=4, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        self.geom_var = tk.StringVar()
        ttk.Combobox(
            frame,
            textvariable=self.geom_var,
            values=["octahedral", "tetrahedral", "square planar"],
            state="readonly"
        ).grid(row=5, column=0, columnspan=2, sticky='ew', padx=5)

        ttk.Label(frame, text="Bond Length (for bulky ligands):").grid(row=6, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        self.bond_length_var = tk.StringVar(value="1.0")
        ttk.Entry(frame, textvariable=self.bond_length_var).grid(row=7, column=0, columnspan=2, sticky='ew', padx=5)

        ttk.Button(frame, text="Draw Complex", command=self.draw_complex).grid(row=8, column=0, columnspan=2, sticky='ew', pady=10)

        self.canvas = tk.Canvas(frame, bg='white')
        self.canvas.grid(row=9, column=0, columnspan=2, sticky='nsew', padx=5, pady=5)
        frame.rowconfigure(9, weight=1)


    def init_orbital_tab(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Orbital Splitting")

        # Configure columns to expand evenly for centering
        for i in range(3):
            frame.columnconfigure(i, weight=1)

        self.metal_orb_var = tk.StringVar()
        self.oxstate_var = tk.StringVar()
        self.selected_ligand = tk.StringVar()
        self.geom_orb_var = tk.StringVar()
        self.lig_list_accum = {} 


        ttk.Label(frame, text="Metal:").grid(row=1, column=0, columnspan=3, sticky='ew', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.metal_orb_var).grid(row=2, column=0, columnspan=3, sticky='ew', padx=5)

        ttk.Label(frame, text="Oxidation State of Metal:").grid(row=3, column=0, columnspan=3, sticky='ew', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.oxstate_var).grid(row=4, column=0, columnspan=3, sticky='ew', padx=5)

        ttk.Label(frame, text="Select Ligand:").grid(row=5, column=0, columnspan=3, sticky='ew', padx=5, pady=5)
        self.ligand_dropdown = ttk.Combobox(frame, textvariable=self.selected_ligand, values=self.ligand_names, state="readonly")
        self.ligand_count_var = tk.StringVar(value="1")
        ttk.Entry(frame, textvariable=self.ligand_count_var, width=5).grid(row=6, column=2, padx=5)
        self.ligand_dropdown.grid(row=6, column=0, columnspan=2, sticky='ew', padx=5)
        ttk.Button(frame, text="Add Ligand", command=self.add_ligand).grid(row=6, column=2, sticky='ew', padx=5)

        ttk.Label(frame, text="Ligands Selected:").grid(row=7, column=0, columnspan=3, sticky='ew', padx=5)
        self.lig_list_display = tk.StringVar()
        ttk.Entry(frame, textvariable=self.lig_list_display, state="readonly", width=60).grid(row=8, column=0, columnspan=3, sticky='ew', padx=5)

        ttk.Label(frame, text="Geometry:").grid(row=9, column=0, columnspan=3, sticky='ew', padx=5, pady=5)
        ttk.Combobox(
            frame,
            textvariable=self.geom_orb_var,
            values=["octahedral", "tetrahedral", "square planar"],
            state="readonly"
        ).grid(row=10, column=0, columnspan=3, sticky='ew', padx=5)

        ttk.Button(frame, text="Compute Orbital Splitting", command=self.compute_orbital).grid(row=11, column=0, columnspan=3, sticky='ew', pady=10)

    def init_info_tab(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Metal & Ligand Info")

        for i in range(2):
            frame.columnconfigure(i, weight=1)

        self.lookup_var = tk.StringVar()
        ttk.Label(frame, text="Enter Metal or Ligand Name:").grid(row=0, column=0, columnspan=2, sticky='ew', pady=10)
        ttk.Entry(frame, textvariable=self.lookup_var).grid(row=1, column=0, columnspan=2, sticky='ew', padx=10)

        ttk.Button(frame, text="Lookup Info", command=self.lookup_info).grid(row=2, column=0, columnspan=2, sticky='ew', padx=20, pady=10)

        self.info_display = tk.Text(frame, height=20, bg='#3a3a3a', fg='white', wrap='word')
        self.info_display.grid(row=3, column=0, columnspan=2, sticky='nsew', padx=10, pady=10)
        frame.rowconfigure(3, weight=1)

    def lookup_info(self):
        import io, sys, os
        from metalandligandsinfo import find_closest_ligand, print_ligand_info, print_metal_info, load_ligand_data

        name = self.lookup_var.get().strip()
        self.info_display.delete("1.0", tk.END)

        if not name:
            self.info_display.insert(tk.END, "⚠️ Please enter a name.")
            return

        # Load ligand data as RDKit Mol objects
        sdf_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "all_ligands.sdf"))
        data = load_ligand_data(sdf_path)

        buf = io.StringIO()
        old_stdout = sys.stdout
        sys.stdout = buf

        try:
            if name in data:
                print_ligand_info(data[name], name)
            elif any(m.get("symbol", "").upper() == name.upper() for m in METALS):
                print_metal_info(name)
            else:
                suggestion = find_closest_ligand(name, data)
                if suggestion:
                    print_ligand_info(data[suggestion], suggestion)
                else:
                    print("❌ Name not recognized as ligand or metal.")
        except Exception as e:
            print(f"❌ Error during lookup: {e}")
        finally:
            sys.stdout = old_stdout

        self.info_display.insert(tk.END, buf.getvalue())



    def add_ligand(self):
        ligand = self.selected_ligand.get()
        try:
            count = int(self.ligand_count_var.get())
        except ValueError:
            count = 1
        if ligand:
            self.lig_list_accum[ligand] = self.lig_list_accum.get(ligand, 0) + count
            display = ", ".join(f"{k} × {v}" for k, v in self.lig_list_accum.items())
            self.lig_list_display.set(display)


    def draw_complex(self):
        try:
            metal = self.metal_var.get()
            geometry = self.geom_var.get().replace(' ', '_')
            ligand_counts = {
                name.strip(): int(cnt)
                for part in self.ligand_var.get().split(',')
                if ':' in part for name, cnt in [part.split(':')]
            }
            
            bond_length = float(self.bond_length_var.get() or "1.0")
            img = create_complex_from_ligand_dict(metal, ligand_counts, geometry, bond_length=bond_length)
            resized_img = img.resize((400, 400), Image.Resampling.LANCZOS)
            self.photo = ImageTk.PhotoImage(resized_img)
            self.canvas.update_idletasks()
            w, h = self.canvas.winfo_width(), self.canvas.winfo_height()
            self.canvas.delete("all")
            self.canvas.create_image(w/2, h/2, image=self.photo, anchor='center')

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def compute_orbital(self):
        import io, sys, os
        old_stdout = sys.stdout
        buf = io.StringIO()

        try:
            metal_sym = self.metal_orb_var.get()
            ox_state = int(self.oxstate_var.get())
            lig_counts = self.lig_list_accum
            geom = self.geom_orb_var.get()

            # Correct path: all_ligands.sdf is in the repo root
            sdf_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "all_ligands.sdf"))
            data = orbital.load_ligand_data(sdf_path)
            metal = orbital.get_metal_data(metal_sym)

            deltas, charges = [], []
            for name, count in lig_counts.items():
                info = data.get(name) or data.get(orbital.find_closest_ligand(name, data)) or {"Delta": 20000, "Charge": 0}
                deltas.extend([info["Delta"]] * count)
                charges.extend([info["Charge"]] * count)

            d_elec = orbital.get_d_electron_count(metal["atomic_number"], ox_state)
            avg_delta = sum(deltas) / len(deltas)
            spin = orbital.predict_spin_state(avg_delta, orbital.get_pairing_energy(metal["block"]))

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
