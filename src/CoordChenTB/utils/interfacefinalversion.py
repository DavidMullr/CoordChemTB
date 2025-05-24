import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from PIL import Image, ImageTk
from complexdrawerfinal import (
    create_complex_from_ligand_dict,
    load_ligands_from_folder,
    LIGAND_FOLDER,
    get_donor_atom_indices,
    predict_4coordinate_geometry
)
from rdkit.Chem import MolToSmiles
import complexorbitalssplitting as orbital
from metals_db import METALS

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


class CoordinationGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Coordination Chemistry Toolkit")

        style = ttk.Style(self)
        style.theme_use('clam')
        style.configure('TNotebook.Tab', background='#444', foreground='white')
        style.map(
            'TNotebook.Tab',
            background=[('selected', '#666')],
            foreground=[('selected', 'white')]
        )
        style.configure('TLabel', background='#333', foreground='white')
        style.configure('TFrame', background='#333')
        style.configure('TEntry', fieldbackground='#555', foreground='white')
        style.configure('TButton', background='#444', foreground='white')

        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill='both', expand=True)

        self.sdf_file_list = []
        self.lig_map = load_ligands_from_folder(LIGAND_FOLDER)
        self.ligand_names = sorted(self.lig_map.keys())
        self.metal_symbols = [m['symbol'] for m in METALS]
        self.lookup_values = self.ligand_names + self.metal_symbols

        self.init_complex_tab()
        self.init_orbital_tab()
        self.init_info_tab()

    def init_complex_tab(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Complex Builder")

        for i in range(2):
            frame.columnconfigure(i, weight=1)

        # Metal selection row
        ttk.Label(frame, text="Metal:").grid(row=0, column=0, sticky='ew', padx=5, pady=5)
        metals = [m['symbol'] for m in METALS]
        self.metal_var = tk.StringVar()
        metal_cb = ttk.Combobox(frame, textvariable=self.metal_var, values=metals, state="readonly")
        metal_cb.grid(row=1, column=0, sticky='ew', padx=5)
        
        # Oxidation state selection
        ttk.Label(frame, text="Oxidation State:").grid(row=0, column=1, sticky='ew', padx=5, pady=5)
        self.oxstate_var = tk.StringVar()
        self.oxstate_cb = ttk.Combobox(frame, textvariable=self.oxstate_var, state="readonly")
        self.oxstate_cb.grid(row=1, column=1, sticky='ew', padx=5)
        
        # Update oxidation states when metal changes
        def update_oxidation_states(*args):
            metal = self.metal_var.get()
            if metal:
                metal_info = next((m for m in METALS if m['symbol'] == metal), None)
                if metal_info:
                    states = metal_info.get('oxidation_states', [2])
                    self.oxstate_cb['values'] = states
                    if states:
                        self.oxstate_var.set(str(states[0]))
        self.metal_var.trace_add('write', update_oxidation_states)

        # Ligands input
        ttk.Label(frame, text="Ligands (Name:Count, e.g. py:2, h2o:4):")\
            .grid(row=2, column=0, columnspan=2, sticky='ew', padx=5, pady=5)
        self.ligand_var = tk.StringVar()
        ttk.Entry(frame, textvariable=self.ligand_var, width=60)\
            .grid(row=3, column=0, columnspan=2, sticky='ew', padx=5)

        self.ligand_sel_vars   = [tk.StringVar() for _ in range(6)]
        self.ligand_count_vars = [tk.StringVar(value="1") for _ in range(6)]

        for i in range(6):
            ttk.Combobox(
                frame,
                textvariable=self.ligand_sel_vars[i],
                values=self.ligand_names,
                state="readonly"
            ).grid(row=4 + i, column=0, sticky="ew", padx=5, pady=2)

            ttk.Entry(
                frame,
                textvariable=self.ligand_count_vars[i],
                width=5
            ).grid(row=4 + i, column=1, sticky="ew", padx=5, pady=2)

        # Bond length and draw button
        ttk.Label(frame, text="Bond Length (for bulky ligands, adjust):")\
            .grid(row=10, column=0, columnspan=2, sticky="ew", padx=5, pady=5)
        self.bond_length_var = tk.StringVar(value="1.0")
        ttk.Entry(frame, textvariable=self.bond_length_var)\
            .grid(row=11, column=0, columnspan=2, sticky="ew", padx=5)

        ttk.Button(frame, text="Draw Complex", command=self.draw_complex)\
            .grid(row=12, column=0, columnspan=2, sticky="ew", pady=10)

        # Canvas for drawing
        self.canvas = tk.Canvas(frame, bg="white")
        self.canvas.grid(row=13, column=0, columnspan=2, sticky="nsew", padx=5, pady=5)
        frame.rowconfigure(13, weight=1)

    def init_orbital_tab(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Orbital Splitting")

        for i in range(4):
            frame.columnconfigure(i, weight=1)

        self.metal_orb_var = tk.StringVar()
        self.orb_oxstate_var = tk.StringVar()
        self.selected_ligand = tk.StringVar()
        self.geom_orb_var = tk.StringVar()
        self.lig_list_accum = {}

        ttk.Label(frame, text="Metal:")\
            .grid(row=1, column=0, columnspan=3, sticky='ew', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.metal_orb_var)\
            .grid(row=2, column=0, columnspan=3, sticky='ew', padx=5)

        ttk.Label(frame, text="Oxidation State of Metal:")\
            .grid(row=3, column=0, columnspan=3, sticky='ew', padx=5, pady=5)
        ttk.Entry(frame, textvariable=self.orb_oxstate_var)\
            .grid(row=4, column=0, columnspan=3, sticky='ew', padx=5)

        ttk.Label(frame, text="Select Ligand:")\
            .grid(row=5, column=0, columnspan=3, sticky='ew', padx=5, pady=5)
        self.ligand_dropdown = ttk.Combobox(
            frame, textvariable=self.selected_ligand, values=self.ligand_names, state="readonly"
        )
        self.ligand_dropdown.grid(row=6, column=0, columnspan=2, sticky='ew', padx=5)
        self.ligand_count_var = tk.StringVar(value="1")
        ttk.Entry(frame, textvariable=self.ligand_count_var, width=5)\
            .grid(row=6, column=2, padx=5)
        ttk.Button(frame, text="Add Ligand", command=self.add_ligand)\
            .grid(row=6, column=3, sticky='ew', padx=5)

        ttk.Label(frame, text="Ligands Selected:")\
            .grid(row=7, column=0, columnspan=4, sticky='ew', padx=5)
        self.lig_list_display = tk.StringVar()
        ttk.Entry(frame, textvariable=self.lig_list_display, width=60)\
            .grid(row=8, column=0, columnspan=4, sticky='ew', padx=5)

        ttk.Label(frame, text="Geometry:")\
            .grid(row=9, column=0, columnspan=4, sticky='ew', padx=5, pady=5)
        ttk.Combobox(
            frame,
            textvariable=self.geom_orb_var,
            values=["octahedral", "tetrahedral", "square planar"],
            state="readonly"
        ).grid(row=10, column=0, columnspan=4, sticky='ew', padx=5)

        ttk.Button(frame, text="Compute Orbital Splitting", command=self.compute_orbital)\
            .grid(row=11, column=0, columnspan=4, sticky='ew', pady=10)

    def init_info_tab(self):
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Metal & Ligand Info")

        for i in range(2):
            frame.columnconfigure(i, weight=1)

        self.lookup_var = tk.StringVar()
        ttk.Label(frame, text="Select Metal or Ligand:")\
            .grid(row=0, column=0, columnspan=2, sticky='ew', pady=10)
        ttk.Combobox(
            frame,
            textvariable=self.lookup_var,
            values=self.lookup_values,
            state="readonly"
        ).grid(row=1, column=0, columnspan=2, sticky='ew', padx=10)

        ttk.Button(frame, text="Lookup Info", command=self.lookup_info)\
            .grid(row=2, column=0, columnspan=2, sticky='ew', padx=20, pady=10)

        self.info_display = tk.Text(frame, height=20, bg='#3a3a3a', fg='white', wrap='word')
        self.info_display.grid(row=3, column=0, columnspan=2, sticky='nsew', padx=10, pady=10)
        frame.rowconfigure(3, weight=1)

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
            # 1. Get metal and oxidation state
            metal = self.metal_var.get()
            ox_state = int(self.oxstate_var.get() or "0")

            # 2. Build ligand_counts dict
            ligand_counts = {}
            for sel_var, cnt_var in zip(self.ligand_sel_vars, self.ligand_count_vars):
                name = sel_var.get()
                if name:
                    try:
                        cnt = int(cnt_var.get())
                    except ValueError:
                        cnt = 1
                    denticity = get_donor_atom_indices(
                        load_ligands_from_folder(LIGAND_FOLDER)[name]
                    )
                    max_allowed = 6 if len(denticity) == 1 else 3
                    if cnt > max_allowed:
                        cnt = max_allowed
                        cnt_var.set(str(max_allowed))
                    ligand_counts[name] = ligand_counts.get(name, 0) + cnt

            # 3. Determine geometry
            lig_map = load_ligands_from_folder(LIGAND_FOLDER)
            total_sites = sum(
                len(get_donor_atom_indices(lig_map[name])) * cnt
                for name, cnt in ligand_counts.items()
            )
            if total_sites == 4:
                metal_smiles = f"[{metal}]"
                ligand_smiles_list = [
                    MolToSmiles(lig_map[name]) for name, cnt in ligand_counts.items() for _ in range(cnt)
                ]
                geometry = predict_4coordinate_geometry(metal_smiles, ligand_smiles_list)
            elif total_sites == 6:
                geometry = "octahedral"
            else:
                raise ValueError(f"Unsupported total donor sites: {total_sites}")

            # 4. Draw the complex
            bond_length = float(self.bond_length_var.get() or "1.0")
            img = create_complex_from_ligand_dict(
                metal,
                ligand_counts,
                geometry,
                bond_length=bond_length,
                oxidation_state=ox_state
            )

            # 5. Resize & display
            resized_img = img.resize((400, 400), Image.Resampling.LANCZOS)
            self.photo = ImageTk.PhotoImage(resized_img)
            self.canvas.delete("all")
            w, h = self.canvas.winfo_width(), self.canvas.winfo_height()
            x = (w - 400) // 2
            y = (h - 400) // 2
            self.canvas.create_image(x, y, anchor='nw', image=self.photo)

        except Exception as e:
            messagebox.showerror("Error", str(e))

    def lookup_info(self):
        import io, sys, os
        from metalandligandsinfo import find_closest_ligand, print_ligand_info, print_metal_info, load_ligand_data

        name = self.lookup_var.get().strip()
        self.info_display.delete("1.0", tk.END)

        if not name:
            self.info_display.insert(tk.END, "⚠️ Please enter a name.")
            return

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

    def compute_orbital(self):
        import io, sys, os
        old_stdout = sys.stdout
        buf = io.StringIO()

        try:
            metal_sym = self.metal_orb_var.get()
            ox_state = int(self.orb_oxstate_var.get())
            lig_counts = self.lig_list_accum
            geom = self.geom_orb_var.get()

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
