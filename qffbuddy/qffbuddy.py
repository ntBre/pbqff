import tkinter as tk
from tkinter import ttk
import argparse
import tkinter.filedialog as fd
import subprocess
import sys
import json
import io
import re
import os
import tkinter.messagebox as msg
import tkinter.scrolledtext as st
from idlelib.tooltip import Hovertip

parser = argparse.ArgumentParser(
    prog="qffbuddy",
    description="helper for generating pbqff input files",
)
parser.add_argument("-p", "--pbqff", help="path to pbqff executable")
parser.add_argument(
    "infile", help="start the configuration based on infile", default=None, nargs="?"
)
args = parser.parse_args()

RE_MOLPRO = re.compile("molpro", re.IGNORECASE)

MOLPRO_F12TZ = """memory,1,g
gthresh,energy=1.d-12,zero=1.d-22,oneint=1.d-22,twoint=1.d-22;
gthresh,optgrad=1.d-8,optstep=1.d-8;
nocompress;

geometry={
{{.geom}}
basis={
default,cc-pVTZ-f12
}
set,charge={{.charge}}
set,spin=0
hf,accuracy=16,energy=1.0d-10
{CCSD(T)-F12,thrden=1.0d-12,thrvar=1.0d-10}
{optg,grms=1.d-8,srms=1.d-8}

pbqff=energy(2)
show[1,f20.12],pbqff
"""

MOLPRO_F12TZCCR = """memory,1,g
gthresh,energy=1.d-12,zero=1.d-22,oneint=1.d-22,twoint=1.d-22;
gthresh,optgrad=1.d-8,optstep=1.d-8;
nocompress;

geometry={
{{.geom}}
basis={
default,cc-pCVTZ-f12
}
set,charge={{.charge}}
set,spin=0
hf,accuracy=16,energy=1.0d-10
{CCSD(T)-F12,thrden=1.0d-12,thrvar=1.0d-10,nocheck;core}
{optg,grms=1.d-8,srms=1.d-8}
etz=energy

basis=cc-pvtz-dk
hf,accuracy=16,energy=1.0d-10
{CCSD(T),thrden=1.0d-12,thrvar=1.0d-10,nocheck;}
edk=energy

basis=cc-pvtz-dk
dkroll=1
hf,accuracy=16,energy=1.0d-10
{CCSD(T),thrden=1.0d-12,thrvar=1.0d-10,nocheck;}
edkr=energy

pbqff=etz(2)+edkr-edk
show[1,f20.12],pbqff
"""

MOPAC_TEMPLATE = "scfcrt=1.D-21 aux(precision=14 comp xp xs xw) PM6 SINGLET THREADS=1"

ATOMS = [
    "X",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
]


TEMPLATES = [MOLPRO_F12TZ]


def make_radio_buttons(pairs, var, parent, **kwargs):
    frame = tk.Frame(parent)
    col = 0
    for text, value in pairs:
        ttk.Radiobutton(frame, text=text, variable=var, value=value, **kwargs).grid(
            column=col, row=0, padx=5
        )
        col += 1
    return frame


class Application(ttk.Frame):
    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        parent.title("qffbuddy")
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(0, weight=1)
        self.parent = parent

        self.geometry_input()

        self.optimize = tk.BooleanVar()
        check = ttk.Checkbutton(
            self,
            text="does it need to be optimized?",
            variable=self.optimize,
        ).grid(column=1, row=4, sticky=tk.W)

        ttk.Label(self, text="Charge").grid(column=1, row=5, padx=10, sticky="E")
        self.charge = tk.IntVar()
        name = ttk.Entry(self, textvariable=self.charge).grid(
            column=2,
            row=5,
        )

        ttk.Label(self, text="Step size in Å").grid(column=1, row=7, sticky="E")
        self.step_size = tk.DoubleVar(value=0.005)
        name = ttk.Entry(self, textvariable=self.step_size).grid(column=2, row=7)

        ttk.Label(self, text="Sleep interval in sec").grid(column=1, row=9, sticky="E")
        self.sleep_int = tk.IntVar(value=2)
        ttk.Entry(self, textvariable=self.sleep_int).grid(column=2, row=9)

        ttk.Label(self, text="Max jobs to submit at once").grid(
            column=1, row=11, sticky="E"
        )
        self.job_limit = tk.IntVar(value=1024)
        ttk.Entry(self, textvariable=self.job_limit).grid(column=2, row=11)

        ttk.Label(self, text="Jobs per chunk").grid(column=1, row=13, sticky="E")
        self.chunk_size = tk.IntVar(value=1)
        name = ttk.Entry(self, textvariable=self.chunk_size).grid(column=2, row=13)

        ttk.Label(self, text="Checkpoint interval (0 to disable)").grid(
            column=1, row=14, stick="E"
        )
        self.check_int = tk.IntVar(value=100)
        ttk.Entry(self, textvariable=self.check_int).grid(column=2, row=14)

        self.coord_select()

        self.chem_prog()

        row = self.queue_system()

        self.template_input()

        self.hybrid_input()

        ttk.Label(self, text="Generated filename").grid(
            column=1, sticky=tk.W, row=row + 13
        )
        self.infile = tk.StringVar(value="pbqff.toml")
        ttk.Entry(self, textvariable=self.infile).grid(
            column=2, row=row + 13, sticky=tk.W
        )

        frame = tk.Frame(self)
        button = ttk.Button(frame, text="Generate", command=self.generate)
        button.grid(column=0, row=0, sticky="E", padx=5)

        button = ttk.Button(frame, text="Run", command=self.run)
        button.grid(column=1, row=0, padx=5)

        button = ttk.Button(frame, text="Exit", command=parent.destroy)
        button.grid(column=2, row=0, sticky="W", padx=5)

        frame.grid(column=2, pady=10, sticky="W")

    def geometry_input(self):
        self.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))
        ttk.Label(self, text="enter your geometry in Å:").grid(
            column=1, row=1, sticky=tk.W
        )
        self.geometry = st.ScrolledText(self, width=80, height=10, undo=True)
        self.geometry.grid(column=1, row=2, columnspan=2, sticky=tk.W)

    def coord_select(self):
        "select the type of coordinates to use for the QFF"
        label = ttk.Label(self, text="Coordinate type")
        label.grid(column=1, sticky="E")
        Hovertip(
            label,
            """The coordinate system to use for the QFF. SIC requires a template INTDER input
file named 'intder.in' in the current directory to specify the internal
coordinates. Cartesian displaces directly along the x, y, and z axes. Normal
uses a Cartesian HFF to generate normal coordinates and uses those for the rest
of the QFF.

Currently, finite differences are always used for Cartesians, and an ANPASS
fitting is always used for SICs. An additional checkbox will appear to toggle
between these for Normals.""",
        )
        self.coord_type = tk.StringVar(value="sic")
        frame = make_radio_buttons(
            [("SIC", "sic"), ("Cartesian", "cart"), ("Normal", "normal")],
            self.coord_type,
            self,
            command=self.toggle_normal,
        )
        self.findiff = tk.BooleanVar(value=False)

        frame.grid(column=2, row=15)

    def toggle_normal(self):
        "check if `self.coord_type` is normal and hide/show the findiff prompt accordingly"
        if self.coord_type.get() == "normal":
            self.findiff_button = ttk.Checkbutton(
                self,
                text="finite differences?",
                variable=self.findiff,
            )
            Hovertip(
                self.findiff_button,
                "Compute the force constants directly with finite differences \
instead of performing a fitting with ANPASS",
            )
            self.findiff_button.grid(row=15, column=3, sticky=tk.W)
        else:
            self.findiff.set(False)
            self.findiff_button.grid_remove()

    def chem_prog(self, column=2, row=16):
        "select the chemistry program to use"
        ttk.Label(self, text="Chemistry program").grid(column=1, sticky="E")
        self.program = tk.StringVar(value="molpro")
        f = make_radio_buttons(
            [("Molpro", "molpro"), ("Mopac", "mopac")],
            self.program,
            self,
            command=self.default_template,
        )
        f.grid(column=column, row=row)

    def queue_system(self):
        "select the queuing system to use"
        l = ttk.Label(self, text="Queuing System")
        l.grid(column=1, sticky="E")
        row = l.grid_info()["row"]
        self.queue = tk.StringVar(value="pbs")
        f = make_radio_buttons(
            [("PBS", "pbs"), ("Slurm", "slurm"), ("Local", "local")], self.queue, self
        )
        f.grid(column=2, row=row)

        self.show_queue_template = tk.BooleanVar(value=False)
        b = ttk.Checkbutton(
            self,
            text="custom template?",
            variable=self.show_queue_template,
            command=self.toggle_queue_template,
        )
        Hovertip(
            b,
            """pbqff includes default templates for the supported queue types, but you can also
include a custom template.""",
        )
        b.grid(column=3, row=row)

        return row

    def toggle_queue_template(self):
        return

    def template_input(self):
        ttk.Label(self, text="enter your template input file:").grid(
            column=1, sticky=tk.W
        )
        self.template = st.ScrolledText(self, width=80, height=10, undo=True)
        self.default_template()
        self.template.grid(column=1, columnspan=2, sticky="W")

    def run(self):
        "run pbqff with the current input file"
        res = msg.askyesno(
            "Run pbqff and exit",
            "are you sure? this will overwrite any existing output",
            icon="warning",
        )
        if res:
            self.generate()
            os.system(f"{args.pbqff} -t 8 -o {self.infile.get()} & disown -h")
            self.parent.destroy()

    def default_template(self):
        self.template.delete("1.0", "end")
        if self.program.get() == "molpro":
            self.template.insert("1.0", MOLPRO_F12TZ)
        elif self.program.get() == "mopac":
            self.template.insert("1.0", MOPAC_TEMPLATE)
        else:
            self.template.insert("1.0", "default template")

    def fill_geometry(self, new_value):
        "clear the geometry input box and fill with `new_value`"
        self.geometry.delete("1.0", "end")
        self.geometry.insert("1.0", new_value)

    def fill_template(self, new_value, hybrid=False):
        "clear the template input box and fill with `new_value`"
        if hybrid:
            self.hybrid_template.delete("1.0", "end")
            self.hybrid_template.insert("1.0", new_value)
        else:
            self.template.delete("1.0", "end")
            self.template.insert("1.0", new_value)

    def hybrid_input(self):
        self.is_hybrid = tk.BooleanVar()
        butt = ttk.Checkbutton(
            self,
            text="do you want to use a hybrid method?",
            variable=self.is_hybrid,
            command=self.toggle_hybrid,
        )
        butt.grid(column=1, stick="W")
        Hovertip(
            butt,
            """Use the main template input file above for the harmonic force constants and a
different template for the cubic and quartic force constants. Currently this
option is only supported for finite difference normal coordinate QFFs.""",
        )
        self.hybrid_label = ttk.Label(
            self, text="enter your template input file for cubics and quartics:"
        )
        self.hybrid_row = butt.grid_info()["row"] + 1
        self.hybrid_template = st.ScrolledText(self, width=80, height=10, undo=True)

    def toggle_hybrid(self):
        if self.is_hybrid.get():
            self.hybrid_label.grid(row=self.hybrid_row, column=1, sticky=tk.W)
            self.hybrid_template.grid(
                row=self.hybrid_row + 1, column=1, sticky=tk.W, columnspan=2
            )
        else:
            self.hybrid_label.grid_remove()
            self.hybrid_template.grid_remove()

    def generate(self):
        with open(self.infile.get(), "w") as out:
            if self.optimize.get():
                opt = "true"
            else:
                opt = "false"
            out.write(
                f"""geometry = \"\"\"
{self.geometry.get('1.0', 'end').strip()}
\"\"\"
optimize = {opt}
charge = {self.charge.get()}
step_size = {self.step_size.get()}
sleep_int = {self.sleep_int.get()}
job_limit = {self.job_limit.get()}
chunk_size = {self.chunk_size.get()}
coord_type = \"{self.coord_type.get()}\"
template = \"\"\"{self.template.get("1.0", "end").strip()}\"\"\"
program = \"{self.program.get()}\"
queue = \"{self.queue.get()}\"
check_int = {self.check_int.get()}
"""
            )


class MenuBar(tk.Menu):
    def __init__(self, parent, *args, **kwargs):
        tk.Menu.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.parent["menu"] = self

        self.menu_file = tk.Menu(self)
        self.add_cascade(menu=self.menu_file, label="File")
        self.menu_file.add_command(label="Open Input File", command=self.open_file)
        self.menu_file.add_command(label="Import Geometry", command=self.import_geom)

        self.menu_templates = tk.Menu(self)
        self.add_cascade(menu=self.menu_templates, label="Templates")

        self.main_templates = tk.Menu(self.menu_templates)
        self.menu_templates.add_cascade(menu=self.main_templates, label="Main")

        self.build_template_menus(self.main_templates)

        self.hybrid_templates = tk.Menu(self.menu_templates)
        self.menu_templates.add_cascade(menu=self.hybrid_templates, label="Hybrid")

        self.build_template_menus(self.hybrid_templates, hybrid=True)

    def build_template_menus(self, menu_var, hybrid=False):
        menu_var.add_command(
            label="F12-DZ",
            command=lambda: app.fill_template(MOLPRO_F12TZ.replace("TZ", "DZ"), hybrid),
        )
        menu_var.add_command(
            label="F12-TZ", command=lambda: app.fill_template(MOLPRO_F12TZ, hybrid)
        )
        menu_var.add_command(
            label="F12-DZ-cCR",
            command=lambda: app.fill_template(
                MOLPRO_F12TZCCR.replace("TZ", "DZ"), hybrid
            ),
        )
        menu_var.add_command(
            label="F12-TZ-cCR",
            command=lambda: app.fill_template(MOLPRO_F12TZCCR, hybrid),
        )

    def import_geom(self):
        infile = fd.askopenfile(parent=self.parent)
        if infile is None:
            return
        else:
            infile = infile.name
        with open(infile, "r") as inp:
            s = inp.read()
            if RE_MOLPRO.search(s):
                geom = self.parse_molpro_geom(s)
                if geom is not None:
                    app.fill_geometry(geom)

    def parse_molpro_geom(self, s):
        skip = 0
        in_geom = False
        geom = io.StringIO()
        for line in s.split("\n"):
            if "Current geometry" in line:
                skip = 1
                in_geom = True
            elif skip > 0:
                skip -= 1
            elif in_geom and len(line) == 0:
                return geom.getvalue()
            elif in_geom:
                geom.write(line + "\n")

    def open_file(self):
        infile = fd.askopenfile(parent=self.parent)
        if infile is None:
            return
        else:
            infile = infile.name
        self.parse_infile(infile)

    def parse_infile(self, infile):
        s = subprocess.run([args.pbqff, "-j", infile], capture_output=True)
        if s.returncode != 0:
            sys.exit(f"pbqff failed to parse {infile}: {s.stderr}")
        d = json.loads(s.stdout)

        geom = d["geometry"]
        if "Zmat" in geom:
            app.fill_geometry(geom["Zmat"])
        elif "Xyz" in geom:
            s = io.StringIO()
            for atom in geom["Xyz"]:
                s.write(
                    "%2s%14.10f%14.10f%14.10f\n"
                    % (ATOMS[atom["atomic_number"]], atom["x"], atom["y"], atom["z"])
                )
            app.fill_geometry(s.getvalue())

        app.optimize.set(d["optimize"])
        app.charge.set(d["charge"])
        app.step_size.set(d["step_size"])
        app.coord_type.set(d["coord_type"].lower())
        app.check_int.set(d["check_int"])
        app.program.set(d["program"].lower())
        app.queue.set(d["queue"].lower())
        app.sleep_int.set(d["sleep_int"])
        app.job_limit.set(d["job_limit"])
        app.chunk_size.set(d["chunk_size"])
        app.fill_template(d["template"])


if __name__ == "__main__":
    root = tk.Tk()
    root.option_add("*tearOff", tk.FALSE)
    app = Application(root, padding="3 3 12 12")
    menu = MenuBar(root)
    if args.infile is not None:
        menu.parse_infile(args.infile)
    root.mainloop()
