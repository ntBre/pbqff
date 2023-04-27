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


class MyHovertip(Hovertip):
    def __init__(self, anchor_widget, text, hover_delay=500):
        super().__init__(anchor_widget, text, hover_delay)


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
    def __init__(self, parent, root, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.root = root
        self.root.title("qffbuddy")
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        self.parent = parent

        self.grid(column=0, row=0, sticky=(tk.N, tk.W))

        self.main_panel = tk.Frame(self)
        self.main_panel.grid(column=1, sticky="nw")

        self.geometry_input()

        self.optimize_prompt()

        self.charge_input()

        self.step_size_input()
        self.sleep_int_input()
        self.job_limit_input()
        self.chunk_size_input()
        self.check_int_input()

        self.coord_type_input()

        self.program_input()

        row = self.queue_input()

        self.queue_template_panel.grid(column=1, sticky=tk.W)

        self.template_input()

        self.hybrid_input()

        self.hybrid_panel.grid(column=1, sticky=tk.W)

        self.build_lower_panel(row)

        self.lower_panel.grid(column=1)

    def build_lower_panel(self, row):
        """build the lower panel of the main interface, including the
        generated filename input andthe final buttons"""
        self.lower_panel = tk.Frame(self)
        ttk.Label(self.lower_panel, text="Generated filename").grid(
            column=0, sticky=tk.W, row=row + 13
        )
        self.infile = tk.StringVar(value="pbqff.toml")
        ttk.Entry(self.lower_panel, textvariable=self.infile).grid(
            column=2, row=row + 13, sticky=tk.W
        )
        frame = tk.Frame(self.lower_panel)
        button = ttk.Button(frame, text="Generate", command=self.generate)
        button.grid(column=0, row=0, sticky="E", padx=5)

        button = ttk.Button(frame, text="Run", command=self.run)
        button.grid(column=1, row=0, padx=5)

        button = ttk.Button(frame, text="Exit", command=self.root.destroy)
        button.grid(column=2, row=0, sticky="W", padx=5)

        frame.grid(column=2, pady=10, sticky="W")

    def step_size_input(self):
        l = ttk.Label(self.main_panel, text="Step size in Å")
        l.grid(column=1, row=7, sticky="E")
        MyHovertip(l, """Displacement size to use in the QFF""")
        self.step_size = tk.DoubleVar(value=0.005)
        name = ttk.Entry(self.main_panel, textvariable=self.step_size).grid(
            column=2, row=7
        )

    def sleep_int_input(self):
        l = ttk.Label(self.main_panel, text="Sleep interval in sec")
        l.grid(column=1, row=9, sticky="E")
        MyHovertip(l, "pbqff will wait this long before polling the running jobs again")
        self.sleep_int = tk.IntVar(value=2)
        ttk.Entry(self.main_panel, textvariable=self.sleep_int).grid(column=2, row=9)

    def job_limit_input(self):
        l = ttk.Label(self.main_panel, text="Max jobs to submit at once")
        l.grid(column=1, row=11, sticky="E")
        MyHovertip(
            l,
            """The maximum number of jobs to submit at one time. "Jobs" here refers to pbqff
jobs, not queue system jobs. The number of queued jobs will be job_limit
(this value) divided by chunk_size (Jobs per chunk).""",
        )
        self.job_limit = tk.IntVar(value=1024)
        ttk.Entry(self.main_panel, textvariable=self.job_limit).grid(column=2, row=11)

    def chunk_size_input(self):
        l = ttk.Label(self.main_panel, text="Jobs per chunk")
        l.grid(column=1, row=13, sticky="E")
        MyHovertip(
            l, "The number of pbqff jobs to group into a single queue submission"
        )
        self.chunk_size = tk.IntVar(value=1)
        name = ttk.Entry(self.main_panel, textvariable=self.chunk_size).grid(
            column=2, row=13
        )

    def check_int_input(self):
        l = ttk.Label(self.main_panel, text="Checkpoint interval (0 to disable)")
        l.grid(column=1, row=14, stick="E")
        MyHovertip(l, "The number of polling cycles between checkpoints")
        self.check_int = tk.IntVar(value=100)
        ttk.Entry(self.main_panel, textvariable=self.check_int).grid(column=2, row=14)

    def charge_input(self):
        l = ttk.Label(self.main_panel, text="Charge")
        MyHovertip(l, """Molecular charge""")
        l.grid(column=1, row=5, padx=10, sticky="E")
        self.charge = tk.IntVar()
        c = ttk.Entry(self.main_panel, textvariable=self.charge)
        c.grid(
            column=2,
            row=5,
        )

    def optimize_prompt(self):
        self.optimize = tk.BooleanVar()
        check = ttk.Checkbutton(
            self.main_panel,
            text="does it need to be optimized?",
            variable=self.optimize,
        )
        check.grid(column=1, row=4, sticky=tk.W)
        MyHovertip(
            check,
            """If checked, pbqff will optimize the geometry in the selected quantum chemistry
program.""",
        )

    def geometry_input(self):
        l = ttk.Label(self.main_panel, text="enter your geometry in Å:")
        l.grid(column=1, row=1, sticky=tk.W)
        MyHovertip(
            l,
            """Input your molecular geometry in a format recognized by your quantum chemistry
program of choice. pbqff recognizes Z-matrices and XYZ geometries in general.""",
        )
        self.geometry = st.ScrolledText(self.main_panel, width=80, height=10, undo=True)
        self.geometry.grid(column=1, row=2, columnspan=2, sticky=tk.W)

    def coord_type_input(self):
        "select the type of coordinates to use for the QFF"
        label = ttk.Label(self.main_panel, text="Coordinate type")
        label.grid(column=1, sticky="E")
        MyHovertip(
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
        self.findiff_frame = make_radio_buttons(
            [("SIC", "sic"), ("Cartesian", "cart"), ("Normal", "normal")],
            self.coord_type,
            self.main_panel,
            command=self.toggle_normal,
        )
        self.findiff = tk.BooleanVar(value=False)
        self.findiff_button = ttk.Checkbutton(
            self.findiff_frame,
            text="finite differences?",
            variable=self.findiff,
        )
        MyHovertip(
            self.findiff_button,
            "Compute the force constants directly with finite differences \
instead of performing a fitting with ANPASS",
        )

        self.findiff_frame.grid(column=2, row=15)

    def set_findiff(self, value):
        "set `self.findiff` to `value` and update the display accordingly"
        self.findiff.set(value)
        if value:
            self.findiff_button.grid(row=0, column=4, sticky=tk.W)
        else:
            self.findiff_button.grid_remove()

    def toggle_normal(self):
        "check if `self.coord_type` is normal and hide/show the findiff prompt accordingly"
        if self.coord_type.get() == "normal":
            self.findiff_button.grid(row=0, column=4, sticky=tk.W)
        else:
            self.findiff.set(False)
            self.findiff_button.grid_remove()

    def program_input(self, column=2, row=16):
        "select the chemistry program to use"
        l = ttk.Label(self.main_panel, text="Chemistry program")
        l.grid(column=1, sticky="E")
        MyHovertip(
            l,
            """The chemistry program to use for geometry optimizations and single-point energy
computations""",
        )
        self.program = tk.StringVar(value="molpro")
        f = make_radio_buttons(
            [("Molpro", "molpro"), ("Mopac", "mopac")],
            self.program,
            self.main_panel,
            command=self.default_template,
        )
        f.grid(column=column, row=row)

    def queue_input(self):
        "select the queuing system to use"
        l = ttk.Label(self.main_panel, text="Queuing System")
        MyHovertip(
            l,
            """The queuing system to use. PBS and Slurm will submit jobs to the corresponding
queues and monitor them using the qstat and squeue utilities. A Local queue will
run jobs directly on the current computer using the bash shell as a mock queue.""",
        )
        l.grid(column=1, sticky="E")
        row = l.grid_info()["row"]
        self.queue = tk.StringVar(value="pbs")
        f = make_radio_buttons(
            [("PBS", "pbs"), ("Slurm", "slurm"), ("Local", "local")],
            self.queue,
            self.main_panel,
        )
        f.grid(column=2, row=row)

        self.queue_template_panel = tk.Frame(self)

        self.show_queue_template = tk.BooleanVar(value=False)
        b = ttk.Checkbutton(
            self.queue_template_panel,
            text="custom queue template?",
            variable=self.show_queue_template,
            command=self.toggle_queue_template,
        )
        MyHovertip(
            b,
            """pbqff includes default templates for the supported queue types, but you can also
include a custom template.""",
        )
        b.grid(column=0, sticky=tk.W)

        self.queue_template_label = ttk.Label(
            self.queue_template_panel, text="queue template:"
        )
        self.queue_template = st.ScrolledText(
            self.queue_template_panel, width=80, height=10, undo=True
        )

        return row

    def toggle_queue_template(self):
        if self.show_queue_template.get():
            self.queue_template_label.grid(column=0, sticky=tk.W)
            self.queue_template.grid(column=0, sticky=tk.W, columnspan=2)
        else:
            self.queue_template_label.grid_remove()
            self.queue_template.grid_remove()

    def template_input(self):
        l = ttk.Label(self.main_panel, text="enter your template input file:")
        l.grid(column=1, sticky=tk.W)
        MyHovertip(
            l,
            """The template input file for the quantum chemistry program of choice. Some
programs support special string substitutions. Molpro, for example, supports
{{.geom}} replacement for the geometry and {{.charge}} for the molecular charge.
Mopac expects the full template to be given here, verbatim. See the Templates
menu for some examples.""",
        )
        self.template = st.ScrolledText(self.main_panel, width=80, height=10, undo=True)
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
            self.root.destroy()

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

    def fill_queue_template(self, new_value):
        self.queue_template.delete("1.0", "end")
        self.queue_template.insert("1.0", new_value)

    def hybrid_input(self):
        self.is_hybrid = tk.BooleanVar()
        self.hybrid_panel = tk.Frame(self)
        butt = ttk.Checkbutton(
            self.hybrid_panel,
            text="do you want to use a hybrid method?",
            variable=self.is_hybrid,
            command=self.toggle_hybrid,
        )
        butt.grid(column=1, sticky=tk.W)
        MyHovertip(
            butt,
            """Use the main template input file above for the harmonic force constants and a
different template for the cubic and quartic force constants. Currently this
option is only supported for finite difference normal coordinate QFFs.""",
        )
        self.hybrid_label = ttk.Label(
            self.hybrid_panel,
            text="enter your template input file for cubics and quartics:",
        )
        self.hybrid_row = butt.grid_info()["row"] + 1
        self.hybrid_template = st.ScrolledText(
            self.hybrid_panel, width=80, height=10, undo=True
        )

    def toggle_hybrid(self):
        if self.is_hybrid.get():
            self.hybrid_label.grid(column=1, sticky=tk.W)
            self.hybrid_template.grid(column=1, sticky=tk.W, columnspan=2)
        else:
            self.hybrid_label.grid_remove()
            self.hybrid_template.grid_remove()

    def generate(self):
        with open(self.infile.get(), "w") as out:
            opt = str(self.optimize.get()).lower()
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
findiff = {str(self.findiff.get()).lower()}
template = \"\"\"{self.template.get("1.0", "end").strip()}\"\"\"
program = \"{self.program.get()}\"
queue = \"{self.queue.get()}\"
check_int = {self.check_int.get()}
"""
            )
            if self.show_queue_template.get():
                out.write(
                    f"""
queue_template = \"\"\"{self.queue_template.get("1.0", "end")}\"\"\"
"""
                )
            if self.is_hybrid.get():
                out.write(
                    f"""
hybrid_template = \"\"\"{self.hybrid_template.get("1.0", "end")}\"\"\"
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

        self.queue_templates = tk.Menu(self.menu_templates)
        self.menu_templates.add_cascade(menu=self.queue_templates, label="Queue")

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

        app.set_findiff(d["findiff"] is not None and d["findiff"])
        app.toggle_normal()

        it = d["hybrid_template"]
        if it is not None and it != "":
            app.is_hybrid.set(True)
            app.toggle_hybrid()
            app.fill_template(it, hybrid=True)

        it = d["queue_template"]
        if it is not None and it != "":
            app.show_queue_template.set(True)
            app.toggle_queue_template()
            app.fill_queue_template(it)


if __name__ == "__main__":
    root = tk.Tk()
    root.option_add("*tearOff", tk.FALSE)
    # from https://stackoverflow.com/a/3092341/12935407
    canvas = tk.Canvas(root)
    app = Application(canvas, root, padding="3 3 12 12")
    vsb = tk.Scrollbar(root, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=vsb.set)

    vsb.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)
    canvas.create_window((0, 0), window=app, anchor="nw")

    app.bind(
        "<Configure>", lambda event: canvas.configure(scrollregion=canvas.bbox("all"))
    )
    canvas.bind("<Button-4>", lambda event: canvas.yview_scroll(-1, "units"))
    canvas.bind("<Button-5>", lambda event: canvas.yview_scroll(1, "units"))

    menu = MenuBar(root)
    if args.infile is not None:
        menu.parse_infile(args.infile)
    root.mainloop()
