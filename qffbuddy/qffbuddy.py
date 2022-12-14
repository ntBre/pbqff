import tkinter as tk
from tkinter import ttk
import argparse
import tkinter.filedialog as fd
import subprocess
import sys

parser = argparse.ArgumentParser(
    prog="qffbuddy",
    description="helper for generating pbqff input files",
)
parser.add_argument("-p", "--pbqff", help="path to pbqff executable")
args = parser.parse_args()

MOLPRO_TEMPLATE = """***,default f12-tz molpro template
memory,1,g
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

MOPAC_TEMPLATE = "scfcrt=1.D-21 aux(precision=14) PM6 SINGLET THREADS=1"


TEMPLATES = [MOLPRO_TEMPLATE]


def make_radio_buttons(pairs, var, parent, **kwargs):
    for (text, value) in pairs:
        ttk.Radiobutton(parent, text=text, variable=var, value=value, **kwargs).grid(
            column=3
        )


class Application(ttk.Frame):
    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        parent.title("qffbuddy")
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(0, weight=1)
        self.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))
        ttk.Label(self, text="enter your geometry in Å:").grid(
            column=3, row=1, sticky=tk.W
        )
        self.geom = tk.Text(self, width=40, height=10)
        self.geom.grid(column=3, row=2)

        self.optimize = tk.BooleanVar()
        check = ttk.Checkbutton(
            self,
            text="does it need to be optimized?",
            variable=self.optimize,
            onvalue="metric",
            offvalue="imperial",
        ).grid(column=3, row=4)

        ttk.Label(self, text="Charge").grid(column=3, row=5)
        self.charge = tk.IntVar()
        name = ttk.Entry(self, textvariable=self.charge).grid(column=3, row=6)

        ttk.Label(self, text="Step size in Å").grid(column=3, row=7)
        self.step_size = tk.DoubleVar(value=0.005)
        name = ttk.Entry(self, textvariable=self.step_size).grid(column=3, row=8)

        ttk.Label(self, text="Sleep interval in sec").grid(column=3, row=9)
        self.sleep_int = tk.IntVar(value=2)
        ttk.Entry(self, textvariable=self.sleep_int).grid(column=3, row=10)

        ttk.Label(self, text="Max jobs to submit at once").grid(column=3, row=11)
        self.job_limit = tk.IntVar(value=1024)
        ttk.Entry(self, textvariable=self.job_limit).grid(column=3, row=12)

        ttk.Label(self, text="Jobs per chunk").grid(column=3, row=13)
        self.chunk_size = tk.IntVar(value=1)
        name = ttk.Entry(self, textvariable=self.chunk_size).grid(column=3, row=14)

        ttk.Label(self, text="Coordinate type").grid(column=3)
        choices = [("sic", "sic"), ("cart", "cart"), ("normal", "normal")]
        self.coord_type = tk.StringVar(value="sic")
        make_radio_buttons(choices, self.coord_type, self)

        ttk.Label(self, text="Chemistry program").grid(column=3)
        self.program = tk.StringVar(value="molpro")
        make_radio_buttons(
            [("Molpro", "molpro"), ("Mopac", "mopac")],
            self.program,
            self,
            command=self.default_template,
        )

        ttk.Label(self, text="enter your template input file").grid(
            column=3, sticky=tk.W
        )
        self.template = tk.Text(self, width=80, height=10)
        self.default_template()
        self.template.grid(column=3)

        ttk.Label(self, text="Queuing System").grid(column=3)
        self.queue = tk.StringVar(value="pbs")
        make_radio_buttons([("PBS", "pbs"), ("Slurm", "slurm")], self.queue, self)

        ttk.Label(self, text="Generated filename").grid(column=3)
        self.infile = tk.StringVar(value="pbqff.toml")
        ttk.Entry(self, textvariable=self.infile).grid(column=3)

        button = ttk.Button(self, text="Generate", command=self.generate)
        button.grid(column=3)

        button = ttk.Button(self, text="exit", command=parent.destroy)
        button.grid(column=3)

    def default_template(self):
        self.template.delete("1.0", "end")
        if self.program.get() == "molpro":
            self.template.insert("1.0", MOLPRO_TEMPLATE)
        elif self.program.get() == "mopac":
            self.template.insert("1.0", MOPAC_TEMPLATE)
        else:
            self.template.insert("1.0", "default template")

    def generate(self):
        with open(self.infile.get(), "w") as out:
            if self.optimize.get():
                opt = "true"
            else:
                opt = "false"
            out.write(
                f"""geometry = \"\"\"
{self.geom.get('1.0', 'end').strip()}
\"\"\"
optimize = {opt}
charge = {self.charge.get()}
step_size = {self.step_size.get()}
sleep_int = {self.sleep_int.get()}
job_limit = {self.job_limit.get()}
chunk_size = {self.chunk_size.get()}
coord_type = \"{self.coord_type.get()}\"
template = \"\"\"
{self.template.get("1.0", "end").strip()}
\"\"\"
program = \"{self.program.get()}\"
queue = \"{self.queue.get()}\"
"""
            )


class MenuBar(tk.Menu):
    def __init__(self, parent, *args, **kwargs):
        self.parent = parent
        tk.Menu.__init__(self, parent, *args, **kwargs)
        self.menu_file = tk.Menu(self)
        self.add_cascade(menu=self.menu_file, label="File")
        parent["menu"] = self
        self.menu_file.add_command(label="Open...", command=self.open_file)

    # TODO need to come up with a Config class that I can serialize and
    # deserialize from file. I mean I already have one, but I need it in python
    def open_file(self):
        infile = fd.askopenfile(parent=self.parent).name
        s = subprocess.run([args.pbqff, "-j", infile], capture_output=True)
        if s.returncode != 0:
            sys.exit(f"pbqff failed to parse {infile}: {s.stderr}")
        print(f"you tried to open '{infile}' with output:")
        print(s.stdout)


if __name__ == "__main__":
    root = tk.Tk()
    root.option_add("*tearOff", tk.FALSE)
    Application(root, padding="3 3 12 12")
    MenuBar(root)
    root.mainloop()
