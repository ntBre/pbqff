import tkinter as tk
from tkinter import ttk


def make_radio_buttons(pairs, var, parent):
    for (text, value) in pairs:
        ttk.Radiobutton(parent, text=text, variable=var, value=value).grid(column=3)


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
            [("Molpro", "molpro"), ("Mopac", "mopac")], self.program, self
        )

        # TODO callback on program that updates the default template
        ttk.Label(self, text="enter your template input file").grid(
            column=3, sticky=tk.W
        )
        self.template = tk.Text(self, width=40, height=10)
        self.template.grid(column=3)
        self.template.insert("1.0", "example template")

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
coord_type = {self.coord_type.get()}
template = \"\"\"
{self.template.get("1.0", "end").strip()}
\"\"\"
program = {self.program.get()}
queue = {self.queue.get()}
"""
            )


if __name__ == "__main__":
    root = tk.Tk()
    Application(root, padding="3 3 12 12")
    root.mainloop()
