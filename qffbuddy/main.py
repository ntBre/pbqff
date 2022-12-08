from tkinter import *
from tkinter import ttk


root = Tk()
root.title("qffbuddy")

parent = ttk.Frame(root, padding="3 3 12 12")
parent.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

ttk.Label(parent, text="enter your geometry in Å:").grid(column=3, row=1, sticky=W)
geom = Text(parent, width=40, height=10)
geom.grid(column=3, row=2)

optimize = BooleanVar()
check = ttk.Checkbutton(
    parent,
    text="does it need to be optimized?",
    variable=optimize,
    onvalue="metric",
    offvalue="imperial",
).grid(column=3, row=4)

ttk.Label(parent, text="Charge").grid(column=3, row=5)
charge = IntVar()
name = ttk.Entry(parent, textvariable=charge).grid(column=3, row=6)

ttk.Label(parent, text="Step size in Å").grid(column=3, row=7)
step_size = DoubleVar(value=0.005)
name = ttk.Entry(parent, textvariable=step_size).grid(column=3, row=8)

ttk.Label(parent, text="Sleep interval in sec").grid(column=3, row=9)
sleep_int = IntVar(value=2)
ttk.Entry(parent, textvariable=sleep_int).grid(column=3, row=10)

ttk.Label(parent, text="Max jobs to submit at once").grid(column=3, row=11)
max_jobs = IntVar(value=1024)
ttk.Entry(parent, textvariable=max_jobs).grid(column=3, row=12)

ttk.Label(parent, text="Jobs per chunk").grid(column=3, row=13)
chunk_size = IntVar(value=1)
name = ttk.Entry(parent, textvariable=chunk_size).grid(column=3, row=14)


def make_radio_buttons(pairs, var, parent):
    for (text, value) in pairs:
        ttk.Radiobutton(parent, text=text, variable=var, value=value).grid(column=3)


ttk.Label(parent, text="Coordinate type").grid(column=3)
choices = [("sic", "sic"), ("cart", "cart"), ("normal", "normal")]
coord_type = StringVar()
make_radio_buttons(choices, coord_type, parent)

ttk.Label(parent, text="Chemistry program").grid(column=3)
program = StringVar(parent)
make_radio_buttons([("Molpro", "molpro"), ("Mopac", "mopac")], program, parent)

# TODO callback on program that updates the default template
ttk.Label(parent, text="enter your template input file").grid(column=3, sticky=W)
template = Text(parent, width=40, height=10)
template.grid(column=3)
template.insert("1.0", "example template")

ttk.Label(parent, text="Queuing System").grid(column=3)
queue = StringVar(parent)
make_radio_buttons([("PBS", "pbs"), ("Slurm", "slurm")], queue, parent)

root.mainloop()
