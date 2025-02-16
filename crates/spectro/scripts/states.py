#!/usr/bin/python

from sys import argv
from dataclasses import dataclass


@dataclass
class State:
    non_deg: [int]
    deg_vt: [int]
    deg_vl: [int]


# the spectro output file to read
infile = argv[1]

# the full number of vibrations, not counting degeneracies
nvib = int(argv[2])

search = False
states = {}
with open(infile, "r") as infile:
    state = 0
    for line in infile:
        if "VIBRATIONAL ENERGY AND PROPERTIES" in line:
            search = True
        elif search and "NON-DEG (Vs)" in line:
            sp = line.split()
            state = int(sp[0])
            if state not in states:
                states[state] = State([], [], [])
            states[state].non_deg.extend([int(x) for x in sp[6:]])
        elif search and "DEGEN   (Vt)" in line:
            sp = line.split()
            states[state].deg_vt.extend([int(x) for x in sp[3:]])
        elif search and "DEGEN   (Vl)" in line:
            sp = line.split()
            states[state].deg_vl.extend([int(x) for x in sp[3:]])
        elif "ROTATIONAL ENERGY LEVEL" in line:
            search = False


print("vec![")
for s in sorted(states):
    state = states[s]
    if all(x == 0 for x in state.non_deg + state.deg_vt + state.deg_vl):
        # ground state
        print(f"I1st(vec!{[0 for x in range(nvib)]}),")
    elif all(x == 0 for x in state.deg_vt + state.deg_vl):
        # nondeg states
        tmp = [0 for x in range(nvib)]
        if sum(state.non_deg) == 1:
            idx = state.non_deg.index(1)
            tmp[idx] = 1
            print(f"I1st(vec!{tmp}),")
        elif 2 in state.non_deg:
            idx = state.non_deg.index(2)
            tmp[idx] = 2
            print(f"I1st(vec!{tmp}),")
        else:
            for (i, v) in enumerate(state.non_deg):
                tmp[i] = v
            print(f"I1st(vec!{tmp}),")
    elif all(x == 0 for x in state.non_deg):
        # doubly-degenerate states
        tmp = [(0, 0) for x in range(nvib)]
        # only look at vt because they always go together I think
        if sum(state.deg_vt) == 1:
            idx = state.deg_vt.index(1)
            tmp[idx] = (1, 1)
            print(f"I2st(vec!{tmp}),")
        elif 2 in state.deg_vt:
            idx = state.deg_vt.index(2)
            tmp[idx] = (2, 0)
            print(f"I2st(vec!{tmp}),")
        else:
            for (i, v) in enumerate(state.deg_vt):
                tmp[i] = (v, v)
            print(f"I2st(vec!{tmp}),")
    else:
        # nondeg-deg combinations
        tmp1 = [0 for x in range(nvib)]
        i = state.non_deg.index(1)
        tmp1[i] = 1
        tmp2 = [(0, 0) for x in range(nvib)]
        i = state.deg_vt.index(1)
        tmp2[i] = (1, 1)
        print(
            f"""
        I12st {{
            i1st: vec!{tmp1},
            i2st: vec!{tmp2},
        }},
"""
        )
print("]")
