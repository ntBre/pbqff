import numpy as np


def print_mat(mat, rows, cols):
    for row in range(0, rows):
        for col in range(0, cols):
            print(f"{mat[row,col]:12.8f}", end="")
        print()


# fxm = np.loadtxt("../testfiles/ph3/fxm_full", dtype=np.double)
# want = np.loadtxt("../testfiles/ph3/pre_bdegnl_lxm", dtype=np.double)

fxm = np.loadtxt("../testfiles/c2h-/step_fxm", dtype=np.double)

(ROWS, COLS) = fxm.shape

w, v = np.linalg.eigh(fxm)

v = np.fliplr(v)

print("got=")
print_mat(v, ROWS, COLS)

# print("want=")
# print_mat(want, ROWS, COLS)

# print("diff=")
# print_mat(v - want, ROWS, COLS)

# EPS = 1e-6

# for col in range(0, COLS):
#     got = v[:, col]
#     iwant = want[:, col]
#     if (norm1 := np.linalg.norm(got - iwant)) > EPS:
#         got = -got
#         if (norm2 := np.linalg.norm(got - iwant)) > EPS:
#             print(f"diff of {min(norm1, norm2)} in col {col}")
#             assert False
