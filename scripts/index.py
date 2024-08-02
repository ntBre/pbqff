def fc2_index(n, a, b):
    sp = [a, b]
    sp.sort()
    [a, b] = sp
    return n * (a - 1) + b - 1


def fc3_index(i, j, k):
    sp = [i, j, k]
    sp.sort()
    return int(
        sp[0] + (sp[1] - 1) * sp[1] / 2 + (sp[2] - 1) * sp[2] * (sp[2] + 1) / 6 - 1
    )


def fc4_index(i, j, k, l):
    sp = [i, j, k, l]
    sp.sort()
    return int(
        sp[0]
        + (sp[1] - 1) * sp[1] / 2
        + (sp[2] - 1) * sp[2] * (sp[2] + 1) / 6
        + (sp[3] - 1) * sp[3] * (sp[3] + 1) * (sp[3] + 2) / 24
        - 1
    )


def index(n, a, b, c, d):
    match (c, d):
        case (0, 0):
            #return fc2_index(n, a, b)
            None
        case (_, 0):
            #return fc3_index(a, b, c)
            None
        case (_, _):
            return fc4_index(a, b, c, d)


ncoords = 15

for i in range(1, ncoords + 1):
    for j in range(1, i + 1):
        for k in range(0, j + 1):
            for l in range(0, k + 1):
                idx = index(ncoords, i, j, k, l)
                if idx is not None:
                    print(i, j, k, l, idx)
