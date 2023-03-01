import numpy as np


def sequence_alignment(x, y, match=3, mismatch=1, gap_open=1, gap_extend=0.5):
    nx = len(x)
    ny = len(y)
    # Optimal score at each possible pair of characters.
    F = np.zeros((nx + 1, ny + 1))
    F[:, 0] = np.linspace(0, -nx * gap_open, nx + 1)
    F[0, :] = np.linspace(0, -ny * gap_open, ny + 1)
    # Pointers to trace through an optimal aligment.
    P = np.zeros((nx + 1, ny + 1))
    P[:, 0] = 3
    P[0, :] = 4
    # Gap lengths for x and y sequences.
    F_gap_x = np.zeros((nx + 1, ny + 1))
    F_gap_x[:, 0] = np.inf
    F_gap_y = np.zeros((nx + 1, ny + 1))
    F_gap_y[0, :] = np.inf
    # Temporary scores.
    t = np.zeros(3)
    for i in range(nx):
        for j in range(ny):
            if x[i] == y[j]:
                t[0] = F[i, j] + match
            else:
                t[0] = F[i, j] - mismatch
            t[1] = F[i, j + 1] - gap_open - (F_gap_x[i, j] - 1) * gap_extend
            t[2] = F[i + 1, j] - gap_open - (F_gap_y[i, j] - 1) * gap_extend
            tmax = np.max(t)
            F[i + 1, j + 1] = tmax
            if t[0] == tmax:
                P[i + 1, j + 1] += 2
            if t[1] == tmax:
                P[i + 1, j + 1] += 3
                F_gap_x[i + 1, j + 1] = F_gap_x[i, j] + 1
            if t[2] == tmax:
                P[i + 1, j + 1] += 4
                F_gap_y[i + 1, j + 1] = F_gap_y[i, j] + 1
    # Trace through an optimal alignment.
    i = nx
    j = ny
    rx = []
    ry = []
    while i > 0 or j > 0:
        if P[i, j] in [2, 5, 6, 9]:
            rx.append(x[i - 1])
            ry.append(y[j - 1])
            i -= 1
            j -= 1
        elif P[i, j] in [3, 5, 7, 9]:
            rx.append(x[i - 1])
            ry.append('-')
            i -= 1
        elif P[i, j] in [4, 6, 7, 9]:
            rx.append('-')
            ry.append(y[j - 1])
            j -= 1
    # Reverse the strings.
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    return F[nx, ny], '\n'.join([rx, ry])


def read_sequences(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        x = lines[0].strip()
        y = lines[1].strip()
        return x, y


x, y = read_sequences('test3.seq')

score, alignment = sequence_alignment(x, y)
print(f"the alignment score is:", score)  # prints the optimal alignment score
print(alignment)  # prints the aligned strings
