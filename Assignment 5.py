import numpy as np

nx = 41
ny = 41

l = 1.0
h = 1.0


dx = l / (nx - 1)
dy = h / (ny - 1)

u = 1  # (m/s)


l1_target = 1e-6




def l_norm(new, old):
    norm = np.sum(np.abs(new - old))

    return norm


def Lap2D(w, p, dx, dy, l1_target):

    l1_norm = 1
    iteration = 0

    wn = np.empty_like(w)
    pn = np.empty_like(p)

    while l1_norm > l1_target:
        wn = w.copy()
        pn = p.copy()

        w[1: -1, 1: -1] = (0.25 * (wn[1:-1, 2:] + wn[1:-1, :-2] +
                                   wn[2:, 1:-1] + wn[:-2, 1:-1]))


        w[-1, :] = (- 0.5 / dy**2 * (8 * p[-2, :] - p[-3, :]) -
                    3 * 1 / dy)
        w[0, :] = - 0.5 / dy**2 * (8 * p[1, :] - p[2, :])
        w[:, -1] = - 0.5 / dx**2 * (8 * p[:, -2] - p[:, -3])
        w[:, 0] = - 0.5 / dx**2 * (8 * p[:, 1] - p[:, 2])

        p[1:-1, 1:-1] = (1 / (2 * (dx**2 + dy**2)) *
                         ((pn[1:-1, 2:] + pn[1:-1, :-2]) * dy**2 +
                          (pn[2:, 1:-1] + pn[:-2, 1:-1]) * dx**2 +
                          w[1:-1, 1:-1] * dx**2 * dy**2))


        l1_norm = max(l_norm(w, wn), l_norm(p, pn))

    return w, p


if __name__ == '__main__':
    x = np.linspace(0, l, nx)
    y = np.linspace(0, h, ny)

    wi = np.zeros((ny, nx))
    pi = np.zeros((ny, nx))

    w, p = Lap2D(wi.copy(), pi.copy(), dx, dy, l1_target)

    print(np.round(np.max(np.abs(p)), 4))
    print(np.round(np.max(np.abs(w)), 4))
    print(np.round(p[32, ::8], 4))
