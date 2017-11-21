import numpy


nx = 81
dx = 0.25
dt = 0.0002
gamma = 1.4


x = numpy.linspace(-10, 10, nx)
rho_IC = numpy.empty(nx)
u_IC = numpy.empty(nx)
p_IC = numpy.empty(nx)

for i in range(len(x)):
    if x[i]<0:
        rho_IC[i] = 1.0
        u_IC[i] = 0.0
        p_IC[i] = 100000
    else:
        rho_IC[i] = 0.125
        u_IC[i] = 0.0
        p_IC[i] = 10000.0

rho = rho_IC
u = u_IC
p = p_IC


T = 0.01
nt = int(T/dt)+1


def compute_eT(rho, u, p, gamma):

    e = p / (gamma - 1) / rho

    return e + 0.5 * u ** 2

def get_u(rho, u, p):

    eT = compute_eT(rho, u, p, gamma)

    return numpy.array([rho, rho * u, rho * eT])


def f(U):
    rho = U[0, :]
    p = U[-1, :]
    eT = compute_eT(rho, u, p, gamma)

    return numpy.array([U[1], \
                        U[1] ** 2 / U[0] + (gamma - 1) * (U[-1] - 0.5 * U[1] ** 2 / U[0]), \
                        (U[-1] + (gamma - 1) * (U[-1] - 0.5 * U[1] ** 2 / U[0])) * U[1] / U[0]])


u_1 = get_u(rho, u, p)
u_nn = numpy.empty_like(u_1)
u_result = numpy.empty_like(u_1)
F = numpy.empty_like(u_1)
F_nn = numpy.empty_like(u_1)

for t in range(1, nt):
    F = f(u_1)

    u_nn[:, :-1] = 0.5 * (u_1[:, 1:] + u_1[:, :-1]) \
                   - (dt / 2 / dx) * (F[:, 1:] - F[:, :-1])

    u_nn[:, -1] = u_1[:, -1]

    F_nn = f(u_nn)

    u_result[:, 1:-1] = u_1[:, 1:-1] - dt / dx * (F_nn[:, 1:-1] - F_nn[:, :-2])

    u_result[:, 0] = u_1[:, 0]
    u_result[:, -1] = u_1[:, -1]
    u_1 = u_result.copy()

for i in range(nt):
    if x[i] == 2.5:
        print(i)
        print(u_result[1, i]/u_result[0, i])
        print(u_result[0, i])
        print(u_result[-1,  i])
        print((gamma - 1)*(u_result[-1,i]-0.5*u_result[1, i]**2/u_result[0, i]))
