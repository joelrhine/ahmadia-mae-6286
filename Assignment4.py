import numpy
from matplotlib import pyplot
import matplotlib.cm as cm
%matplotlib inline


n = 192
Du, Dv, F, k = 0.00016, 0.00008, 0.035, 0.065 # Bacteria 1
dh = 5/(n-1)
T = 8000
dt = .9 * dh**2 / (4*max(Du,Dv))
nt = int(T/dt)


uvinitial = numpy.load('c:/users/joel/desktop/uvinitial.npz')
U_initial = uvinitial['U']
V_initial = uvinitial['V']


U_initial


def assignment4(nt,U_initial,V_initial,Du, Dv, F, k):
    U = U_initial.copy()
    V = V_initial.copy()

    for t in range(1,nt):
        Un = U.copy()
        Vn = V.copy()
        U[1:-1,1:-1] = Un[1:-1,1:-1]+\
                   Du*dt/(dh**2)*((Un[1:-1,2:]-2*Un[1:-1,1:-1]+Un[1:-1,:-2])+\
                                  (Un[2:,1:-1]-2*Un[1:-1,1:-1]+Un[:-2,1:-1]))-\
                  dt*Un[1:-1,1:-1]*Vn[1:-1,1:-1]**2+\
                  F*(1-Un[1:-1,1:-1])*dt
        V[1:-1,1:-1] = Vn[1:-1,1:-1]+\
                   Dv*dt/(dh**2)*(Vn[1:-1,2:]-2*Vn[1:-1,1:-1]+Vn[1:-1,:-2]+Vn[2:,1:-1]-2*Vn[1:-1,1:-1]+Vn[:-2,1:-1])+\
                   dt*Un[1:-1,1:-1]*Vn[1:-1,1:-1]**2-\
                   (F+k)*Vn[1:-1,1:-1]*dt
    #BC
        U[0,1:-1] = U[1,1:-1]
        U[-1,1:-1] = U[-2,1:-1]
        U[1:-1,0] = U[1:-1,1]
        U[1:-1,-1] = U[1:-1,-2]
        U[0,0] = U[1,1]
        U[0,-1] = U[1,-2]
        U[-1,0] = U[-2,1]
        U[-1,-1] = U[-2,-2]

    #corner
        V[0,1:-1] = V[1,1:-1]
        V[-1,1:-1] = V[-2,1:-1]
        V[1:-1,0] = V[1:-1,1]
        V[1:-1,-1] = V[1:-1,-2]
        V[0,0] = V[1,1]
        V[0,-1] = V[1,-2]
        V[-1,0] = V[-2,1]
        V[-1,-1] = V[-2,-2]

    a = [U, V]
    return a


# k = 0.00016, 0.00008, 0.035, 0.065 # Bacteria 1 
U1 = assignment4(nt,U_initial, V_initial,0.00016, 0.00008, 0.035, 0.065 )[0]
V1 = assignment4(nt,U_initial, V_initial,0.00016, 0.00008, 0.035, 0.065 )[1]


U1[100,::40]
