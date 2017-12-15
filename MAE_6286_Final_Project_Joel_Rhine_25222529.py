
import numpy as np
import matplotlib.pyplot as plt

L = 0.1 # (m), Thickness of the wall
n = 20 # (Constant), Number of equal sections
alpha = 0.0001 # (Constant), Thermal Diffusivity
T0 = 0 # (C), Temperature of sections between the two surfaces
T1s = 40 # (C), Temperature of Surface 1
T2s = 20 # (C), Temperature of Surface 2
dx = L / n # (Constant), Width of each Section
t_final = 60 # (s), Time cycle of 60 secs
dt = 0.1 # (Constant), Time Step

x = np.linspace(dx / 2, L - dx / 2, n) # Defining the positions of nodes at the center of each section

T_2 = np.ones(n) * T0 # Initializing Temperatures T_2 to 0
T_4 = np.ones(n) * T0 # Initializing Temperatures T_4 to 0

dTdt_2 = np.empty(n)
dTdt_4 = np.empty(n)

t = np.arange(0, t_final, dt) # Initializing Time Cycle

# Calculation of Temperature values over Time

for j in range(1, len(t)):
    plt.clf()

    # Second Order Accuracy

    dTdt_2[0] = alpha * ((T1s - T_2[0]) / dx ** 2 + (T_2[1] - T_2[0]) / dx ** 2)
    dTdt_2[n - 1] = alpha * ((T_2[n - 2] - T_2[n - 1]) / dx ** 2 + (T2s - T_2[n - 1]) / dx ** 2)

    for i in range(1, n - 1):
        dTdt_2[i] = alpha * (-(T_2[i] - T_2[i - 1]) / dx ** 2 + (T_2[i + 1] - T_2[i]) / dx ** 2)

    # Fourth Order Accuracy

    dTdt_4[0] = dTdt_2[0]
    dTdt_4[n-1] = dTdt_2[n-1]

    dTdt_4[1] = dTdt_2[1]
    dTdt_4[n-2] = dTdt_2[n-2]
    for i in range(2, n-2):
        dTdt_4[i] = alpha*((-T_2[i+2] + 16*T_2[i+1] - 30*T_2[i] + 16*T_2[i-1] - T_2[i-2])/(12*dx**2))


    T_2 += dTdt_2 * dt

    T_4 += dTdt_4 * dt

# Plotting of Values

    plt.grid()
    plt.plot(x, T_2, color = 'r', marker = 'o', ls = ' ')
    plt.plot(x, T_4, color = 'b', ls = '--')
    plt.legend(['2nd Order Accuracy','4th Order Accuracy'])
    plt.xlabel('Distance in m')
    plt.ylabel('Temperature in Clesius')
    plt.axis([0, L, 0, 50])
    plt.pause(0.05)


# Difference in Accuracy of Equilibrium Results

T = np.empty(n)
T = (T_2 - T_4)
print(T)
