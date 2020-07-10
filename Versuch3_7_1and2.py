import numpy as np
import scipy.integrate as sci
import matplotlib.pyplot as plt
import Function_Blocks as FB
import time

# zu plottende Aufgabe, Arbeitspunkt und Strecke auswählen
Aufgabennummer = "1a"



if Aufgabennummer == "1a":
    A1 = True
    genaeherteStrecke = False
    gekoppelt = False
    fuehrungssprung = True
elif Aufgabennummer == "1b":
    A1 = True
    genaeherteStrecke = False
    gekoppelt = True
    fuehrungssprung = True
elif Aufgabennummer == "1c":
    A1 = False
    genaeherteStrecke = False
    gekoppelt = True
    fuehrungssprung = True
elif Aufgabennummer == "1d":
    A1 = True
    genaeherteStrecke = False
    gekoppelt = False
    fuehrungssprung = False
elif Aufgabennummer == "1e":
    A1 = True
    genaeherteStrecke = False
    gekoppelt = True
    fuehrungssprung = False
elif Aufgabennummer == "2a":
    A1 = False
    genaeherteStrecke = True
    gekoppelt = True
    fuehrungssprung = True
else:
    A1 = False
    genaeherteStrecke = True
    gekoppelt = True
    fuehrungssprung = False


class Parameters(object):
    pass

# Simulationsparameter
sim_para = Parameters()  # instance of class Parameters
sim_para.t0 = 0        # start time
sim_para.tf = 20   # final time
sim_para.h = 0.001

sim_para.w0 = 0.0
sim_para.t_step = 1.0
sim_para.w1 = 1.0
sim_para.w2 = 0.0
sim_para.z11 = 1.0
sim_para.z21 = 0.0
sim_para.deltaT = 0.001

tt = np.arange(sim_para.t0, sim_para.tf + sim_para.h, sim_para.h)

#Erzeugen der Eingänge
w1 =FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.w1, sim_para.deltaT)
w2 =FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.w2, sim_para.deltaT)
z11 = FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.z11, sim_para.deltaT)
z21 = FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.z21, sim_para.deltaT)

# Listen in denen Ergebnisse für das spätere Plotten gespeichert werden

# Parameterfestlegung
if A1:
    Ki_p11 = 0.4
    Ki_p12 = -0.8
    Kp_p21 = -0.2
    Ti_p21 = 1
    Ki_p22 = 1.2
else:                   #A2
    Ki_p11 = 0.4
    Ki_p12 = -1.28
    Kp_p21 = -0.32
    Ti_p21 = 1
    Ki_p22 = 1.2

if genaeherteStrecke:
    Kp_r11 = 2.972
    Ti_r11 = 9.950
    Kp_r22 = 0.991
    Ti_r22 = 9.950
else:                   #resultierende Strecke
    Kp_r11 = 1.699
    Ti_r11 = 1.291
    Kp_r22 = 0.566
    Ti_r22 = 1.291

h = sim_para.h

# starte Timer zur Simulationszeitberechnung
time1 = time.perf_counter()
# Simulation


# weitere Aufgaben unter diese if Bedingungen
def ode(t, x):
    x3, x4, x5, x6, x7, x8 = x
    if fuehrungssprung:
        w1_cur = w1.linear(t)
        z11_cur = z11.zero()
    else:
        z11_cur = z11.linear(t)
        w1_cur = w1.zero()

    w2_cur = w2.zero()
    e_R11 = w1_cur - (x5+x7)
    e_R22 = w2_cur - (x6+x8)

    xdot3 = e_R11*Kp_r11/Ti_r11
    xdot4 = e_R22*Kp_r22/Ti_r22

    m1 = x3 + Kp_r11*e_R11
    m2 = x4 + Kp_r22*e_R22

    xdot5 = (m1+z11_cur) * Ki_p11
    if gekoppelt:
        xdot6 = -x6/Ti_p21+m1*Kp_p21/Ti_p21
        xdot7 = m2 * Ki_p12
    else:
        xdot6 = 0
        xdot7 = 0

    xdot8 = m2 * Ki_p22

    xdot=np.array([xdot3, xdot4, xdot5, xdot6, xdot7, xdot8])
    return xdot

x0 = np.array([0, 0, 0, 0, 0, 0])
solv = sci.solve_ivp(lambda t, x:  ode(t, x), (sim_para.t0, sim_para.tf), x0, t_eval=tt)
x3, x4, x5, x6, x7, x8 = solv.y

w1_traj = []
z11_traj = []

m1_traj = x3
m2_traj = x4

x1_traj = np.add(x5,x7)
x2_traj = np.add(x8,x6)
for x in tt:
    if fuehrungssprung:
        w1_traj.append(w1.linear(x))
        z11_traj.append(z11.zero())
    else:
        w1_traj.append(w1.zero())
        z11_traj.append(z11.linear(x))

# plot data
print("Simulationszeit: " + str(time.perf_counter()-time1))
fig1, (ax1_1, ax1_2) = plt.subplots(2)

ax1_1.plot(tt, w1_traj, color='blue', MarkerSize=1)
ax1_2.plot(tt, z11_traj)
# axis label
ax1_1.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Eingangssignal w1")
ax1_2.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Störsignal z11")
# format x axis
plt.tight_layout()
# adding grid
ax1_1.xaxis.grid(True)
ax1_1.yaxis.grid(True)
ax1_2.xaxis.grid(True)
ax1_2.yaxis.grid(True)
filename ="Grafiken/V7_1_" + Aufgabennummer + "_Eingaenge"
plt.savefig(filename, format="svg")
fig2, (ax2_1, ax2_2) = plt.subplots(2)

ax2_1.plot(tt, m1_traj, color='blue', MarkerSize=1)
ax2_2.plot(tt, m2_traj)
# axis label
ax2_1.set(xlabel="Zeit in s", title="m1")
ax2_2.set(xlabel="Zeit in s", title="m2")

# format x axis
plt.tight_layout()
# adding grid
ax2_1.xaxis.grid(True)
ax2_1.yaxis.grid(True)
ax2_2.xaxis.grid(True)
ax2_2.yaxis.grid(True)
filename ="Grafiken/V7_1_" + Aufgabennummer + "_m"
plt.savefig(filename, format="svg")
fig3, (ax3_1, ax3_2) = plt.subplots(2)

ax3_1.plot(tt, x1_traj)
ax3_2.plot(tt, x2_traj)
# axis label
ax3_1.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Ausgang x1")
ax3_2.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Ausgang x2")
# format x axis
plt.tight_layout()
# adding grid
ax3_1.xaxis.grid(True)
ax3_1.yaxis.grid(True)
ax3_2.xaxis.grid(True)
ax3_2.yaxis.grid(True)
filename ="Grafiken/V7_1_" + Aufgabennummer + "_Ausgaenge"
plt.savefig(filename, format="svg")
#plot anzeigen
plt.show()


