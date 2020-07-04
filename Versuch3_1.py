import numpy as np
from numpy import cos, sin, tan
import scipy.integrate as sci
import matplotlib.pyplot as plt
import Function_Blocks as FB
import time

# zu plottende Aufgabe, Arbeitspunkt und Strecke auswählen
Aufgabennummer = 1
A1 = True
normaleStrecke = True


class Parameters(object):
    pass

# Simulationsparameter
sim_para = Parameters()  # instance of class Parameters
sim_para.t0 = 0.0          # start time
sim_para.tf = 20.0     # final time
sim_para.h = 0.01

sim_para.w0 = 0.0
sim_para.t_step = 1.0
sim_para.w1 = 1.0
sim_para.w2 = 0.0
sim_para.z11 = 0.0
sim_para.z21 = 0.0
sim_para.deltaT = 0.001

tt = np.arange(sim_para.t0, sim_para.tf + sim_para.h, sim_para.h)

#Erzeugen der Eingänge
w1 =FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.w1, sim_para.deltaT)
w2 =FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.w2, sim_para.deltaT)
z11 = FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.z11, sim_para.deltaT)
z21 = FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.z21, sim_para.deltaT)

# Listen in denen Ergebnisse für das spätere Plotten gespeichert werden
m1_traj=  []
m2_traj = []
x1_traj = []
x2_traj = []
w1_traj = []
w2_traj = []



# Parameterfestlegung
if A1:
    Ti_p11 = 0.4
    Ti_p12 = -0.8
    Kp_p21 = -0.2
    Ti_p21 = 1
    Ti_p22 = 1.2
else:
    Ti_p11 = 0.4
    Ti_p12 = -1.28
    Kp_p21 = -0.32
    Ti_p21 = 1
    Ti_p22 = 1.2

if normaleStrecke:
    Kp_r11 = 1.699
    Ti_r11 = 1.291
    Kp_r22 = 0.566
    Ti_r22 = 1.291
else:
    Kp_r11 = 2.972
    Ti_r11 = 9.950
    Kp_r22 = 0.991
    Ti_r22 = 9.950

"""
Initialisierung der Blöcke
Übergabeparameter bestehen aus folgenden Blockparametern: Zeit-/Verstärkungskonstanten und activation Flag
activation Flag muss nicht übergeben werden (default ist True)
um Block zu berechnen dessen .calc() Methode aufrufen: 
calc_Parameter: Eingabe vom Block, aktueller Zeitpunkt, aktueller Zeitpunkt+1 (--> aktueller Zeitpunkt + Schrittweite)
"""
R11_block = FB.IT1_Glied(Kp_r11, Ti_r11)
R22_block = FB.IT1_Glied(Kp_r22, Ti_r22)
P11_block = FB.I_Glied(Ti_p11)
P21_block = FB.PT1_Glied(Kp_p21, Ti_p21)
P12_block = FB.I_Glied(Ti_p12)
P22_block = FB.I_Glied(Ti_p22)
calc = FB.Simple_Calc()

x1 = 0.0
x2 = 0.0
h = sim_para.h

# starte Timer zur Simulationszeitberechnung
time1 = time.perf_counter()
# Simulation
if Aufgabennummer == 1:
    for x in range(len(tt)):
        t = tt[x]
        w1_cur = w1.linear(t)
        w2_cur = w2.zero()
        w1_traj.append(w1_cur)
        w2_traj.append(w2_cur)
        z11_cur = z11.zero()
        z21_cur = z21.zero()

        e1 = calc.sub(w1_cur, x1)
        e2 = calc.sub(w2_cur, x2)
        m1 = R11_block.calc(e1, t, t+h)
        m1_traj.append(m1)

        m2 = R22_block.calc(e2, t, t+h)
        m2_traj.append(m2)

        u11 = calc.add(z11_cur, m1)
        u21 = calc.add(z21_cur, m1)

        y11 = P11_block.calc(u11, t, t+h)
        y21 = P21_block.calc(u21, t, t+h)
        u12 = m2
        u22 = m2
        y12 = P12_block.calc(u12, t, t+h)
        y22 = P22_block.calc(u22, t, t+h)

        x1 = calc.add(y11, y12)
        x2 = calc.add(y22, y21)
        x1_traj.append(x1)
        x2_traj.append(x2)

# weitere Aufgaben unter diese if Bedingungen
elif Aufgabennummer == 2:
    pass
elif Aufgabennummer == 3:
    pass


# plot data
print("Simulationszeit: " + str(time.perf_counter()-time1))
fig1, (ax1_1, ax1_2) = plt.subplots(2)
ax1_1.plot(tt, w1_traj, '-o', color='blue', MarkerSize=2)
ax1_2.plot(tt, w2_traj)
# axis label
ax1_1.set(xlabel="Zeit", title="Eingangssignal w1")
ax1_2.set(xlabel="Zeit", title="Eingangssignal w2")
# format x axis
plt.tight_layout()
# adding grid
ax1_1.xaxis.grid(True)
ax1_1.yaxis.grid(True)
ax1_2.xaxis.grid(True)
ax1_2.yaxis.grid(True)

fig2, (ax2_1, ax2_2) = plt.subplots(2)
ax2_1.plot(tt, m1_traj, '-o', color='blue', MarkerSize=2)
ax2_2.plot(tt, m2_traj)
# axis label
ax2_1.set(xlabel="Zeit", title="m1")
ax2_2.set(xlabel="Zeit", title="m2")

# format x axis
plt.tight_layout()
# adding grid
ax2_1.xaxis.grid(True)
ax2_1.yaxis.grid(True)
ax2_2.xaxis.grid(True)
ax2_2.yaxis.grid(True)

fig1, (ax3_1, ax3_2) = plt.subplots(2)
ax3_1.plot(tt, x1_traj)
ax3_2.plot(tt, x2_traj)
# axis label
ax3_1.set(xlabel="Zeit", title="Ausgang x1")
ax3_2.set(xlabel="Zeit", title="Ausgang x2")
# format x axis
plt.tight_layout()
# adding grid
ax3_1.xaxis.grid(True)
ax3_1.yaxis.grid(True)
ax3_2.xaxis.grid(True)
ax3_2.yaxis.grid(True)
#plot anzeigen
plt.show()



