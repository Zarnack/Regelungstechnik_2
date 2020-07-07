import numpy as np
from numpy import cos, sin, tan
import scipy.integrate as sci
import matplotlib.pyplot as plt
import Function_Blocks as FB
import time

# zu plottende Aufgabe, Arbeitspunkt und Strecke auswählen
Aufgabennummer = 4
A1 = True               
normaleStrecke = True   

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
else:                   #A2
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
else:                   #resultierende Strecke
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
# R11_block = FB.IT1_Glied(Kp_r11, Ti_r11)
# R22_block = FB.IT1_Glied(Kp_r22, Ti_r22)
# P11_block = FB.I_Glied(Ti_p11)
# P21_block = FB.PT1_Glied(Kp_p21, Ti_p21)
# P12_block = FB.I_Glied(Ti_p12)
# P22_block = FB.I_Glied(Ti_p22)
# V11_block = FB.P_Glied(1)
# V22_block = FB.P_Glied(1)
# V12_block = FB.P_Glied(3.2)
# V21_block = FB.DT1_Glied(0.27, 1.0)
calcus = FB.Simple_Calc()


test = []
timer = []

h = sim_para.h

# starte Timer zur Simulationszeitberechnung
time1 = time.perf_counter()
# Simulation


# weitere Aufgaben unter diese if Bedingungen
def ode(t, x):
    if Aufgabennummer == 1:
        pass

    elif Aufgabennummer == 2:
        pass

    elif Aufgabennummer == 3:
        pass

    elif Aufgabennummer == 4:
        pass


x0 = np.array([0,0,0,0,0,0,0])
solv = sci.solve_ivp(lambda t, x:  ode(t, x), (sim_para.t0, sim_para.tf), x0, min_step=1, t_eval=tt)
x1 = solv.y.T[:,0]



def plot():
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


print("Simulationszeit: " + str(time.perf_counter()-time1))
fig1, (ax1_1) = plt.subplots(1)
ax1_1.plot(tt, x1, '-o', color='blue', MarkerSize=2)
# axis label
ax1_1.set(xlabel="Zeit", title="Eingangssignal w1")
# format x axis
plt.tight_layout()
# adding grid
ax1_1.xaxis.grid(True)
ax1_1.yaxis.grid(True)
plt.show()
