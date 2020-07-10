import numpy as np
import scipy.integrate as sci
import matplotlib.pyplot as plt
import Function_Blocks as FB
import time

# zu plottende Aufgabe, Arbeitspunkt und Strecke auswählen
A1 = True
normaleStrecke = True
gekoppelt = False


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
sim_para.z12 = 1.0
sim_para.z21 = 1.0
sim_para.deltaT = 0.001

tt = np.arange(sim_para.t0, sim_para.tf + sim_para.h, sim_para.h)

#Erzeugen der Eingänge
w1 =FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.w1, sim_para.deltaT)
w2 =FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.w2, sim_para.deltaT)
z12 = FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.z12, sim_para.deltaT)
z21 = FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.z21, sim_para.deltaT)


# Listen in denen Ergebnisse für das spätere Plotten gespeichert werden
m1_traj=  []
m2_traj = []
x1_traj = []
x2_traj = []
w1_traj = []
w2_traj = []
z12_traj = []
z21_traj = []
zp = []


# Parameterfestlegung
if A1:
    Ki_p11 = 0.4    #k1
    Ki_p12 = -0.8   #k3
    Kp_p21 = -0.2   #k4
    T_p21 = 1       #T
    Ki_p22 = 1.2    #k2
else:                   #A2
    Ki_p11 = 0.4
    Ki_p12 = -1.28
    Kp_p21 = -0.32
    T_p21 = 1
    Ki_p22 = 1.2

if normaleStrecke:
    Kp_r11 = 1.699
    Ti_r11 = 1.291
    Kp_r22 = 0.566
    Ti_r22 = 1.291
else:                   #aus resultierende Strecke
    Kp_r11 = 2.972
    Ti_r11 = 9.950
    Kp_r22 = 0.991
    Ti_r22 = 9.950


calcus = FB.Simple_Calc()


test = []
timer = []

h = sim_para.h

# starte Timer zur Simulationszeitberechnung
time1 = time.perf_counter()
# Simulation


# weitere Aufgaben unter diese if Bedingungen
def ode(t, x):
    x1, x2, x3, x4, x5, x6 = x

    w1_cur = w1.zero()
    w2_cur = w2.zero()

    e_R11 = w1_cur + (x3+x4)
    e_R22 = w2_cur + (x1+x2)

    xdot5 = e_R11
    xdot6 = e_R22

    m1 = x5*Kp_r11/Ti_r11 + Kp_r11*e_R11
    m2 = x6*Kp_r22/Ti_r22 + Kp_r22*e_R22

    xdot2 = (m2+z12.linear(t))*Ki_p12
    xdot3 =  1/T_p21*(-x3+Kp_p21*(m1+z21.linear(t)))

    if gekoppelt:
       xdot1 = m1*Ki_p11
       xdot4 = m2*Ki_p22
    else:
        xdot1 = 0
        xdot4 = 0

    xdot=np.array([xdot1, xdot2, xdot3, xdot4, xdot5, xdot6])
    return xdot
   
 
x0 = np.array([0, 0, 0, 0, 0, 0])
solv = sci.solve_ivp(lambda t, x:  ode(t, x), (sim_para.t0, sim_para.tf), x0, t_eval=tt)
x1_traj = solv.y.T[:,0] + solv.y.T[:,1] #2 param wahl x
x2_traj = solv.y.T[:,2] + solv.y.T[:,3]

for tx in tt:
    w1_traj.append(w1.zero())
    w2_traj.append(w2.zero())
    z12_traj.append(z12.linear(tx))
    z21_traj.append(z21.linear(tx))


for tx in range(len(tt)):
    m1 = solv.y.T[tx,4] + Kp_r11*(w1_traj[tx] + (solv.y.T[tx,2]+solv.y.T[tx,3]))
    m1_traj.append(m1) 
    m2 = solv.y.T[tx,5] + Kp_r22*(w2_traj[tx] + (solv.y.T[tx,0]+solv.y.T[tx,1]))
    m2_traj.append(m2)
    

print("Simulationszeit: " + str(time.perf_counter()-time1))
fig1, (ax1_1, ax1_2, axz12, axz21) = plt.subplots(4)

ax1_1.plot(tt, w1_traj) # '-o', color='blue', MarkerSize=2
ax1_2.plot(tt, w2_traj)
# axis label
ax1_1.set(xlabel="Zeit", title="Eingangssignal w1")
ax1_2.set(xlabel="Zeit", title="Eingangssignal w2")
plt.tight_layout()
# format x axis
# adding grid
ax1_1.xaxis.grid(True)
ax1_1.yaxis.grid(True)
ax1_2.xaxis.grid(True)
ax1_2.yaxis.grid(True)


axz12.plot(tt, z12_traj) # '-o', color='blue', MarkerSize=2
axz21.plot(tt, z21_traj)


# axis label
axz12.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Störsignal z12")
axz21.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Störsignal z21")
plt.tight_layout()
# format x axis
# adding grid
axz12.xaxis.grid(True)
axz12.yaxis.grid(True)
axz21.xaxis.grid(True)
axz21.yaxis.grid(True)

filename ="Grafiken/V7_4_1_Eingaenge"
plt.savefig(filename, format="svg")
fig2, (ax2_1, ax2_2) = plt.subplots(2)
plt.tight_layout()
ax2_1.plot(tt, m1_traj)
ax2_2.plot(tt, m2_traj)

# axis label
ax2_1.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Signal m1")
ax2_2.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Signal m2")
plt.tight_layout()
# format x axis
# adding grid
ax2_1.xaxis.grid(True)
ax2_1.yaxis.grid(True)
ax2_2.xaxis.grid(True)
ax2_2.yaxis.grid(True)

filename ="Grafiken/V7_4_1_m"
plt.savefig(filename, format="svg")
fig3, (ax3_1, ax3_2) = plt.subplots(2)
plt.tight_layout()
ax3_1.plot(tt, x1_traj)
ax3_2.plot(tt, x2_traj)

# axis label
ax3_1.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Ausgangssignal x1")
ax3_2.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Ausgangssignal x2")
plt.tight_layout()
# format x axis

# adding grid
ax3_1.xaxis.grid(True)
ax3_1.yaxis.grid(True)
ax3_2.xaxis.grid(True)
ax3_2.yaxis.grid(True)
filename ="Grafiken/V7_4_1_Ausgaenge"
plt.savefig(filename, format="svg")
plt.show()