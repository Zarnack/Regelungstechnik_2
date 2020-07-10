import numpy as np
import scipy.integrate as sci
import matplotlib.pyplot as plt
import Function_Blocks as FB
import time

# zu plottende Aufgabe, Arbeitspunkt und Strecke auswählen

A1 = False
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
z11_traj = []
z21_traj = []



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

zp = []
m1_traj = []
m2_traj = []
w1_traj = []
w2_traj = []
# weitere Aufgaben unter diese if Bedingungen
def ode(t, x):

    x1, x2, x3, x4, e1_dt, e2_dt, x7 = x

    r1 = x1 + x2
    e1 = w1.jump(t) - r1
    x_dot5 = e1     #Hilfsgröße für m1 -> Berechung |e2*dt

    r2 = x3 + x4
    e2 = w2.zero() - r2
    x_dot6 = e2

    m1 = Kp_r11*e1 + Kp_r11/Ti_r11*e1_dt    #Lösung PI-Regler im Zeitbereich
    m2 = Kp_r22*e2 + Kp_r22/Ti_r22*e2_dt

    # Eingang DT1-Glied ohne Abhängigkeit von sich selbst
    a1 = 1/(1-(Ki_p12 * Kp_p21)/(Ki_p11 * Ki_p22*T_p21)) * (m1-Ki_p12/Ki_p11 * m2 - (Ki_p12 * Kp_p21)/(Ki_p11 * Ki_p22 * T_p21) * x7)

    x_dot7 = 1/T_p21*(a1-x7)        #DGL PT1-Glied
    y = x_dot7*(-Kp_p21/Ki_p22)     #Lösung DT1-Glied

    a2 = m2 + y

    x_dot1 = Ki_p11*(a1+z11.zero())                    #DGL I-Glied
    x_dot2 = Ki_p12*(a2)                               #DGL I-Glied
    x_dot3 = 1/T_p21*(-x3+Kp_p21*(a1+z21.zero()))      #DGL PT1-Glied
    x_dot4 = Ki_p22*(a2)                               #DGL I-Glied


    x_dot = np.array([x_dot1, x_dot2, x_dot3, x_dot4, x_dot5, x_dot6, x_dot7])
    return x_dot


x0 = np.array([0,0,0,0,0,0,0])
solv = sci.solve_ivp(lambda t, x:  ode(t, x), (sim_para.t0, sim_para.tf), x0, t_eval=tt)
x1_traj = solv.y.T[:,0] + solv.y.T[:,1] #2 param wahl x
x2_traj = solv.y.T[:,2] + solv.y.T[:,3]

for tx in tt:
    w1_traj.append(w1.linear(tx))
    w2_traj.append(w2.zero())
    z11_traj.append(z11.zero())
    z21_traj.append(z21.zero())


for tx in range(len(tt)):
    m1 = Kp_r11*(w1_traj[tx] - solv.y.T[tx,0] - solv.y.T[tx,1]) + Kp_r11/Ti_r11*solv.y.T[tx,4]    #Lösung PI-Regler im Zeitbereich
    m1_traj.append(m1) 
    m2 = Kp_r22*(w1_traj[tx] - solv.y.T[tx,2] - solv.y.T[tx,3]) + Kp_r22/Ti_r22*solv.y.T[tx,5]
    m2_traj.append(m2)

    

print("Simulationszeit: " + str(time.perf_counter()-time1))
fig1, (ax1_1, ax1_2, axz11, axz21) = plt.subplots(4)

ax1_1.plot(tt, w1_traj) # '-o', color='blue', MarkerSize=2
ax1_2.plot(tt, w2_traj)
# axis label
ax1_1.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Eingangssignal w1")
ax1_2.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Eingangssignal w2")
plt.tight_layout()
# format x axis
# adding grid
ax1_1.xaxis.grid(True)
ax1_1.yaxis.grid(True)
ax1_2.xaxis.grid(True)
ax1_2.yaxis.grid(True)


axz11.plot(tt, z11_traj) 
axz21.plot(tt, z21_traj)


# axis label
axz11.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Störsignal z11")
axz21.set(xlabel="Zeit in s", ylabel="Spannung in V", title="Störsignal z21")
plt.tight_layout()
# format x axis
# adding grid
axz11.xaxis.grid(True)
axz11.yaxis.grid(True)
axz21.xaxis.grid(True)
axz21.yaxis.grid(True)


fig2, (ax2_1, ax2_2, ax3_1, ax3_2) = plt.subplots(4)
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

plt.show()