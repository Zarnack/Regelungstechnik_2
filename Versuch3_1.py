import numpy as np
from numpy import cos, sin, tan
import scipy.integrate as sci
import matplotlib.pyplot as plt
import Function_Blocks as FB

# LISTING_START ParaClass

class Parameters(object):
    pass

# LISTING_START SimuPara
# Simulation parameter
sim_para = Parameters()  # instance of class Parameters
sim_para.t0 = 0.0          # start time
sim_para.tf = 10.0     # final time
sim_para.h = 0.001

sim_para.w0 = 0.0
sim_para.t_step = 0.0
sim_para.w1 = 2.0
sim_para.deltaT = 0.1


#step = Inputfunctions(sim_para.t_step, sim_para.w0, sim_para.w1, sim_para.t_step+sim_para.delta_T)
input_step =FB.Inputfunction(sim_para.t_step, sim_para.w0, sim_para.w1, sim_para.deltaT)
tt = np.arange(sim_para.t0, sim_para.tf + sim_para.h, sim_para.h)
#u = step.jump_v(tt)

pt1_block = FB.PT1_Glied(1, 1)
I_block = FB.I_Glied(2)
sub = FB.Simple_Calc()
d_block = FB.D_Glied(1)

# simulation
x_traj = []
y_out= [0]
u_traj = []
for x in range(len(tt)):
    t = tt[x]
    u_cur = input_step.linear(t)
    u_traj.append(u_cur)

    error = sub.sub(u_cur, y_out[0])


    y_pt1 = pt1_block.calc(error, t, t + sim_para.h)
    #y_I = I.calc(y_pt1, t, t+sim_para.h)
    y_out = d_block.calc(y_pt1[0])
    x_traj.append(y_out[0])


fig1, (ax1, ax2) = plt.subplots(2)
ax1.plot(tt, x_traj, '-o', color='blue', MarkerSize=2)
ax2.plot(tt, u_traj)
# axis label
ax1.set(xlabel="Zeit", title="Sprungantwort")
ax2.set(xlabel="Zeit", title="Eingangssignal")
ax1.set_ylabel("", color='blue')
# format x axis
plt.tight_layout()
# adding grid
ax1.xaxis.grid(True)
ax1.yaxis.grid(True)
plt.show()
