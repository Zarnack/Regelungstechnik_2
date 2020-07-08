import numpy as np
import scipy.integrate as sci


class PT1_Glied():
    def __init__(self, Kp, T, activate=True):
        self.Kp = Kp
        self.T = T

    def calc(self, u, x):
        dxdt = -x/self.T+u*self.Kp/self.T
        return dxdt


class I_Glied():
    def __init__(self, Ki, activate=True):
        self.Ki = Ki

    def calc(self, u):
        dxdt = self.Ki*u
        return dxdt


class DT1_Glied():
    def __init__(self, Tp, Td, activate=True):
        self.Td = Td
        self.Tp = Tp
        self.pt1 = PT1_Glied(1, self.Tp)
        self.x1 = 0
        self.activate = activate

    def calc(self, u, x):
        if self.activate:
            x1 = self.pt1.calc(u, x)
            y = (-x1/self.Tp+u/self.Tp) * self.Td
            return y
        else:
            return 0


class D_Glied:
    def __init__(self, Kd, t_init, activate=True, u_init=None):
        self.Kd = Kd
        self.u_past = u_init
        self.t_past = t_init
        self.activate = activate

    def calc(self, u, t):
        if self.activate:
            if self.u_past is None:
                self.u_past = u
            u_grad = [self.u_past, u]
            try:
                du = np.gradient(u_grad, t-self.t_past)
            except ValueError:
                du = [0]
                self.u_past = u
                self.t_past = t

            return du[0] * self.Kd
        else:
            return 0


class IT1_Glied():
    def __init__(self, Kp, Ti, activate=True, state_sv=0):
        self.Kp = Kp
        self.Ti = Ti
        self.P = P_Glied(self.Kp)
        self.I = I_Glied(self.Kp / self.Ti)

    def calc(self, u):
        y1 = self.P.calc(u)
        y2 = self.I.calc(u)
        return y1


class P_Glied:
    def __init__(self, Kp, activate=True):
        self.Kp = Kp
        self.activate = activate

    def calc(self, u):
        if self.activate == true:  
            return self.Kp*u
        else:
            return 0


class Simple_Calc:
    def __init__(self):
        pass

    def sub(self, u1, u2):
        return u1 - u2

    def add(self, u1, u2):
        return u1 + u2

    def multiply(self, u1, Kp):
        return u1*Kp


class Inputfunction:
    def __init__(self, t0=0.0, w0=0.0, w1=1.0, deltaT=0.1):
        self.t0 =t0
        self.w1 = w1
        self.w0 =w0
        self.deltaT = deltaT

    def jump(self, t):
        return self.w0 if t < self.t0 else self.w1

    def jump_v(self, t):
        u1 = np.ones(t.size)*self.w0
        u1[t >= self.t0] = self.w1
        return u1

    def linear(self, t):
        if t < self.t0:
            return self.w0
        if self.t0 <= t <= self.t0+self.deltaT:
            return 1/2*(self.w1-self.w0)*np.tanh(2*np.pi/self.deltaT*(t-self.t0))+1/2*(self.w1+self.w0)
        if t > self.t0+self.deltaT:
            return self.w1

    def zero(self):
        return 0