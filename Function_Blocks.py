import numpy as np
import scipy.integrate as sci

class Calculate():
    def __init__(self, state_sv):
        self.state_sv = [state_sv]

    def ode(self, t, x, u):
        pass

    def calc(self, u, t0, t1):
        solv = sci.solve_ivp(self.ode, (t0, t1), self.state_sv, args=(u,))
        self.state_sv = solv.y.T[1]
        return solv.y.T[0][0]


class PT1_Glied(Calculate):
    def __init__(self, Kp, T, state_sv=0):
        super().__init__(state_sv)
        self.Kp = Kp
        self.T = T

    def ode(self, t, x, u):
        x1 = x
        dxdt = np.array([-x1/self.T+u*self.Kp/self.T])
        return dxdt


class I_Glied(Calculate):
    def __init__(self, Ki, state_sv=0):
        super().__init__(state_sv)
        self.Ki = Ki

    def ode(self, t, x, u):
        dxdt = np.array([self.Ki*u])
        return dxdt


class D_Glied:
    def __init__(self, Kd, u_init=None):
        self.Kd = Kd
        self.u_past = 0
        self.u_past_2 = 0

    def calc(self, u, h):
        if self.u_past is None:
            self.u_past = u

        #u_grad = [self.u_past, u]
        #du = np.gradient(u_grad)
        du = (u - self.u_past)/h
        self.u_past_2 = self.u_past
        self.u_past = u

        return du * self.Kd


class P_Glied:
    def __init__(self, Kp):
        self.Kp = Kp

    def calc(self, u):
        return self.Kp*u


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
        #return 0.5*(self.w1-self.w0)*(np.sign(t-self.t0)+1)+self.w0

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