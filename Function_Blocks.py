import numpy as np
import scipy.integrate as sci


class Calculate:
    def __init__(self, activate, state_sv):
        self.state_sv = [state_sv]
        self.activate = activate

    def ode(self, t, x, u):
        pass

    def calc(self, u, t0, t1):
        if self.activate:
            solv = sci.solve_ivp(self.ode, (t0, t1), self.state_sv, args=(u, ))
            self.state_sv = solv.y.T[1:, 0]
            return solv.y.T[0][0]
        else:
            return 0


class PT1_Glied(Calculate):
    def __init__(self, Kp, T, activate=True, state_sv=0):
        super().__init__(activate, state_sv)
        self.Kp = Kp
        self.T = T

    def ode(self, t, x, u):
        x1 = x
        dxdt = np.array([-x1/self.T+u*self.Kp/self.T])
        return dxdt


class I_Glied(Calculate):
    def __init__(self, Ki, activate=True, state_sv=0):
        super().__init__(activate, state_sv)
        self.Ki = Ki

    def ode(self, t, x, u):
        dxdt = np.array([self.Ki*u])
        return dxdt


class DT1_Glied():
    def __init__(self, Tp, Td, activate=True):
        self.Td = Td
        self.Tp = Tp
        self.pt1 = PT1_Glied(1, self.Tp)
        self.x1 = 0
        self.activate = activate

    def calc(self, u, t, t2):
        if self.activate:
            x1 = self.pt1.calc(u, t, t + t2)
            y = (-x1/self.Tp+u/self.Tp) * self.Td
            return y
        else:
            return 0


class D_Glied:
    def __init__(self, Kd, activate=True, u_init=None):
        self.Kd = Kd
        self.u_past = u_init
        self.activate = activate

    def calc(self, u, h):
        if self.activate:
            if self.u_past is None:
                self.u_past = u
            u_grad = [self.u_past, u]
            du = np.gradient(u_grad, h)
            self.u_past = u
            return du[0] * self.Kd
        else:
            return 0


class IT1_Glied(Calculate):
    def __init__(self, Kp, Ti, activate=True, state_sv=0):
        super().__init__(activate, state_sv)
        self.Kp = Kp
        self.Ti = Ti
        self.P = P_Glied(self.Kp)
        self.I = I_Glied(self.Kp / self.Ti)

    def calc(self, u, t, t2):
        y1 = self.P.calc(u)
        y2 = self.I.calc(u, t, t2)
        return y1+y2


class P_Glied:
    def __init__(self, Kp, activate=True):
        self.Kp = Kp
        self.activate = activate

    def calc(self, u):
        if self.activate:
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