import numpy as np


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
        return 0.0