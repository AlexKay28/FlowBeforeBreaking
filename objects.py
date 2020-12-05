import math
from scipy.integrate import quad
import numpy as np
import numpy as np
import sympy

from input_params import *

class Deffect:
    """
    Состояние трещины
    """
    def __init__(self, a, c):
        self.a = a
        self.c = c

    @property
    def double_c(self):
        return 2 * self.c

class Steel:
    """
    Параметры стали
    """
    def __init__(self, C, m, Rp02_min, Rp02_max, Rm_min, Rm_max, E_module, mu):
        self.C  = C
        slef.m  = m
        self.mu = mu
        self.Rp02_min = Rp02_min # Pa
        self.Rp02_max = Rp02_max # Pa
        self.Rm_min   = Rm_min   # Pa
        self.Rm_max   = Rm_max   # Pa
        self.E_module = E_module # Pa

class Problem:
    """
    Состояние системы в целом. Параметры окружающей среды и трещина
    """
    def __init__(self, deffect, t, Dout, p, T, M, C, m):
        self.deffect = deffect
        self.t = t
        self.Dout = Dout
        self.Din  = Dout - 2*t
        self.p = p
        self.T = T
        self.M = M
        self.C = C
        self.m = m

    @property
    def Rout(self):
        return self.Dout/2

    @property
    def Rin(self):
        return self.Din/2

    @property
    def ro(self):
        if self.T == 20:
            return 7856
        elif self.T == 150:
            return 7819
        elif self.T == 285:
            return 7775
        raise ValueError('Cant calculate the density.')

class Solve:
    """
    Объект решения
    """
    def __init__(self, problem, deffect, type_):
        self.problem = problem
        self.deffect = deffect
        self.type = type_

    def gamma_c(self):
        return 1

    def gamma_d(self):
        return (1.1 + 0.35 * (self.deffect.a / self.problem.t)**2) * (self.deffect.a/self.deffect.c)**0.5

    def m_c(self):
        a = self.deffect.a
        c = self.deffect.c
        t = self.problem.t
        Rin = self.problem.Rin
        ro = self.problem.ro
        if self.type == 'Кольцевой эффект':
            m = 1 - a / (t * Rin)**0.5 * (1 - (a / c)**0.5)
        elif self.type == 'Продольный эффект':
            numerator = 1 + 4 * (a / ro) * (1 - (a/c)**0.5)
            denomerator = 1 + 5 * a/Rin * (1 - (a/c)**0.5) * (1 + 2 * (a/t)**2)
            m = numerator / denomerator
        return m

    def m_d(self):
        if self.type == 'Кольцевой эффект':
            m = 1
        elif self.type == 'Продольный эффект':
            m = 1
        return m

    def calc_Y(self, point_name):
        a = self.deffect.a
        c = self.deffect.c
        t = self.problem.t
        y_simple = (2 - 0.82 * a/c) / ((1 - (0.89 - 0.57*(a/c)**0.5)**3 * (a/t)**1.5)**3.25)
        if point_name == 'C':
            Y = y_simple * self.m_c() * self.gamma_c()
        elif point_name == 'D':
            Y = y_simple * self.m_d() * self.gamma_d()
        else:
            raise ValueError('Idk point name...')
        return Y


    def sig_m(self):
        if self.type == 'Кольцевой эффект':
            sig_m = self.problem.p * self.problem.Rin**2/(self.problem.Rout**2 - self.problem.Rin**2)
        elif self.type == 'Продольный эффект':
            sig_m = self.problem.p * self.problem.Rin / self.problem.t
        return sig_m

    @staticmethod
    def I_moment(radius, t):
        return math.pi * radius**3 * t

    def sig_b(self):
        """
        Интегральное
        """
        if self.type == 'Кольцевой эффект':
            #fun = lambda y, M, I: M/I * y
            #c = self.problem.Din/self.problem.Dout
            I = self.I_moment(self.problem.Rout-self.problem.t/2, self.problem.t)
            #sig_b = quad(fun, self.problem.Rin, self.problem.Rout, args=(self.problem.M*1e-6, I))[0]
            sig_b = self.problem.M*1e-6 / I * self.problem.Rout
        elif self.type == 'Продольный эффект':
            sig_b = 0
        return sig_b

    def sig_a(self):
        return self.sig_m()  + self.sig_b() * self.problem.Rin/self.problem.Rout

    def sig_c(self):
        return self.sig_m()  + self.sig_b() * (self.problem.Rin + self.deffect.a)/self.problem.Rout

    def calc_Sig_eq(self, point_name):
        a = self.deffect.a
        c = self.deffect.c
        t = self.problem.t

        sig_a = self.sig_a()
        sig_c = self.sig_c()

        if point_name == 'C':
            sig_eq = 0.39*sig_a + 0.61*sig_c - 0.11*a/c*(sig_a - sig_c) + 0.28*a/t*(1-(a/c)**0.5)*(sig_a - sig_c)
        elif point_name == 'D':
            sig_eq = 0.82*sig_a + 0.18*sig_c
        else:
            raise ValueError('Idk point name...')
        return sig_eq

    @staticmethod
    def calc_K(Y, Sig_eq, a):
        return Y * Sig_eq * a ** 0.5

    def changing_per_iter(self):
        res = {
            'C': None,
            'D': None,
            'delK_C': None,
            'delK_D': None,
        }
        for point in res.keys():
            if len(point) != 1:
                continue
            Sig_eq = self.calc_Sig_eq(point)
            Y = self.calc_Y(point)
            K = self.calc_K(Y, Sig_eq, self.deffect.a)
            change = C * K**self.problem.m
            res[point] = change
            res['delK_'+point] = K

        return res
