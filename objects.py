import math
from scipy.integrate import quad
import numpy as np
import sympy
from scipy.optimize import fsolve, bisect, newton_krylov, diagbroyden

import warnings
warnings.filterwarnings("error")

class Deffect:
    """
    Состояние трещины
    """
    def __init__(self, a, c):
        self.a = a # m
        self.c = c # m

    @property
    def double_c(self):
        return 2 * self.c # m

class Steel:
    """
    Параметры стали
    """
    def __init__(self, C, m, T, Rp02_min, Rp02_max, Rm_min, Rm_max, E_module, mu, name):
        self.C  = C
        self.m  = m
        self.mu = mu
        self.Rp02_min = Rp02_min # Pa
        self.Rp02_max = Rp02_max # Pa
        self.Rm_min   = Rm_min   # Pa
        self.Rm_max   = Rm_max   # Pa
        self.E_module = E_module # Pa
        self.name = name
        self.T = T # cel

    @property
    def ro(self):
        if self.name == "Сталь 20":
            if self.T == 20:
                return 7856 # kg/m3
            elif self.T == 150:
                return 7819 # kg/m3
            elif self.T == 285:
                return 7775 # kg/m3
        elif self.name == "Сталь 16ГС":
           return 7850
        raise ValueError('Cant calculate the density.')

    @property
    def structure_type(self):
        if self.name == "Сталь 20":
            return None
        elif self.name == "Сталь 16ГС":
            return None
        raise ValueError('Cant define structure.')


class Problem:
    """
    Состояние системы в целом. Параметры окружающей среды и трещина
    """
    def __init__(self, deffect, steel, t, Dout, p, N=None, Nx=0, Ny=0, Nz=0, M=None, Mx=0, My=0, Mz=0,):
        self.deffect = deffect
        self.steel = steel
        self.t = t # m
        self.Dout = Dout # m
        self.Din  = Dout - 2*t # m
        self.p = p # Pa
        self.M = M if M else (Mx**2 + My**2 + Mz**2)**0.5 # N*m
        self.N = N if N else (Nx**2 + Ny**2 + Nz**2)**0.5 # N*m
        self.C = steel.C
        self.m = steel.m

    @property
    def Rout(self):
        return self.Dout/2

    @property
    def Rin(self):
        return self.Din/2

    @property
    def Rm(self):
        return (self.Rin + self.Rout)/2

class Solver:
    """
    Объект решения
    """
    def __init__(self, problem, deffect, steel, type_):
        self.problem = problem
        self.deffect = deffect
        self.steel = steel
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
        ro = self.steel.ro
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
        if self.type == 'Кольцевой дефект':
            section = math.pi * (self.problem.Rout ** 2 - self.problem.Rin ** 2)
            sig_m = self.problem.p * self.problem.Rin**2/(self.problem.Rout**2 - self.problem.Rin**2) + self.problem.N / section
        elif self.type == 'Продольный дефект':
            sig_m = self.problem.p * self.problem.Rin / self.problem.t
        else:
            raise ValueError('WRONG DEFFECT NAME!')
        return sig_m

    @staticmethod
    def I_moment(radius, t):
        return math.pi * radius**3 * t

    def sig_b(self):
        """
        Интегральное
        """
        if self.type == 'Кольцевой дефект':
            I = self.I_moment(self.problem.Rout-self.problem.t/2, self.problem.t)
            sig_b = self.problem.M/ I * self.problem.Rout
        elif self.type == 'Продольный дефект':
            sig_b = 0
        else:
            raise ValueError('WRONG DEFFECT NAME!')
        return sig_b

    def sig_a(self):
        sig_a = self.sig_m()  + self.sig_b() * self.problem.Rin/self.problem.Rout
        return sig_a

    def sig_c(self):
        sig_c = self.sig_m()  + self.sig_b() * (self.problem.Rin + self.deffect.a)/self.problem.Rout
        return sig_c

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
        print(Y, Sig_eq, a)
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
            K = self.calc_K(Y, Sig_eq*1e-6, self.deffect.a)
            change = self.steel.C * (K)**self.problem.m
            res[point] = change
            res['delK_'+point] = K
            print(res)
        return res

class Find2cc:

    def __init__(self, Solver, Deffect, Problem, Steel, problem_type):
        self.solver = Solver
        self.problem = Problem
        self.steel = Steel
        self.deffect = Deffect
        self.problem_type = problem_type

    @property
    def sig_f(self):
        deffect_type = self.solver.type
        if deffect_type == 'Кольцевой дефект':
            return 0.42*(self.steel.Rp02_min + self.steel.Rm_min)
        if deffect_type == 'Продольный дефект':
            if self.steel.name in ['Сталь 20', 'Сталь 16ГС']:
                return 0.5*(self.steel.Rp02_min + self.steel.Rm_min)
        raise ValueError('Steel not found')

    def get_solution(self, fun):
        """ Solve equation """
        init_guess = 0.3
        a, b = 0, 1.5
        try:
            interested_variable = fsolve(fun, init_guess)[0]
            #interested_variable = newton_krylov(fun, init_guess, verbose=False)
            #interested_variable = diagbroyden(fun, init_guess, verbose=False)
            #interested_variable = round(interested_variable, 3)
        except (RuntimeWarning, ValueError, Exception)  as e:
            interested_variable = None
        return interested_variable

    def _tetta(self, dva_cc):
        return dva_cc/(2*self.problem.Rin)

    def _betta(self, dva_cc, alpha, sig_m, sig_f):
        return math.pi/2 * (1 - alpha * self._tetta(dva_cc)/np.pi - sig_m/sig_f)

    def _ka(self, dva_cc, f):
        tetta = self._tetta(dva_cc)
        ka_nominator = lambda dva_cc: 1 - f*(tetta/np.pi + np.sin(2*tetta)/(2*np.pi) - 2*np.sin(tetta)/np.pi*np.cos(tetta))
        denominator =  lambda dva_cc: (1 - f * tetta/np.pi)*(1 - f*(tetta/np.pi + np.sin(2*tetta)/(2*np.pi))) - (f**2)*(2*(np.sin(tetta)**2/(np.pi**2)))
        return ka_nominator(dva_cc)/denominator(dva_cc)

    def _kb(self, dva_cc, f):
        tetta = self._tetta(dva_cc)
        kb_nominator = lambda dva_cc: np.sin(tetta)/np.pi * f + (1 - f*tetta/np.pi)*np.cos(tetta)
        denominator =  lambda dva_cc: (1 - f * tetta/np.pi)*(1 - f*(tetta/np.pi + np.sin(2*tetta)/(2*np.pi)) - (f**2)*(2*(np.sin(tetta)**2/(np.pi**2))))
        return kb_nominator(dva_cc)/denominator(dva_cc)

    def check_possibility_to_use_formula(self, tetta, betta):
        bool1 = betta >= 0 and betta <= np.pi / 2
        bool2 = (tetta <= np.pi/2) if self.problem_type == 'ППН' else True
        bool3 = (tetta + betta)<=np.pi if self.problem_type == 'ППН' else (-tetta + betta)<=np.pi
        bool4 = (2*tetta <= np.pi) if (self.solver.sig_m()==0 or self.solver.sig_b()==0) else True
        return (bool1 and bool2 and bool3 and bool4)

    def dva_cc(self, alpha):
        sig_f = self.sig_f
        sig_m = self.solver.sig_m()
        sig_b = self.solver.sig_b()
        init_guess = 0.3
        if self.problem_type == "ППН":
            func = lambda dva_cc: 2/np.pi * sig_f * (2 * np.sin(self._betta(dva_cc, alpha, sig_m, sig_f)) - alpha * np.sin(self._tetta(dva_cc))) - sig_b
            dva_cc_solution = self.get_solution(func)
        elif self.problem_type == "ЛРН":
            func = lambda dva_cc: self._ka(dva_cc, alpha)*sig_m + self._kb(dva_cc, alpha)*sig_b - sig_f
            dva_cc_solution = self.get_solution(func)
        else:
            raise ValueError('problem type not exists')

        check_usability = None
        if dva_cc_solution:
            check_usability = self.check_possibility_to_use_formula(self._tetta(dva_cc_solution), self._betta(dva_cc_solution, alpha, sig_m, sig_f))
        return dva_cc_solution, check_usability

    def find_a(self, dva_cc):
        sig_f = self.sig_f
        Rin = self.problem.Rin
        sig_m = self.solver.sig_m()
        sig_b = self.solver.sig_b()
        func = lambda a: 2/math.pi * sig_f * (2 * math.sin(math.pi/2 * (1 + a * - dva_cc/(2*Rin*math.pi) - sig_m/sig_f)) - a * math.sin(dva_cc/(2*Rin)) ) - sig_b
        a = self.get_solution(func)
        return a

class B_K_method:
    E = 1.83e11
    Mu = 0.3
    def __init__(self, deffect, problem, solver, find2cc):
        self.deffect = deffect
        self.problem = problem
        self.solver = solver
        self.find2cc = find2cc
        self.sig_m = self.solver.sig_m()
        self.sig_b = self.solver.sig_b()

    def get_E(self):
        return self.E

    def get_Mu(self):
        return self.Mu

    def get_sig_f(self):
        return self.find2cc.sig_f

    @property
    def sig_eqv(self):
        return self.sig_m + self.sig_b

    @property
    def gamma(self):
        x = 1/2**0.5 * (self.sig_eqv / self.find2cc.sig_f)
        return (1-x**3) / (1-x**2)**2

    def lambda_(self, c):
        return (12*(1 - self.Mu**2))**0.25 * c / (self.problem.Rm * self.problem.t)**0.5

    def alpha(self, c):
        lambda_ = self.lambda_(c)
        if lambda_ <= 5:
            return (1 + 0.117*lambda_**2)**0.5
        elif lambda_ <= 8 and lambda_ > 5:
            return 1 + 0.1*lambda_ + 0.16*lambda_**2
        else:
            return 1 + 0.1*lambda_ + 0.16*lambda_**2
            raise ValueError(f'lambda: {lambda_}>8 ')

    def A0(self, c):
        return 7.54 * (self.sig_eqv / self.E) * c**2

    def get_COA(self, c):
        if c == 0:
            c = 1e-15
        #print(self.alpha(c), self.gamma, self.A0(c))
        return self.alpha(c) * self.gamma * self.A0(c)

    def get_COD(self, c):
        if c == 0:
            c = 1e-15
        COA = self.get_COA(c)
        return 2*COA/(np.pi*c)
