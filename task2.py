import numpy as np
import scipy
import seaborn as sns
from objects import Deffect, Steel, Problem, Solver, Find2cc

import pandas as pd
import matplotlib.pyplot as plt

def find_2cc(steel, deffect, problem, problem_type, defect_type, alpha):
    solver  = Solver(problem, deffect, steel, defect_type)
    find2cc = Find2cc(solver, deffect, problem, steel, problem_type)
    dva_cc, usability  = find2cc.dva_cc(alpha)
    return solver.sig_b(), solver.sig_m(), dva_cc, usability

def find_alpha(steel, deffect, problem, effect_type, problem_type, dcc):
    solver  = Solver(problem, deffect, steel, effect_type)
    find2cc = Find2cc(solver, deffect, problem, steel, problem_type)
    dva_cc  = find2cc.find_a(dcc)
    return dva_cc

def calculate_interception(effect_type, problem_type, dva_cc, *args):
    """ Not using """
    i = 0
    while deffect.a < t and 2*deffect.c < dva_cc:
        solver  = Solver(problem, deffect, steel, effect_type)
        changes = solver.changing_per_iter()
        da, dc = changes['C'], changes['D']
        deffect.a += da
        deffect.c += dc
        i += 1
        print(da, dc)
    return (deffect.a, 2*deffect.c, i)

def calc_params_for_task2(steel, deffect, problem, defect_type, problem_type, alpha=1, dc_max=0.8):
    sig_b, sig_m, dva_cc, usability = find_2cc(steel, deffect, problem, problem_type, defect_type, alpha)
    points = []
    for dcc in np.linspace(dva_cc, dva_cc+0.2, num=20):
        alpha = find_alpha(steel, deffect, problem, defect_type, problem_type, dcc)
        points.append((dcc , alpha))
    return sig_b, sig_m, dva_cc, points, usability

def create_picture(points, dva_cc, dva_c_max, index):
    #print(points)
    x = [0] + [p[0] for p in points]
    y = [1] + [p[1] for p in points]
    plt.figure(figsize=(10, 8))
    df = pd.DataFrame({'2c': x, 'a//t': y})
    plot = sns.lineplot(data=df, x='2c', y='a//t', color='red')
    plot.set(xlim=(0, dva_c_max+0.2))
    plot.set(ylim=(0, 1.3))
    plt.axvline(x[1])
    plt.title(f'model: 2cc={dva_cc}')
    plt.savefig(f"figures2/picture{index+1}.png")
    plt.close()

def main(Dout=273e-3,
         t=16e-3,
         p=6.9e6,
         steel_type='Сталь 20',
         Nz=10e3,
         Mx=6.94e4,
         My=-13.507e3,
         dva_c_max=0.8,
         C=1.5e-10,
         m=3.1,
         Rp02_min=1.83e8,
         Rp02_max=2.08e8,
         Rm_min=3.66e8,
         Rm_max=4.65e8,
         E_module=1.83e11,
         mu=0.3,
         T=285,
         defect_type='Кольцевой дефект',
         problem_type='ППН',
         index=None):
    # start point data for deffect
    a_0 = 0.2*t
    c_0 = 0.5*t
    steel = Steel(C, m, T, Rp02_min, Rp02_max, Rm_min, Rm_max, E_module, mu, steel_type)
    deffect = Deffect(a_0, c_0)
    problem = Problem(deffect, steel, t, Dout, p, Nz=Nz, Mx=Mx, My=My)
    sig_b, sig_m, dva_cc, points, usability = calc_params_for_task2(steel, deffect, problem, defect_type, problem_type)
    if problem_type=='ППН' and defect_type=='Кольцевой дефект':
        create_picture(points, dva_cc, dva_c_max, index)
    return [sig_b, sig_m, dva_cc, usability]

if __name__ == '__main__':
    while 1:
        asked_question = input('Are u sure tpu wanted launch from main?[yN]')
        if asked_question == 'y':
            main()
            break
        elif asked_question == 'N':
            break
        print('pls, input correct answer')
