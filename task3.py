import numpy as np
import scipy
import seaborn as sns
from objects import Deffect, Steel, Problem, Solver, Find2cc, B_K_method

import pandas as pd
import matplotlib.pyplot as plt

def get_COA_COD_graphs(bkmethod, dva_c_max):
    x = np.linspace(0, dva_c_max/2, 20)
    y_COA = [bkmethod.get_COA(i) for i in x]
    y_COD = [bkmethod.get_COD(i) for i in x]

    plt.figure(figsize=(10, 8))
    result_df = pd.DataFrame({
        'Длина трещины, м': x,
        'COA': y_COA, 'COD': y_COD
        })
    plot = sns.lineplot(data=result_df, x='Длина трещины, м', y='COA')
    plot = sns.lineplot(data=result_df, x='Длина трещины, м', y='COD')
    plt.title(f'model: dva_c_max={dva_c_max}')
    plt.savefig(f"COA-COD.png")


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
    solver  = Solver(problem, deffect, steel, defect_type)
    find2cc = Find2cc(solver, deffect, problem, steel, problem_type)

    bkmethod = B_K_method(deffect, problem, solver, find2cc)
    print(bkmethod.get_COD(0))
    print(get_COA_COD_graphs(bkmethod, dva_c_max))

if __name__ == '__main__':
    main()
