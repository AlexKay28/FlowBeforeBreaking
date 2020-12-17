import numpy as np
import scipy
import seaborn as sns
from objects import Deffect, Steel, Problem, Solver, Find2cc, B_K_method, Flow_Q

import pandas as pd
import matplotlib.pyplot as plt

def main(Dout=273e-3,
        t=16e-3,
        p=6.9e6,
        steel_type='Сталь 20',
        Nz=13e3,
        Mx=7.153e3,
        My=-9.608e3,
        dva_c_max=0.43,
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
    a_0 = t # трещина сквозная
    c_0 = 0

    steel = Steel(C, m, T, Rp02_min, Rp02_max, Rm_min, Rm_max, E_module, mu, steel_type)
    deffect = Deffect(a_0, c_0)
    problem = Problem(deffect, steel, t, Dout, p, Nz=Nz, Mx=Mx, My=My)
    solver  = Solver(problem, deffect, steel, defect_type)
    find2cc = Find2cc(solver, deffect, problem, steel, problem_type)

    x = np.linspace(0, 0.2, 15)
    y_COA, y_COD = [], []
    for i in x:
        steel = Steel(C, m, T, Rp02_min, Rp02_max, Rm_min, Rm_max, E_module, mu, steel_type)
        deffect = Deffect(a_0, i)
        problem = Problem(deffect, steel, t, Dout, p, Nz=Nz, Mx=Mx, My=My)
        solver  = Solver(problem, deffect, steel, defect_type)
        find2cc = Find2cc(solver, deffect, problem, steel, problem_type)
        bkmethod = B_K_method(deffect, problem, solver, find2cc)

        print(i)

        y_COA.append(bkmethod.get_COA())
        y_COD.append(bkmethod.get_COD())
        print(y_COA[-1], y_COD[-1])

        flow = Flow_Q(deffect, problem, solver, find2cc, bkmethod)
        print('Qc:', flow.get_Qc() * 60, 'Qld:', flow.get_Qld(1.9, 5))

    # plt.figure(figsize=(10, 8))
    # result_df = pd.DataFrame({
    #     'Длина трещины, м': x,
    #     'COA': y_COA, 'COD': y_COD
    #     })
    # sns.lineplot(data=result_df, x='Длина трещины, м', y='COA')
    # sns.lineplot(data=result_df, x='Длина трещины, м', y='COD')
    # plt.legend(('COA', 'COD'))
    # plt.title(f'model: dva_c_max={dva_c_max}')
    # plt.ylabel('Площадь раскрытия (м2)')
    # plt.savefig(f"COA-COD.png")




if __name__ == '__main__':
    main()