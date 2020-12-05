import json
from input_params import *

from objects import Deffect, Steel, Problem, Solve
from input_params import *
"""
Deffect - оъект трещины
Problem - обхект состояния среды
Solve - объект поставновки и решения задачи
"""


def single_story(a_0, c_0, type_):
    logs = {}
    effect = {'a':None, 'c':None}
    # порождаем деффект
    deffect = Deffect(a_0, c_0)
    a, c = a_0, c_0
    logs = {'C': {}, 'D': {}, 'delK_C':{}, 'delK_D':{}}
    for cycle_name in ['НУЭ', 'ГИ1', 'ГИ2']:
        logs['C'][cycle_name] = [deffect.a] if cycle_name=='НУЭ' else []
        logs['D'][cycle_name] = [deffect.c] if cycle_name=='НУЭ' else []
        logs['delK_C'][cycle_name] = []
        logs['delK_D'][cycle_name] = []

        p = cycles_info[cycle_name]['p']
        T = cycles_info[cycle_name]['T']
        M = cycles_info[cycle_name]['M']
        N = cycles_info[cycle_name]['N']
        problem = Problem(deffect, t, Dout, p, T, M, C, m)

        a_was, c_was = deffect.a, deffect.c
        for cycle in range(N):
            solve = Solve(problem, deffect, type_).changing_per_iter()
            deffect.a += solve['C']
            deffect.c += solve['D']
            logs['C'][cycle_name].append(deffect.a)
            logs['D'][cycle_name].append(deffect.c)
            logs['delK_C'][cycle_name].append(solve['delK_C'])
            logs['delK_D'][cycle_name].append(solve['delK_D'])

            if deffect.a > t:
                break
        a_became, c_became = deffect.a, deffect.c

        effect['a'] = a_became
        effect['c'] = c_became
    return effect, logs

def few_stories(a_0, c_0):
    logs = {}
    effect = {
        'Кольцевой эффект':  {None},
        'Продольный эффект': {None}
    }
    for type_ in effect.keys():
        effect[type_], logs[type_] = single_story(a_0, c_0, type_)
    return effect, logs

def while_loop_till_destroy(a_0, c_0, type_):
    a, c = a_0, c_0
    number_of_cycles = 0
    a_story, c_story = [], []
    while a < t:
        number_of_cycles += 1
        effect, logs = single_story(a, c, type_)
        a, c = effect['a'], effect['c']
        a_story.append(a)
        c_story.append(c)
    return a, c, number_of_cycles, a_story, c_story

if __name__ == "__main__":
    # фиксация изменений при испытании за одну историю нагружений
    print(f'Before iteration a and c are: {a_0}, {c_0}')
    effect, logs = few_stories(a_0, c_0)
    for problem_type, changes in effect.items():
        print(f"After iteration for problem type {problem_type} a and 2c became: {changes['a']}, {2*changes['c']}")

    # фиксируем в логах для визуализации
    with open('data/iter_logs.json', 'w') as output_file:
        json.dump(logs, output_file, ensure_ascii=False)

    # определяем максмальное число историй нагружений до сквозного пророста трещины
    stories = {}
    for effect_type in ['Кольцевой эффект', 'Продольный эффект']:
        a, c, number_of_cycles, a_story, c_story = while_loop_till_destroy(a_0, c_0, effect_type)
        print(f'[{effect_type}]Pipe will be broken with: a={a}, 2c={2*c}', end=' ')
        print(f'on number_of_cycles = {number_of_cycles}')
        stories[effect_type] = {'a_story': a_story, 'c_story': c_story}

    # фиксируем в логах для визуализации
    with open(f'data/stories_logs.json', 'w') as output_file:
        json.dump(stories, output_file, ensure_ascii=False)
