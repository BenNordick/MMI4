import hiloop.multistability as multistability
import numpy as np
import pandas as pd
import sobol

def paramsdict(runner):
    return {p: runner[p] for p in runner.ps()}

def makeinitialconditions(runner, count=100, /, base=5., minexp=-4.5, maxexp=1.5, rna_prefix='R_', protein_prefix='X_', protein_param_prefix='l_'):
    exprange = maxexp - minexp
    rnas = sum(1 for f in runner.fs() if f.startswith(rna_prefix))
    rna_ini_combs = np.power(base, sobol.sample(rnas, count, skip=count) * exprange + minexp)
    ini_combs = np.zeros((rna_ini_combs.shape[0], len(runner.fs())))
    next_rna_col = 0
    rna_cols = dict()
    protein_cols = dict()
    for i, species in enumerate(runner.fs()):
        product = species.split('_')[-1]
        if species.startswith(rna_prefix):
            ini_combs[:, i] = rna_ini_combs[:, next_rna_col]
            next_rna_col += 1
            rna_cols[product] = i
        elif species.startswith(protein_prefix):
            protein_cols[product] = i
    def update_protein_inicombs():
        for product, col in protein_cols.items():
            ini_combs[:, col] = runner[f'{protein_param_prefix}{product}'] * ini_combs[:, rna_cols[product]]
    update_protein_inicombs()
    return ini_combs, update_protein_inicombs

def findsteadystate(runner, ini_comb, init_time, max_time, init_step, min_step):
    timestep = init_step
    time = init_time
    for iv, v in enumerate(runner.fs()):
        runner[v] = ini_comb[iv]
    while time < max_time:
        while True:
            points = round(time / timestep) + 1
            try:
                sim_points = runner.simulate(start=0, end=time, points=points)
                attractor = multistability.describe_attractor(sim_points, timestep, runner, False)
                if attractor is not None:
                    return attractor, timestep
                else:
                    break
            except RuntimeError as e:
                if 'Error test failures occurred too many times' in e.args[0] and timestep >= min_step * 5:
                    timestep /= 5
                else:
                    return None
        time *= 2
    return None, timestep

def findattractors(runner, ini_combs, init_time=100, max_time=1000, init_step=5, min_step=0.2):
    timestep = init_step
    attractors = []
    for ini_comb in ini_combs:
        result = findsteadystate(runner, ini_comb, init_time, max_time, timestep, min_step)
        if result is None:
            continue
        attractor, timestep = result
        if attractor is not None and all(not multistability.equivalent_attractors(attractor, e) for e in attractors):
            attractors.append(attractor)
    return attractors

def findpointattractors(*args, **kwargs):
    return [a for a in findattractors(*args, **kwargs) if not isinstance(a, dict)]

def pointattractorsdf(runner, attractors):
    return pd.DataFrame.from_records(attractors, columns=runner.fs())