import argparse
import multiattractor
import numpy as np
import pandas as pd
import random
import tellurium as te

parser = argparse.ArgumentParser()
parser.add_argument('count', type=int, help='number of psets to find')
parser.add_argument('output', type=str, help='output CSV')
parser.add_argument('--realistic', action='store_true', help='biological MMI rate progression')
args = parser.parse_args()

runner = te.loada('''
JR: -> R_M; k_ZEB1 - R_M
JI1: -> R_1; k_miR1 * ((1 - r_1) + r_1 / (1 + (X_M / K_1_M) ^ n_1_M)) - gm_1 * R_1
JI2: -> R_2; k_miR2 * ((1 - r_2) + r_2 / (1 + (X_M / K_2_M) ^ n_2_M)) - gm_2 * R_2

JP: -> X_M; l_M * R_M - X_M

JC0_1: R_M + R_1 -> C1; 2 * kOn_1 * R_M * R_1 - kOff_1 * C1
JC1_1: C1 + R_1 -> C11; kOn_1 * C1 * R_1 - 2 * kOff_1 * C11
JC0_2: R_M + R_2 -> C2; 2 * kOn_2 * R_M * R_2 - kOff_2 * C2
JC2_2: C2 + R_2 -> C22; kOn_2 * C2 * R_2 - 2 * kOff_2 * C22

JC1_2: C1 + R_2 -> C12; 2 * kOn_2 * C1 * R_2 - kOff_2 * C12
JC12_2: C12 + R_2 -> C122; kOn_2 * C12 * R_2 - 2 * kOff_2 * C122
JC11_2: C11 + R_2 -> C112; 2 * kOn_2 * C11 * R_2 - kOff_2 * C112
JC112_2: C112 + R_2 -> C1122; kOn_2 * C112 * R_2 - 2 * kOff_2 * C1122

JC2_1: C2 + R_1 -> C12; 2 * kOn_1 * C2 * R_1 - kOff_1 * C12
JC12_1: C12 + R_1 -> C112; kOn_1 * C12 * R_1 - 2 * kOff_1 * C112
JC22_1: C22 + R_1 -> C122; 2 * kOn_1 * C22 * R_1 - kOff_1 * C122
JC122_1: C122 + R_1 -> C1122; kOn_1 * C122 * R_1 - 2 * kOff_1 * C1122

JC1_DM: C1 -> R_1; a_1 * C1
JC1_D1: C1 -> R_M; b_1 * gm_1 * C1

JC2_DM: C2 -> R_2; a_2 * C2
JC2_D2: C2 -> R_M; b_2 * gm_2 * C2

JC11_DM: C11 -> 2 R_1; a_11 * C11
JC11_D1: C11 -> C1; 2 * b_11 * gm_1 * C11

JC22_DM: C22 -> 2 R_2; a_22 * C22
JC22_D2: C22 -> C2; 2 * b_22 * gm_2 * C22

JC12_DM: C12 -> R_1 + R_2; a_12 * C12
JC12_D1: C12 -> C2; b_12_1 * gm_1 * C12
JC12_D2: C12 -> C1; b_12_2 * gm_2 * C12

JC112_DM: C112 -> 2 R_1 + R_2; a_112 * C112
JC112_D1: C112 -> C12; 2 * b_112_1 * gm_1 * C112
JC112_D2: C112 -> C11; b_112_2 * gm_2 * C112

JC122_DM: C122 -> R_1 + 2 R_2; a_122 * C122
JC122_D1: C122 -> C22; b_122_1 * gm_1 * C122
JC122_D2: C122 -> C12; 2 * b_122_2 * gm_2 * C122

JC1122_DM: C1122 -> 2 R_1 + 2 R_2; a_1122 * C1122
JC1122_D1: C1122 -> C122; 2 * b_1122_1 * gm_1 * C1122
JC1122_D2: C1122 -> C112; 2 * b_1122_2 * gm_2 * C1122

R_M = 1; R_1 = 0.1; R_2 = 0.1
C1 = 0; C2 = 0; C11 = 0; C22 = 0
C12 = 0; C112 = 0; C122 = 0; C1122 = 0
X_M = 1

k_ZEB1 = 1; k_miR1 = 0.1; k_miR2 = 0.1
gm_1 = 1; gm_2 = 1
kOff_1 = 100; kOff_2 = 100
K_1 = 50000; K_2 = 50000
kOn_1 := kOff_1 * K_1; kOn_2 := kOff_2 * K_2

a_1 = 1; b_1 = 1; a_11 = 1; b_11 = 1
a_2 = 1; b_2 = 1; a_22 = 1; b_22 = 1
a_12 = 1; b_12_1 = 1; b_12_2 = 1
a_112 = 1; b_112_1 = 1; b_112_2 = 1
a_122 = 1; b_122_1 = 1; b_122_2 = 1
a_1122 = 1; b_1122_1 = 1; b_1122_2 = 1

l_M = 1
r_1 = 0.99; K_1_M = 0.1; n_1_M = 1
r_2 = 0.99; K_2_M = 0.1; n_2_M = 1
''')

ini_combs, update_protein_inicombs = multiattractor.makeinitialconditions(runner)

results = []
while len(results) < args.count:
    pset = {}
    pset['k_ZEB1'] = 10 ** np.random.uniform(-0.5, 0.5)
    pset['l_M'] = 10 ** np.random.uniform(0, 1)
    for p in runner.ps():
        if not args.realistic:
            if p.startswith('a_'):
                pset[p] = 2 ** np.random.uniform(-1, 4)
            elif p.startswith('b_'):
                pset[p] = 2 ** np.random.uniform(-4, 1)
        if p.startswith('K_') and len(p) == 3:
            pset[p] = 10 ** np.random.uniform(3, 5)
        elif p.startswith('k_miR'):
            pset[p] = 10 ** np.random.uniform(-1.5, 0)
    if args.realistic:
        pset['a_1'] = 2 ** np.random.uniform(0, 4)
        pset['a_11'] = min(pset['a_1'] * 2 ** np.random.uniform(0, 3), 32)
        pset['a_2'] = 2 ** np.random.uniform(0, 4)
        pset['a_22'] = min(pset['a_2'] * 2 ** np.random.uniform(0, 3), 32)
        pset['a_12'] = min(max(pset['a_1'], pset['a_2']) * 2 ** np.random.uniform(0, 2), 32)
        pset['a_112'] = min(max(pset['a_11'], pset['a_2']) * 2 ** np.random.uniform(0, 2), 32)
        pset['a_122'] = min(max(pset['a_1'], pset['a_22']) * 2 ** np.random.uniform(0, 2), 32)
        pset['a_1122'] = min(max(pset['a_11'], pset['a_22']) * 2 ** np.random.uniform(0, 1), 32)
        pset['b_1'] = 2 ** np.random.uniform(-2, 0)
        pset['b_11'] = pset['b_1'] * 2 ** np.random.uniform(-3, 0)
        pset['b_2'] = 2 ** np.random.uniform(-2, 0)
        pset['b_22'] = pset['b_2'] * 2 ** np.random.uniform(-3, 0)
        pset['b_12_1'] = min(pset['b_1'] * 2 ** np.random.uniform(-1, 1), 1)
        pset['b_12_2'] = min(pset['b_2'] * 2 ** np.random.uniform(-1, 1), 1)
        pset['b_112_1'] = min(pset['b_11'] * 2 ** np.random.uniform(-1, 1), 1)
        pset['b_112_2'] = min(pset['b_2'] * 2 ** np.random.uniform(-1, 1), 1)
        pset['b_122_1'] = min(pset['b_1'] * 2 ** np.random.uniform(-1, 1), 1)
        pset['b_122_2'] = min(pset['b_22'] * 2 ** np.random.uniform(-1, 1), 1)
        pset['b_1122_1'] = max(min(pset['b_112_1'] * 2 ** np.random.uniform(-1, 1), 1), 1/32)
        pset['b_1122_2'] = max(min(pset['b_122_2'] * 2 ** np.random.uniform(-1, 1), 1), 1/32)
    for mi in ['1', '2']:
        pset[f'K_{mi}_M'] = 10 ** np.random.uniform(-3, 1)
        pset[f'n_{mi}_M'] = random.choice([1, 2, 3, 4])
        pset[f'r_{mi}'] = np.random.uniform(0.8, 0.99)
    for (p, value) in pset.items():
        runner[p] = value
    update_protein_inicombs()
    attractors = multiattractor.findpointattractors(runner, ini_combs)
    if len(attractors) >= 4:
        # print(pset)
        # print(multiattractor.pointattractorsdf(runner, attractors))
        pset['attractors'] = len(attractors)
        results.append(pset)
        pd.DataFrame.from_dict(results).to_csv(args.output, index=False)