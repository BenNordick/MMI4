import argparse
import multiattractor
import numpy as np
import pandas as pd
import sobol
import tellurium as te

parser = argparse.ArgumentParser()
parser.add_argument('count', type=int, help='number of psets to find')
parser.add_argument('output', type=str, help='output CSV')
parser.add_argument('--realistic', action='store_true', help='biological MMI rate progression')
args = parser.parse_args()

runner = te.loada('''
JR: -> mRNA; k_ZEB1 - mRNA
JI1: -> miRNA1; k_miR1 - gm_1 * miRNA1
JI2: -> miRNA2; k_miR2 - gm_2 * miRNA2

JC0_1: mRNA + miRNA1 -> C1; 2 * kOn_1 * mRNA * miRNA1 - kOff_1 * C1
JC1_1: C1 + miRNA1 -> C11; kOn_1 * C1 * miRNA1 - 2 * kOff_1 * C11
JC0_2: mRNA + miRNA2 -> C2; 2 * kOn_2 * mRNA * miRNA2 - kOff_2 * C2
JC2_2: C2 + miRNA2 -> C22; kOn_2 * C2 * miRNA2 - 2 * kOff_2 * C22

JC1_2: C1 + miRNA2 -> C12; 2 * kOn_2 * C1 * miRNA2 - kOff_2 * C12
JC12_2: C12 + miRNA2 -> C122; kOn_2 * C12 * miRNA2 - 2 * kOff_2 * C122
JC11_2: C11 + miRNA2 -> C112; 2 * kOn_2 * C11 * miRNA2 - kOff_2 * C112
JC112_2: C112 + miRNA2 -> C1122; kOn_2 * C112 * miRNA2 - 2 * kOff_2 * C1122

JC2_1: C2 + miRNA1 -> C12; 2 * kOn_1 * C2 * miRNA1 - kOff_1 * C12
JC12_1: C12 + miRNA1 -> C112; kOn_1 * C12 * miRNA1 - 2 * kOff_1 * C112
JC22_1: C22 + miRNA1 -> C122; 2 * kOn_1 * C22 * miRNA1 - kOff_1 * C122
JC122_1: C122 + miRNA1 -> C1122; kOn_1 * C122 * miRNA1 - 2 * kOff_1 * C1122

JC1_DM: C1 -> miRNA1; a_1 * C1
JC1_D1: C1 -> mRNA; b_1 * gm_1 * C1

JC2_DM: C2 -> miRNA2; a_2 * C2
JC2_D2: C2 -> mRNA; b_2 * gm_2 * C2

JC11_DM: C11 -> 2 miRNA1; a_11 * C11
JC11_D1: C11 -> C1; b_11 * gm_1 * C11

JC22_DM: C22 -> 2 miRNA2; a_22 * C22
JC22_D2: C22 -> C2; b_22 * gm_2 * C22

JC12_DM: C12 -> miRNA1 + miRNA2; a_12 * C12
JC12_D1: C12 -> C2; b_12_1 * gm_1 * C12
JC12_D2: C12 -> C1; b_12_2 * gm_2 * C12

JC112_DM: C112 -> 2 miRNA1 + miRNA2; a_112 * C112
JC112_D1: C112 -> C12; b_112_1 * gm_1 * C112
JC112_D2: C112 -> C11; b_112_2 * gm_2 * C112

JC122_DM: C122 -> miRNA1 + 2 miRNA2; a_122 * C122
JC122_D1: C122 -> C22; b_122_1 * gm_1 * C122
JC122_D2: C122 -> C12; b_122_2 * gm_2 * C122

JC1122_DM: C1122 -> 2 miRNA1 + 2 miRNA2; a_1122 * C1122
JC1122_D1: C1122 -> C122; b_1122_1 * gm_1 * C1122
JC1122_D2: C1122 -> C112; b_1122_2 * gm_2 * C1122

mRNA = 1; miRNA1 = 0.1; miRNA2 = 0.1
C1 = 0; C2 = 0; C11 = 0; C22 = 0
C12 = 0; C112 = 0; C122 = 0; C1122 = 0

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
''')

rna_ini_combs = np.power(5, sobol.sample(3, 100, skip=99) * 5 - 4)
ini_combs = np.zeros((rna_ini_combs.shape[0], len(runner.fs())))
next_rna_col = 0
for i, species in enumerate(runner.fs()):
    if 'RNA' in species:
        ini_combs[:, i] = rna_ini_combs[:, next_rna_col]
        next_rna_col += 1

results = []
while len(results) < args.count:
    pset = {}
    pset['k_ZEB1'] = 10 ** np.random.uniform(-0.5, 0.5)
    for p in runner.ps():
        if not args.realistic:
            if p.startswith('a_'):
                pset[p] = 2 ** np.random.uniform(-1, 4)
            elif p.startswith('b_'):
                pset[p] = 2 ** np.random.uniform(-4, 1)
        if p.startswith('K_'):
            pset[p] = 10 ** np.random.uniform(3, 5)
        elif p.startswith('k_miR'):
            pset[p] = 10 ** np.random.uniform(-1.5, -0.5)
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
        pset['b_11'] = pset['b_1'] * 2 ** np.random.uniform(-2, 0)
        pset['b_2'] = 2 ** np.random.uniform(-2, 0)
        pset['b_22'] = pset['b_2'] * 2 ** np.random.uniform(-2, 0)
        pset['b_12_1'] = min(pset['b_1'] * 2 ** np.random.uniform(-1, 1), 1)
        pset['b_12_2'] = min(pset['b_2'] * 2 ** np.random.uniform(-1, 1), 1)
        pset['b_112_1'] = min(pset['b_11'] * 2 ** np.random.uniform(-1, 1), 1)
        pset['b_112_2'] = min(pset['b_2'] * 2 ** np.random.uniform(-1, 1), 1)
        pset['b_122_1'] = min(pset['b_22'] * 2 ** np.random.uniform(-1, 1), 1)
        pset['b_122_2'] = min(pset['b_1'] * 2 ** np.random.uniform(-1, 1), 1)
        pset['b_1122_1'] = max(min(pset['b_112_1'] * 2 ** np.random.uniform(-1, 1), 1), 1/32)
        pset['b_1122_2'] = max(min(pset['b_122_2'] * 2 ** np.random.uniform(-1, 1), 1), 1/32)
    for (p, value) in pset.items():
        runner[p] = value
    attractors = multiattractor.findpointattractors(runner, ini_combs)
    if len(attractors) >= 3:
        # print(pset)
        # print(multiattractor.pointattractorsdf(runner, attractors))
        pset['attractors'] = len(attractors)
        results.append(pset)
        if len(results) % 5 == 0 or len(results) == args.count:
            pd.DataFrame.from_dict(results).to_csv(args.output, index=False)