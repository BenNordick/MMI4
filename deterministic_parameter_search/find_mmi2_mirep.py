import argparse
import multiattractor
import numpy as np
import pandas as pd
import random
import sobol
import tellurium as te

parser = argparse.ArgumentParser()
parser.add_argument('count', type=int, help='number of psets to find')
parser.add_argument('output', type=str, help='output CSV')
args = parser.parse_args()

runner = te.loada('''
J0: -> mRNA; k_ZEB1 - 2 * kOn * mRNA * miRNA + kOff * C1 + b1 * gm * C1 - mRNA
//            -- TRANSCRIPTION CONTROL ------------------------------------   -- MMI -->
J1: -> miRNA; mu * ((1 - r_I) + r_I / (1 + (X_ZEB1 / K_I_ZEB1) ^ n_I_ZEB1)) - 2 * kOn * mRNA * miRNA + kOff * C1 + a1 * C1 - kOn * C1 * miRNA + 2 * kOff * C2 + 2 * a2 * C2 - gm * miRNA
J2: -> C1; 2 * kOn * mRNA * miRNA - kOff * C1 - a1 * C1 - b1 * gm * C1 - kOn * C1 * miRNA + 2 * kOff * C2 + 2 * b2 * gm * C2
J3: -> C2; kOn * C1 * miRNA - 2 * kOff * C2 - a2 * C2 - 2 * b2 * gm * C2
J4: -> X_ZEB1; l_ZEB1 * mRNA - X_ZEB1

mRNA = 1.0; miRNA = 1.0; C1 = 0.0; C2 = 0.0; X_ZEB1 = 0.0
// MMI
a1 = 5; b1 = 0.5; a2 = 15; b2 = 0.5
gm = 1.0; kOff = 100; K = 5000; k_ZEB1 = 2
kOn := K * kOff
// Repression
mu = 0.07; r_I = 0.93; K_I_ZEB1 = 0.9; n_I_ZEB1 = 5; l_ZEB1 = 2
''')

rna_ini_combs = np.power(5, sobol.sample(2, 100, skip=99) * 5 - 4)
ini_combs = np.zeros((rna_ini_combs.shape[0], len(runner.fs())))
next_rna_col = 0
mrna_col = -1
protein_col = -1
for i, species in enumerate(runner.fs()):
    if 'RNA' in species:
        ini_combs[:, i] = rna_ini_combs[:, next_rna_col]
        next_rna_col += 1
    if species == 'mRNA':
        mrna_col = i
    elif species.startswith('X_'):
        protein_col = i
ini_combs[:, protein_col] = runner['l_ZEB1'] * ini_combs[:, mrna_col]

results = []

while len(results) < args.count:
    pset = {}
    pset['k_ZEB1'] = 10 ** np.random.uniform(-0.5, 0.5)
    pset['l_ZEB1'] = 10 ** np.random.uniform(0, 1)
    pset['mu'] = 10 ** np.random.uniform(-1.5, 0)
    pset['K'] = 10 ** np.random.uniform(3, 5)
    pset['a1'] = 2 ** np.random.uniform(0, 4)
    pset['a2'] = min(pset['a1'] * 2 ** np.random.uniform(0, 3), 32)
    pset['b1'] = 2 ** np.random.uniform(-2, 0)
    pset['b2'] = pset['b1'] * 2 ** np.random.uniform(-2, 0)
    pset['K_I_ZEB1'] = 10 ** np.random.uniform(-2, 1)
    pset['n_I_ZEB1'] = random.choice([1, 2, 3, 4, 5, 6])
    pset['r_I'] = np.random.uniform(0.8, 0.99)
    for p, value in pset.items():
        runner[p] = value
    ini_combs[:, protein_col] = ini_combs[:, mrna_col] * pset['l_ZEB1']
    attractors = multiattractor.findpointattractors(runner, ini_combs)
    if len(attractors) == 3:
        results.append(pset)
        pd.DataFrame.from_dict(results).to_csv(args.output, index=False)
