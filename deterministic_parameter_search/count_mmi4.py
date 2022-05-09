import argparse
import collections
import json
import multiattractor
import numpy as np
import sobol
import tellurium as te

parser = argparse.ArgumentParser()
parser.add_argument('count', type=int, help='number of psets to try')
parser.add_argument('output', type=str, help='output counts file')
args = parser.parse_args()

runner = te.loada('''
JM: -> mRNA; k_ZEB1 - mRNA
JI: -> miRNA; k_miR - gm * miRNA

JC1: mRNA + miRNA -> C1; 4 * kOn * mRNA * miRNA - kOff * C1
JDM1: C1 -> miRNA; a1 * C1
JDI1: C1 -> mRNA; b1 * gm * C1

JC2: C1 + miRNA -> C2; 3 * kOn * C1 * miRNA - 2 * kOff * C2
JDM2: C2 -> 2 miRNA; a2 * C2
JDI2: C2 -> mRNA + miRNA; 2 * b2 * gm * C2

JC3: C2 + miRNA -> C3; 2 * kOn * C2 * miRNA - 3 * kOff * C3
JDM3: C3 -> 3 miRNA; a3 * C3
JDI3: C3 -> mRNA + 2 miRNA; 3 * b3 * gm * C3

JC4: C3 + miRNA -> C4; kOn * C3 * miRNA - 4 * kOff * C4
JDM4: C4 -> 4 miRNA; a4 * C4
JDI4: C4 -> mRNA + 3 miRNA; 4 * b4 * gm * C4

mRNA = 1; miRNA = 1; C1 = 0; C2 = 0; C3 = 0; C4 = 0
k_ZEB1 = 1; k_miR = 1; gm = 1
kOff = 100; K = 100000
kOn := kOff * K
a1 = 1; a2 = 1; a3 = 1; a4 = 1
b1 = 1; b2 = 1; b3 = 1; b4 = 1
''')

rna_ini_combs = np.power(5, sobol.sample(2, 100, skip=63) * 5 - 4)
ini_combs = np.zeros((rna_ini_combs.shape[0], len(runner.fs())))
next_rna_col = 0
for i, species in enumerate(runner.fs()):
    if 'RNA' in species:
        ini_combs[:, i] = rna_ini_combs[:, next_rna_col]
        next_rna_col += 1

results = collections.defaultdict(int)
for test in range(args.count):
    pset = {}
    for p in ['a1', 'a2', 'a3', 'a4', 'b1', 'b2', 'b3', 'b4']:
        pset[p] = 2 ** np.random.uniform(-3, 4)
    pset['K'] = 10 ** np.random.uniform(3, 6)
    pset['k_miR'] = 2 ** np.random.uniform(-4, 1)
    for (p, value) in pset.items():
        runner[p] = value
    attractors = multiattractor.findpointattractors(runner, ini_combs)
    results[len(attractors)] += 1
    if test % 1000 == 0 or test == args.count - 1:
        with open(args.output, 'w') as f:
            json.dump(results, f)
