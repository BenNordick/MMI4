import argparse
import itertools
import hiloop.multistability as multistability
import numpy as np
import pandas as pd
import tellurium as te

parser = argparse.ArgumentParser()
parser.add_argument('count', type=int, help='number of psets to find')
parser.add_argument('output', type=str, help='output CSV')
args = parser.parse_args()

runner = te.loada('''
J0: -> mRNA_A; (1 - r_A) + r_A * (Protein_B / K_A_B)^n_A_B / (1 + (Protein_B / K_A_B)^n_A_B) - mRNA_A
J2: -> mRNA_B; k_B * ((1 - r_B) + r_B * (Protein_A / K_B_A)^n_B_A / (1 + (Protein_A / K_B_A)^n_B_A)) - mRNA_B
J3: -> Protein_A; l_A * mRNA_A - Protein_A
J4: -> Protein_B; l_B * mRNA_B - Protein_B

mRNA_A = 1.0; mRNA_B = 0; Protein_A = 0; Protein_B = 0
r_A = 0.99; r_B = 0.99; K_A_B = 1; K_B_A = 1; n_A_B = 3; n_B_A = 3; k_B = 1
l_A = 1; l_B = 1
''')

pts1d = [0, 0.05, 0.5, 5]
ptsnd = np.tile(pts1d, 2).reshape(2, -1)
mrna_ini_combs = np.array(list(itertools.product(*ptsnd)))
ini_combs = np.zeros((mrna_ini_combs.shape[0], 4))
next_mrna_col = 0
mrna_cols = dict()
protein_cols = dict()
for i, species in enumerate(runner.fs()):
    if 'mRNA' in species:
        ini_combs[:, i] = mrna_ini_combs[:, next_mrna_col]
        next_mrna_col += 1
        mrna_cols[species[-1]] = i
    else:
        protein_cols[species[-1]] = i

domains = {'K': (0.1, 1.5), 'k': (1, 1.5), 'r': [0.9, 0.99], 'n': [0.5, 3.5], 'l': (2, 0.5)}

results = []
while len(results) < args.count:
    pset = dict()
    for p in runner.ps():
        domain = domains[p[0]]
        if isinstance(domain, list):
            value = np.random.uniform(low=domains[p[0]][0], high=domains[p[0]][1])
        else:
            mean, sd = domain
            value = np.random.lognormal(np.log(mean) - (sd ** 2) / 2, sd)
        if p[0] == 'n':
            value = np.round(value)
        runner[p] = value
        pset[p] = value
        if p[0] == 'l':
            ini_combs[:, protein_cols[p[-1]]] = ini_combs[:, mrna_cols[p[-1]]] * value
    attractors = []
    time = 100
    while len(attractors) == 0 and time < 10000:
        attractors = multistability.findattractors(runner, ini_combs, time, 5)
        time *= 2
    if len(attractors) == 2:
        #print(pd.DataFrame.from_records(attractors, columns=runner.fs()))
        #print(pset)
        results.append(pset)
        print(len(results), end='\r')

df = pd.DataFrame.from_dict(results)
df.to_csv(args.output, index=False)