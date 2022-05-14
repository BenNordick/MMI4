import argparse
import multiattractor
import numpy as np
import pandas as pd
import sobol
import tellurium as te

parser = argparse.ArgumentParser()
parser.add_argument('count', type=int, help='number of psets to find')
parser.add_argument('output', type=str, help='output CSV')
args = parser.parse_args()

runner = te.loada('''
function Hminus(B, B0, n)
    1 / (1 + (B / B0)^n)
end
function Hs(B, B0, n, fc)
    Hminus(B, B0, n) + fc * (1 - Hminus(B, B0, n))
end

JI101:  -> R_101;   k_101 * Hs(X_SNAI, K_101_SNAI, n_101_SNAI, fc_101_SNAI) * Hs(X_SLUG, K_101_SLUG, n_101_SLUG, fc_101_SLUG) - d_101 * R_101
JI200:  -> R_200;   k_200 * Hs(X_ZEB, K_200_ZEB, n_200_ZEB, fc_200_ZEB) * Hs(X_SNAI, K_200_SNAI, n_200_SNAI, fc_200_SNAI) * Hs(X_SLUG, K_200_SLUG, n_200_SLUG, fc_200_SLUG) - d_200 * R_200
JRZ:    -> R_ZEB;   k_ZEB * Hs(X_ZEB, K_ZEB_ZEB, n_ZEB_ZEB, fc_ZEB_ZEB) * Hs(X_SNAI, K_ZEB_SNAI, n_ZEB_SNAI, fc_ZEB_SNAI) * Hs(X_SLUG, K_ZEB_SLUG, n_ZEB_SLUG, fc_ZEB_SLUG) - dR_ZEB * R_ZEB
JRSn:   -> R_SNAI;  k_SNAI * Hs(X_SNAI, K_SNAI_SNAI, n_SNAI_SNAI, fc_SNAI_SNAI) * Hs(X_SLUG, K_SNAI_SLUG, n_SNAI_SLUG, fc_SNAI_SLUG) - dR_SNAI * R_SNAI
JRSl:   -> R_SLUG;  k_SLUG * Hs(X_SLUG, K_SLUG_SLUG, n_SLUG_SLUG, fc_SLUG_SLUG) * Hs(X_SNAI, K_SLUG_SNAI, n_SLUG_SNAI, fc_SLUG_SNAI) - dR_SLUG * R_SLUG

JPZ:    -> X_ZEB;   l_ZEB * (R_ZEB + t_ZEB_2 * C_ZEB_2 + t_ZEB_22 * C_ZEB_22 + t_ZEB_1 * C_ZEB_1 + t_ZEB_12 * C_ZEB_12 + t_ZEB_122 * C_ZEB_122) - dP_ZEB * X_ZEB
JPSn:   -> X_SNAI;  l_SNAI * R_SNAI - dP_SNAI * X_SNAI
JPSl:   -> X_SLUG;  l_SLUG * (R_SLUG + t_SLUG * C1_SLUG) - dP_SLUG * X_SLUG

JCZ2:       R_ZEB + R_200 -> C_ZEB_2;       2 * kOn_ZEB_200 * R_ZEB * R_200 - kOff * C_ZEB_2
JCZ2DR:     C_ZEB_2 -> R_200;               dR_ZEB * a_ZEB_2 * C_ZEB_2
JCZ2DI:     C_ZEB_2 -> R_ZEB;               d_200 * b_200_ZEB_2 * C_ZEB_2
JCZ22:      C_ZEB_2 + R_200 -> C_ZEB_22;    kOn_ZEB_200 * C_ZEB_2 * R_200 - 2 * kOff * C_ZEB_22
JCZ22DR:    C_ZEB_22 -> 2 R_200;            dR_ZEB * a_ZEB_22 * C_ZEB_22
JCZ22DI:    C_ZEB_22 -> C_ZEB_2;            2 * d_200 * b_200_ZEB_22 * C_ZEB_22
JCZ1:       R_ZEB + R_101 -> C_ZEB_1;       kOn_ZEB_101 * R_ZEB * R_101 - kOff * C_ZEB_1
JCZ1DR:     C_ZEB_1 -> R_101;               dR_ZEB * a_ZEB_1 * C_ZEB_1
JCZ1DI:     C_ZEB_1 -> R_ZEB;               d_101 * b_101_ZEB_1 * C_ZEB_1
JCZ12:      C_ZEB_1 + R_200 -> C_ZEB_12;    2 * kOn_ZEB_200 * C_ZEB_1 * R_200 - kOff * C_ZEB_12
JCZ21:      C_ZEB_2 + R_101 -> C_ZEB_12;    kOn_ZEB_101 * C_ZEB_2 * R_101 - kOff * C_ZEB_12
JCZ12DR:    C_ZEB_12 -> R_101 + R_200;      dR_ZEB * a_ZEB_12 * C_ZEB_12
JCZ12DI1:   C_ZEB_12 -> C_ZEB_2;            d_101 * b_101_ZEB_12 * C_ZEB_12
JCZ12DI2:   C_ZEB_12 -> C_ZEB_1;            d_200 * b_200_ZEB_12 * C_ZEB_12
JCZ122:     C_ZEB_12 + R_200 -> C_ZEB_122;  kOn_ZEB_200 * C_ZEB_12 * R_200 - 2 * kOff * C_ZEB_122
JCZ221:     C_ZEB_22 + R_101 -> C_ZEB_122;  kOn_ZEB_101 * C_ZEB_22 * R_101 - kOff * C_ZEB_122
JCZ122DR:   C_ZEB_122 -> R_101 + 2 R_200;   dR_ZEB * a_ZEB_122 * C_ZEB_122
JCZ122DI1:  C_ZEB_122 -> C_ZEB_22;          d_101 * b_101_ZEB_122 * C_ZEB_122
JCZ122DI2:  C_ZEB_122 -> C_ZEB_12;          2 * d_200 * b_200_ZEB_122 * C_ZEB_122

JCSl:   R_SLUG + R_200 -> C1_SLUG;  kOn_SLUG * R_SLUG * R_200 - kOff * C1_SLUG
JCSlDR: C1_SLUG -> R_200;           dR_SLUG * a_SLUG * C1_SLUG
JCSlDI: C1_SLUG -> R_SLUG;          d_200 * b_SLUG * C1_SLUG

k_101 = 10; k_200 = 10; k_ZEB = 0.01; k_SNAI = 0.1; k_SLUG = 1
d_101 = 1; d_200 = 1; dR_ZEB = 1; dR_SNAI = 1; dR_SLUG = 1
l_ZEB = 1; l_SNAI = 1.8; l_SLUG = 5
dP_ZEB = 1; dP_SNAI = 1.25; dP_SLUG = 1.1

K_101_SNAI = 1.9; K_101_SLUG = 2.2
n_101_SNAI = 2; n_101_SLUG = 1
fc_101_SNAI = 0.1; fc_101_SLUG = 0.4

K_200_ZEB = 2.2; K_200_SNAI = 1.9; K_200_SLUG = 2.2
n_200_ZEB = 3; n_200_SNAI = 2; n_200_SLUG = 1
fc_200_ZEB = 0.1; fc_200_SNAI = 0.1; fc_200_SLUG = 0.4

K_ZEB_ZEB = 0.25; K_ZEB_SNAI = 1.8; K_ZEB_SLUG = 2
n_ZEB_ZEB = 2; n_ZEB_SNAI = 2; n_ZEB_SLUG = 2
fc_ZEB_ZEB = 7.5; fc_ZEB_SNAI = 10; fc_ZEB_SLUG = 4

K_SNAI_SNAI = 1.8; K_SNAI_SLUG = 2.25
n_SNAI_SNAI = 5; n_SNAI_SLUG = 3
fc_SNAI_SNAI = 0.4; fc_SNAI_SLUG = 0.5

K_SLUG_SLUG = 2.5; K_SLUG_SNAI = 1.8
n_SLUG_SLUG = 4; n_SLUG_SNAI = 1
fc_SLUG_SLUG = 4; fc_SLUG_SNAI = 0.5

t_ZEB_2 = 0.1; t_ZEB_22 = 0.1; t_ZEB_1 = 0.1; t_ZEB_12 = 0.1; t_ZEB_122 = 0.01; t_SLUG = 0.1
a_ZEB_2 = 2; a_ZEB_22 = 4; a_ZEB_1 = 2; a_ZEB_12 = 4; a_ZEB_122 = 8; a_SLUG = 2
b_200_ZEB_2 = 0.5; b_200_ZEB_22 = 0.25; b_200_ZEB_12 = 0.5; b_200_ZEB_122 = 0.25; b_101_ZEB_1 = 0.5; b_101_ZEB_12 = 0.5; b_101_ZEB_122 = 0.5; b_SLUG = 0.5
kOff = 100
KA_ZEB_101 = 10; KA_ZEB_200 = 10; KA_SLUG = 10
kOn_ZEB_101 := kOff * KA_ZEB_101; kOn_ZEB_200 := kOff * KA_ZEB_200; kOn_SLUG := kOff * KA_SLUG
''')

rna_ini_combs = np.power(10, sobol.sample(5, 150, skip=150) * 4.5 - 3)
ini_combs = np.zeros((rna_ini_combs.shape[0], len(runner.fs())))
next_rna_col = 0
rna_cols = dict()
protein_cols = dict()
for i, species in enumerate(runner.fs()):
    product = species.split('_')[-1]
    if species.startswith('R_'):
        ini_combs[:, i] = rna_ini_combs[:, next_rna_col]
        next_rna_col += 1
        rna_cols[product] = i
    elif species.startswith('X_'):
        protein_cols[product] = i

def update_protein_inicombs():
    for product, col in protein_cols.items():
        ini_combs[:, col] = runner[f'l_{product}'] * ini_combs[:, rna_cols[product]]

results = []
while len(results) < args.count:
    pset = {}
    for p in runner.ps():
        if p.startswith('k_'):
            pset[p] = 10 ** np.random.uniform(-1, 0.5 if p[2:].isdigit() else 1)
        elif p.startswith('l_'):
            pset[p] = 10 ** np.random.uniform(-0.2, 0.5)
        elif p.startswith('K_'):
            pset[p] = 10 ** np.random.uniform(-1.5, 2)
        elif p.startswith('a_'):
            pset[p] = 2 ** np.random.uniform(0, 3)
        elif p.startswith('b_'):
            pset[p] = 2 ** np.random.uniform(-3, 0)
        elif p.startswith('t_'):
            pset[p] = np.random.uniform(0, 1)
        elif p.startswith('KA_'):
            pset[p] = 10 ** np.random.uniform(2, 6)
    pset['a_ZEB_22'] = pset['a_ZEB_2'] * np.random.uniform(1, 4)
    pset['a_ZEB_12'] = max((pset['a_ZEB_1'] * pset['a_ZEB_2']) ** np.random.uniform(0.5, 1.5), pset['a_ZEB_1'], pset['a_ZEB_2'])
    pset['a_ZEB_122'] = max((pset['a_ZEB_1'] * pset['a_ZEB_22']) ** np.random.uniform(0.5, 1.5), pset['a_ZEB_1'], pset['a_ZEB_12'])
    pset['b_101_ZEB_12'] = pset['b_101_ZEB_1'] / np.random.uniform(1, 3)
    pset['b_101_ZEB_122'] = pset['b_101_ZEB_12'] * (2 ** np.random.uniform(-1, 0.5))
    pset['b_200_ZEB_22'] = pset['b_200_ZEB_2'] * (2 ** np.random.uniform(-1.5, 0.5))
    pset['b_200_ZEB_12'] = pset['b_200_ZEB_2'] / np.random.uniform(1, 3)
    pset['b_200_ZEB_122'] = pset['b_200_ZEB_22'] / np.random.uniform(1, 3)
    pset['t_ZEB_22'] = pset['t_ZEB_2'] ** np.random.uniform(1, 2)
    pset['t_ZEB_12'] = min((pset['t_ZEB_1'] * pset['t_ZEB_2']) ** np.random.uniform(0.5, 1.5), pset['t_ZEB_1'], pset['t_ZEB_2'])
    pset['t_ZEB_122'] = min((pset['t_ZEB_1'] * pset['t_ZEB_22']) ** np.random.uniform(0.5, 1.5), pset['t_ZEB_12'], pset['t_ZEB_22'])
    for p, value in pset.items():
        runner[p] = value
    update_protein_inicombs()
    attractors = multiattractor.findpointattractors(runner, ini_combs)
    if len(attractors) > 3:
        pset['attractors'] = len(attractors)
        results.append(pset)
        # print(attractors)
        # print(pset)
        if len(results) % 10 == 0 or len(results) == args.count:
            pd.DataFrame.from_dict(results).to_csv(args.output, index=False)
