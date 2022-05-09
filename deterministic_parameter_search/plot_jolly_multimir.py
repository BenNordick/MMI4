import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as mpltick
import multiattractor
import numpy as np
import pandas as pd
import sobol
import tellurium as te

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='input CSV')
parser.add_argument('output', type=str, help='output image')
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
JRC:    -> R_CDH;   k_CDH * Hs(X_SLUG, K_CDH_SLUG, n_CDH_SLUG, fc_CDH_SLUG) * Hs(X_SNAI, K_CDH_SNAI, n_CDH_SNAI, fc_CDH_SNAI) * Hs(X_ZEB, K_CDH_ZEB, n_CDH_ZEB, fc_CDH_ZEB) - dR_CDH * R_CDH

JPZ:    -> X_ZEB;   l_ZEB * (R_ZEB + t_ZEB_2 * C_ZEB_2 + t_ZEB_22 * C_ZEB_22 + t_ZEB_1 * C_ZEB_1 + t_ZEB_12 * C_ZEB_12 + t_ZEB_122 * C_ZEB_122) - dP_ZEB * X_ZEB
JPSn:   -> X_SNAI;  l_SNAI * R_SNAI - dP_SNAI * X_SNAI
JPSl:   -> X_SLUG;  l_SLUG * (R_SLUG + t_SLUG * C1_SLUG) - dP_SLUG * X_SLUG
JPC:    -> X_CDH;   l_CDH * R_CDH - dP_CDH * X_CDH

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
JCZ122:     C_ZEB_12 + R_200 -> C_ZEB_122;  kOn_ZEB_200 * C_ZEB_12 * R_200 - kOff * C_ZEB_122
JCZ221:     C_ZEB_22 + R_101 -> C_ZEB_122;  kOn_ZEB_101 * C_ZEB_22 * R_101 - kOff * C_ZEB_122
JCZ122DR:   C_ZEB_122 -> R_101 + 2 R_200;   dR_ZEB * a_ZEB_122 * C_ZEB_122
JCZ122DI1:  C_ZEB_122 -> C_ZEB_22;          d_101 * b_101_ZEB_122 * C_ZEB_122
JCZ122DI2:  C_ZEB_122 -> C_ZEB_12;          2 * d_200 * b_200_ZEB_122 * C_ZEB_122

JCSl:   R_SLUG + R_200 -> C1_SLUG;  kOn_SLUG * R_SLUG * R_200 - kOff * C1_SLUG
JCSlDR: C1_SLUG -> R_200;           dR_SLUG * a_SLUG * C1_SLUG
JCSlDI: C1_SLUG -> R_SLUG;          d_200 * b_SLUG * C1_SLUG

k_101 = 10; k_200 = 10; k_ZEB = 0.01; k_SNAI = 0.1; k_SLUG = 1; k_CDH = 10
d_101 = 1; d_200 = 1; dR_ZEB = 1; dR_SNAI = 1; dR_SLUG = 1; dR_CDH = 1
l_ZEB = 1; l_SNAI = 1.8; l_SLUG = 5; l_CDH = 1
dP_ZEB = 1; dP_SNAI = 1.25; dP_SLUG = 1.1; dP_CDH = 1

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

// K_CDH_SLUG = 2; K_CDH_SNAI = 2; K_CDH_ZEB = 2
K_CDH_SLUG = 0.5; K_CDH_SNAI = 0.5; K_CDH_ZEB = 1
n_CDH_SLUG = 2; n_CDH_SNAI = 2; n_CDH_ZEB = 2
fc_CDH_SLUG = 0.3; fc_CDH_SNAI = 0.3; fc_CDH_ZEB = 0.75

t_ZEB_2 = 0.1; t_ZEB_22 = 0.1; t_ZEB_1 = 0.1; t_ZEB_12 = 0.1; t_ZEB_122 = 0.01; t_SLUG = 0.1
a_ZEB_2 = 2; a_ZEB_22 = 4; a_ZEB_1 = 2; a_ZEB_12 = 4; a_ZEB_122 = 8; a_SLUG = 2
b_200_ZEB_2 = 0.5; b_200_ZEB_22 = 0.25; b_200_ZEB_12 = 0.5; b_200_ZEB_122 = 0.25; b_101_ZEB_1 = 0.5; b_101_ZEB_12 = 0.5; b_101_ZEB_122 = 0.5; b_SLUG = 0.5
kOff = 100
KA_ZEB_101 = 10; KA_ZEB_200 = 10; KA_SLUG = 10
kOn_ZEB_101 := kOff * KA_ZEB_101; kOn_ZEB_200 := kOff * KA_ZEB_200; kOn_SLUG := kOff * KA_SLUG
''')

rna_ini_combs = np.power(10, sobol.sample(6, 150, skip=150) * 4.5 - 3)
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

psets = pd.read_csv(args.input).to_dict('records')

plt.rcParams['font.size'] = 11
fig, axs = plt.subplots(ncols=4, figsize=(6, 3), sharey=True)

pindex = 0
for pset in psets:
    for p, value in pset.items():
        if p in runner.ps():
            runner[p] = value
    update_protein_inicombs()
    attractors = multiattractor.findpointattractors(runner, ini_combs)
    df = pd.DataFrame.from_records(attractors, columns=runner.fs())
    runner['K_CDH_SLUG'] = df['X_SLUG'].median()
    runner['K_CDH_SNAI'] = df['X_SNAI'].median()
    runner['K_CDH_ZEB'] = df['X_ZEB'].median() / 2
    attractors = multiattractor.findpointattractors(runner, ini_combs)
    df = pd.DataFrame.from_records(attractors, columns=runner.fs()).sort_values('X_ZEB')
    axs[0].plot(df['X_CDH'], df['X_ZEB'], 'o-')
    axs[1].plot(df['X_SLUG'], df['X_ZEB'], 'o-')
    axs[2].plot(df['R_101'], df['X_ZEB'], 'o-')
    axs[3].plot(df['R_200'], df['X_ZEB'], 'o-')
    print(f'--- pset {pindex} ---')
    print(', '.join(f'"{p}" => {runner[p]}' for p in runner.ps()))
    print(runner.fs())
    print(attractors)
    print()
    pindex += 1

axs[0].set_ylabel('Zeb1 protein (AU)', fontsize=14)
axs[0].set_xlabel('E-cadherin', fontsize=14)
axs[1].set_xlabel('Slug', fontsize=14)
axs[2].set_xlabel('miR-101', fontsize=14)
axs[3].set_xlabel('miR-200', fontsize=14)
for ax in axs:
    ax.set_xscale('log', base=(2 if ax is axs[0] else 10))
    ax.set_yscale('log')
    ax.yaxis.set_minor_locator(mpltick.LogLocator(base=10, numticks=10))
    ax.yaxis.set_minor_formatter(mpltick.NullFormatter())
    ax.xaxis.set_minor_locator(mpltick.LogLocator(base=10, numticks=10))
    ax.xaxis.set_minor_formatter(mpltick.NullFormatter())
axs[2].xaxis.set_major_locator(mpltick.FixedLocator([1e-6, 1e-1]))
axs[3].xaxis.set_major_locator(mpltick.FixedLocator([1e-7, 1e-2]))
fig.tight_layout()
fig.savefig(args.output)
