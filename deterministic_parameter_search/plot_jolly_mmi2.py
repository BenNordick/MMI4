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
parser.add_argument('--attractors', type=int, default=5, help='attractor count')
parser.add_argument('--limit', type=int, default=10, help='number of systems to check')
parser.add_argument('--start', type=int, default=0, help='first system to show')
args = parser.parse_args()

runner = te.loada('''
function Hminus(B, B0, n)
    1 / (1 + (B / B0)^n)
end
function Hs(B, B0, n, fc)
    Hminus(B, B0, n) + fc * (1 - Hminus(B, B0, n))
end

JI: -> R_200; k_200 * Hs(X_ZEB, K_200_ZEB, n_200_ZEB, fc_200_ZEB) * Hs(X_SNAI, K_200_SNAI, n_200_SNAI, fc_200_SNAI) * Hs(X_SLUG, K_200_SLUG, n_200_SLUG, fc_200_SLUG) - d_200 * R_200
JRZ: -> R_ZEB; k_ZEB * Hs(X_ZEB, K_ZEB_ZEB, n_ZEB_ZEB, fc_ZEB_ZEB) * Hs(X_SNAI, K_ZEB_SNAI, n_ZEB_SNAI, fc_ZEB_SNAI) * Hs(X_SLUG, K_ZEB_SLUG, n_ZEB_SLUG, fc_ZEB_SLUG) - dR_ZEB * R_ZEB
JRSn: -> R_SNAI; k_SNAI * Hs(X_SNAI, K_SNAI_SNAI, n_SNAI_SNAI, fc_SNAI_SNAI) * Hs(X_SLUG, K_SNAI_SLUG, n_SNAI_SLUG, fc_SNAI_SLUG) - dR_SNAI * R_SNAI
JRSl: -> R_SLUG; k_SLUG * Hs(X_SLUG, K_SLUG_SLUG, n_SLUG_SLUG, fc_SLUG_SLUG) * Hs(X_SNAI, K_SLUG_SNAI, n_SLUG_SNAI, fc_SLUG_SNAI) - dR_SLUG * R_SLUG
JRC: -> R_CDH; k_CDH * Hs(X_SLUG, K_CDH_SLUG, n_CDH_SLUG, fc_CDH_SLUG) * Hs(X_SNAI, K_CDH_SNAI, n_CDH_SNAI, fc_CDH_SNAI) * Hs(X_ZEB, K_CDH_ZEB, n_CDH_ZEB, fc_CDH_ZEB) - dR_CDH * R_CDH

JPZ: -> X_ZEB; l_ZEB * (R_ZEB + t_ZEB_1 * C1_ZEB + t_ZEB_2 * C2_ZEB) - dP_ZEB * X_ZEB
JPSn: -> X_SNAI; l_SNAI * R_SNAI - dP_SNAI * X_SNAI
JPSl: -> X_SLUG; l_SLUG * (R_SLUG + t_SLUG * C1_SLUG) - dP_SLUG * X_SLUG
JPC: -> X_CDH; l_CDH * R_CDH - dP_CDH * X_CDH

JC1Z: R_ZEB + R_200 -> C1_ZEB; 2 * kOn_ZEB * R_ZEB * R_200 - kOff_ZEB * C1_ZEB
JC1ZDR: C1_ZEB -> R_200; dR_ZEB * a_ZEB_1 * C1_ZEB
JC1ZDI: C1_ZEB -> R_ZEB; d_200 * b_ZEB_1 * C1_ZEB
JC2Z: C1_ZEB + R_200 -> C2_ZEB; kOn_ZEB * C1_ZEB * R_200 - 2 * kOff_ZEB * C2_ZEB
JC2ZDR: C2_ZEB -> 2 R_200; dR_ZEB * a_ZEB_2 * C2_ZEB
JC2ZDI: C2_ZEB -> C1_ZEB; 2 * d_200 * b_ZEB_2 * C2_ZEB

JCSl: R_SLUG + R_200 -> C1_SLUG; kOn_SLUG * R_SLUG * R_200 - kOff_SLUG * C1_SLUG
JCSlDR: C1_SLUG -> R_200; dR_SLUG * a_SLUG * C1_SLUG
JCSlDI: C1_SLUG -> R_SLUG; d_200 * b_SLUG * C1_SLUG

k_200 = 10; k_ZEB = 0.01; k_SNAI = 0.1; k_SLUG = 1; k_CDH = 10
d_200 = 1; dR_ZEB = 1; dR_SNAI = 1; dR_SLUG = 1; dR_CDH = 1
l_ZEB = 1; l_SNAI = 1.8; l_SLUG = 5; l_CDH = 1
dP_ZEB = 1; dP_SNAI = 1.25; dP_SLUG = 1.1; dP_CDH = 1

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

K_CDH_SLUG = 0.5; K_CDH_SNAI = 0.5; K_CDH_ZEB = 1
n_CDH_SLUG = 2; n_CDH_SNAI = 2; n_CDH_ZEB = 2
fc_CDH_SLUG = 0.3; fc_CDH_SNAI = 0.3; fc_CDH_ZEB = 0.75

t_ZEB_1 = 0.1; t_ZEB_2 = 0.1; t_SLUG = 0.1
a_ZEB_1 = 2; a_ZEB_2 = 4; a_SLUG = 2
b_ZEB_1 = 0.5; b_ZEB_2 = 0.25; b_SLUG = 0.5
kOff_ZEB = 100; kOff_SLUG = 100
KA_ZEB = 10; KA_SLUG = 10
kOn_ZEB := kOff_ZEB * KA_ZEB; kOn_SLUG := kOff_SLUG * KA_SLUG
''')

rna_ini_combs = np.power(5, sobol.sample(5, 100, skip=99) * 6 - 4.5)
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

base = {'k_200': 0.9496069628526422, 'K_200_ZEB': 0.127345680065309, 'K_200_SNAI': 0.2630962417982151, 'K_200_SLUG': 13.034166192151783, 
        'k_ZEB': 0.7271460964706331, 'K_ZEB_ZEB': 7.792059960170112, 'K_ZEB_SNAI': 0.8060389305197769, 'K_ZEB_SLUG': 0.0457781343265527, 
        'k_SNAI': 0.3970372010577482, 'K_SNAI_SNAI': 0.1778751274114133, 'K_SNAI_SLUG': 10.745375581449895, 'k_SLUG': 0.5687185521885274, 
        'K_SLUG_SLUG': 0.5704651975937381, 'K_SLUG_SNAI': 3.789913067473965, 'l_ZEB': 0.7366281866507639, 't_ZEB_1': 0.3006972304193904, 
        'l_SNAI': 0.331467636646333, 'l_SLUG': 0.7457170892176783, 't_SLUG': 0.1462260067498192, 'a_ZEB_1': 5.04830442566583, 'b_ZEB_1': 0.1302556574684447, 
        'a_SLUG': 8.889856943150424, 'b_SLUG': 0.0928846366194664, 'KA_ZEB': 36.74723689940514, 'KA_SLUG': 4.194993534849837}
for p, value in base.items():
    runner[p] = value
update_protein_inicombs()

psets = pd.read_csv(args.input).to_dict('records')
# One parameter set originally found by extending a tetrastable parameterization of a model very similar to the Jolly lab's
psets.insert(0, {'k_200': 0.5839497832890144, 'k_ZEB': 0.5820778984946993, 'k_SNAI': 0.66542133902077, 'k_SLUG': 0.3826439137691515, 'KA_ZEB': 2056.3180088009767, 'KA_SLUG': 178.40238711823713, 'a_ZEB_2': 12.430561391344789, 'b_ZEB_2': 0.1007638327443842, 't_ZEB_2': 0.0668216778075153, 'attractors': 5})

plt.rcParams['font.size'] = 11
fig, axs = plt.subplots(ncols=3, figsize=(6, 3), sharey=True)

pset_index = 0
for pset in psets:
    if pset['attractors'] != args.attractors:
        continue
    if pset_index < args.start:
        pset_index += 1
        continue
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
    axs[2].plot(df['R_200'], df['X_ZEB'], 'o-')
    
    print(f'--- pset {pset_index} ---')
    print(pset)
    print(attractors)
    print()
    pset_index += 1
    if pset_index >= args.limit:
        break
print(runner.fs())

axs[0].set_ylabel('Zeb1 protein (AU)', fontsize=14)
axs[0].set_xlabel('E-cadherin (AU)', fontsize=14)
axs[0].set_xscale('log', base=2)
axs[1].set_xlabel('Slug (AU)', fontsize=14)
axs[1].set_xscale('log')
axs[1].xaxis.set_minor_locator(mpltick.LogLocator(base=10, numticks=10))
axs[1].xaxis.set_minor_formatter(mpltick.NullFormatter())
axs[2].set_xlabel('miR-200 (AU)', fontsize=14)
axs[2].set_xscale('log')
axs[2].xaxis.set_minor_locator(mpltick.LogLocator(base=10, numticks=10))
axs[2].xaxis.set_minor_formatter(mpltick.NullFormatter())
for ax in axs:
    ax.set_yscale('log')
    ax.yaxis.set_minor_locator(mpltick.LogLocator(base=10, numticks=10))
    ax.yaxis.set_minor_formatter(mpltick.NullFormatter())
fig.tight_layout()
fig.savefig(args.output)
