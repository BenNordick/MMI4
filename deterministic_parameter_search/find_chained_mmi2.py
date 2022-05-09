import argparse
import multiattractor
import numpy as np
import pandas as pd
import tellurium as te

parser = argparse.ArgumentParser()
parser.add_argument('count', type=int, help='number of psets to find')
parser.add_argument('output', type=str, help='output CSV')
args = parser.parse_args()

runner = te.loada('''
// ZEB1 with 2 miR-200 binding sites, transcriptionally activated by SNAIL
JRZ: -> R_ZEB; k_ZEB * (f_ZEB_cons + f_ZEB_SNAI * (X_SNAI / K_ZEB_SNAI)^n_ZEB_SNAI / (1 + (X_SNAI / K_ZEB_SNAI)^n_ZEB_SNAI)) - R_ZEB
JR200: -> R_200; k_200 - gm_200 * R_200
JC1Z: R_ZEB + R_200 -> C1_ZEB; 2 * kOn_ZEB * R_ZEB * R_200 - kOff_ZEB * C1_ZEB
JC2Z: C1_ZEB + R_200 -> C2_ZEB; kOn_ZEB * C1_ZEB * R_200 - 2 * kOff_ZEB * C2_ZEB
JC1ZDR: C1_ZEB -> R_200; a1_ZEB * C1_ZEB
JC1ZDI: C1_ZEB -> R_ZEB; b1_ZEB * gm_200 * C1_ZEB
JC2ZDR: C2_ZEB -> 2 R_200; a2_ZEB * C2_ZEB
JC2ZDI: C2_ZEB -> C1_ZEB; 2 * b2_ZEB * gm_200 * C2_ZEB

// SNAIL with 2 miR-34 binding sites
JRS: -> R_SNAI; k_SNAI - R_SNAI
JR34: -> R_34; k_34 - gm_34 * R_34
JC1S: R_SNAI + R_34 -> C1_SNAI; 2 * kOn_SNAI * R_SNAI * R_34 - kOff_SNAI * C1_SNAI
JC2S: C1_SNAI + R_34 -> C2_SNAI; kOn_SNAI * C1_SNAI * R_34 - 2 * kOff_SNAI * C2_SNAI
JC1SDR: C1_SNAI -> R_34; a1_SNAI * C1_SNAI
JC1SDI: C1_SNAI -> R_SNAI; b1_SNAI * gm_34 * C1_SNAI
JC2SDR: C2_SNAI -> 2 R_34; a2_SNAI * C2_SNAI
JC2SDI: C2_SNAI -> C1_SNAI; 2 * b2_SNAI * gm_34 * C2_SNAI

JPZ: -> X_ZEB; l_ZEB * R_ZEB - X_ZEB
JPS: -> X_SNAI; l_SNAI * R_SNAI - X_SNAI

R_ZEB = 1.0; R_200 = 1.0; C1_ZEB = 0.0; C2_ZEB = 0.0; X_ZEB = 0.0
R_SNAI = 1.0; R_34 = 1.0; C1_SNAI = 0.0; C2_SNAI = 0.0; X_SNAI = 0.0

// MMI
a1_ZEB = 5; b1_ZEB = 0.65; a2_ZEB = 10; b2_ZEB = 0.4
gm_200 = 1.0; kOff_ZEB = 100; K_ZEB = 5000; 
kOn_ZEB := K_ZEB * kOff_ZEB
a1_SNAI = 5; b1_SNAI = 0.65; a2_SNAI = 10; b2_SNAI = 0.4
gm_34 = 1.0; kOff_SNAI = 100; K_SNAI = 5000; 
kOn_SNAI := K_SNAI * kOff_SNAI

// Transcription
k_ZEB = 1; k_SNAI = 1
k_200 = 0.1; k_34 = 0.1
f_ZEB_SNAI = 0.5; f_ZEB_cons := 1 - f_ZEB_SNAI
K_ZEB_SNAI = 0.1; n_ZEB_SNAI = 1

// Translation
l_ZEB = 2; l_SNAI = 2
''')

ini_combs, update_protein_inicombs = multiattractor.makeinitialconditions(runner)

results = []
while len(results) < args.count:
    pset = {}
    pset['f_ZEB_SNAI'] = np.random.uniform(0.4, 0.95)
    pset['K_ZEB_SNAI'] = 10 ** np.random.uniform(-2, 1)
    pset['k_ZEB'] = 10 ** np.random.uniform(0, 0.5)
    pset['k_SNAI'] = 10 ** np.random.uniform(-0.5, 0.5)
    for mi in ['200', '34']:
        pset[f'k_{mi}'] = 10 ** np.random.uniform(-1.5, -0.2)
    for target in ['ZEB', 'SNAI']:
        pset[f'l_{target}'] = 10 ** np.random.uniform(0, 1)
        pset[f'a1_{target}'] = 2 ** np.random.uniform(0, 3)
        pset[f'a2_{target}'] = min(pset[f'a1_{target}'] * 2 ** np.random.uniform(0, 2), 32)
        pset[f'b1_{target}'] = 2 ** np.random.uniform(-2, 0)
        pset[f'b2_{target}'] = pset[f'b1_{target}'] * 2 ** np.random.uniform(-2, 0)
    for (p, value) in pset.items():
        runner[p] = value
    update_protein_inicombs()
    attractors = multiattractor.findpointattractors(runner, ini_combs)
    if len(attractors) >= 4:
        pset['attractors'] = len(attractors)
        df = multiattractor.pointattractorsdf(runner, attractors).sort_values('X_ZEB')
        other_concs = list(df['X_SNAI'])
        pset['correlated'] = other_concs[3] > other_concs[1] and other_concs[2] > other_concs[1] and other_concs[3] > other_concs[0]
        # print(pset)
        # print(df)
        # print()
        results.append(pset)
        if len(results) % 5 == 0 or len(results) == args.count:
            pd.DataFrame.from_dict(results).to_csv(args.output, index=False)

