import argparse
import multiattractor
import numpy as np
import pandas as pd
import random
import tellurium as te

parser = argparse.ArgumentParser()
parser.add_argument('count', type=int, help='number of psets to find')
parser.add_argument('mmi', type=str, help='tetrastable systems CSV')
parser.add_argument('tma', type=str, help='bistable TMA systems CSV')
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

// SNAIL with 2 miR-34 binding sites, transcriptionally activated by TWIST1
JRS: -> R_SNAI; k_SNAI * (f_SNAI_cons + f_SNAI_TWIST * (X_TWIST / K_SNAI_TWIST)^n_SNAI_TWIST / (1 + (X_TWIST / K_SNAI_TWIST)^n_SNAI_TWIST)) - R_SNAI
JR34: -> R_34; k_34 - gm_34 * R_34
JC1S: R_SNAI + R_34 -> C1_SNAI; 2 * kOn_SNAI * R_SNAI * R_34 - kOff_SNAI * C1_SNAI
JC2S: C1_SNAI + R_34 -> C2_SNAI; kOn_SNAI * C1_SNAI * R_34 - 2 * kOff_SNAI * C2_SNAI
JC1SDR: C1_SNAI -> R_34; a1_SNAI * C1_SNAI
JC1SDI: C1_SNAI -> R_SNAI; b1_SNAI * gm_34 * C1_SNAI
JC2SDR: C2_SNAI -> 2 R_34; a2_SNAI * C2_SNAI
JC2SDI: C2_SNAI -> C1_SNAI; 2 * b2_SNAI * gm_34 * C2_SNAI

// TWIST1 transcriptionally activated by SNAIL
JRT: -> R_TWIST; k_TWIST * (f_TWIST_cons + f_TWIST_SNAI * (X_SNAI / K_TWIST_SNAI)^n_TWIST_SNAI / (1 + (X_SNAI / K_TWIST_SNAI)^n_TWIST_SNAI)) - R_TWIST

JPZ: -> X_ZEB; l_ZEB * R_ZEB - X_ZEB
JPS: -> X_SNAI; l_SNAI * R_SNAI - X_SNAI
JPT: -> X_TWIST; l_TWIST * R_TWIST - X_TWIST

R_ZEB = 1.0; R_200 = 1.0; C1_ZEB = 0.0; C2_ZEB = 0.0; X_ZEB = 0.0
R_SNAI = 1.0; R_34 = 1.0; C1_SNAI = 0.0; C2_SNAI = 0.0; X_SNAI = 0.0
R_TWIST = 1.0; X_TWIST = 0.0

// MMI
a1_ZEB = 5; b1_ZEB = 0.65; a2_ZEB = 10; b2_ZEB = 0.4
gm_200 = 1.0; kOff_ZEB = 100; K_ZEB = 5000; 
kOn_ZEB := K_ZEB * kOff_ZEB
a1_SNAI = 5; b1_SNAI = 0.65; a2_SNAI = 10; b2_SNAI = 0.4
gm_34 = 1.0; kOff_SNAI = 100; K_SNAI = 5000; 
kOn_SNAI := K_SNAI * kOff_SNAI

// Transcription
k_ZEB = 1; k_SNAI = 1; k_TWIST = 1
k_200 = 0.1; k_34 = 0.1
f_ZEB_SNAI = 0.5; f_ZEB_cons := 1 - f_ZEB_SNAI
K_ZEB_SNAI = 0.1; n_ZEB_SNAI = 1
f_SNAI_TWIST = 0.5; f_SNAI_cons := 1 - f_SNAI_TWIST
K_SNAI_TWIST = 0.1; n_SNAI_TWIST = 2
f_TWIST_SNAI = 0.5; f_TWIST_cons := 1 - f_TWIST_SNAI
K_TWIST_SNAI = 0.1; n_TWIST_SNAI = 2

// Translation
l_ZEB = 2; l_SNAI = 2; l_TWIST = 2
''')

ini_combs, update_protein_inicombs = multiattractor.makeinitialconditions(runner)
mmi = pd.read_csv(args.mmi).to_dict('records')
tma = pd.read_csv(args.tma).to_dict('records')

results = []
while len(results) < args.count:
    m = random.choice(range(len(mmi)))
    pset = {k: v for (k, v) in mmi[m].items() if k in runner.ps()}
    t = random.choice(range(len(tma)))
    p_tma = tma[t]
    pset['f_SNAI_TWIST'] = p_tma['r_A']
    pset['K_SNAI_TWIST'] = p_tma['K_A_B']
    pset['n_SNAI_TWIST'] = p_tma['n_A_B']
    pset['k_TWIST'] = p_tma['k_B']
    pset['f_TWIST_SNAI'] = p_tma['r_B']
    pset['K_TWIST_SNAI'] = p_tma['K_B_A']
    pset['n_TWIST_SNAI'] = p_tma['n_B_A']
    pset['l_TWIST'] = p_tma['l_B']
    pset['l_SNAI'] = np.random.uniform(pset['l_SNAI'], p_tma['l_A'])
    pset['k_SNAI'] = np.random.uniform(pset['k_SNAI'], 1)
    for (p, value) in pset.items():
        runner[p] = value
    update_protein_inicombs()
    attractors = multiattractor.findpointattractors(runner, ini_combs)
    if len(attractors) > 4:
        pset['attractors'] = len(attractors)
        # df = multiattractor.pointattractorsdf(runner, attractors).sort_values('X_ZEB')
        # print(pset)
        # print(df)
        # print()
        results.append(pset)
        pd.DataFrame.from_dict(results).to_csv(args.output, index=False)
        