import bifurcation
import matplotlib.pyplot as plt
import matplotlib.ticker as mpltick
import multiattractor
import numpy as np
import sobol
import tellurium as te

runner = te.loada('''
JF: $F -> mRNA; 0

J0: -> mRNA; k_ZEB1 - 2 * kOn * mRNA * miRNA + kOff * C1 + b1 * gm * C1 - mRNA
//            -- TRANSCRIPTION CONTROL ------------------------------------   -- MMI -->
J1: -> miRNA; mu * ((1 - r_I) + r_I / (1 + (X_ZEB1 / K_I_ZEB1) ^ n_I_ZEB1)) - 2 * kOn * mRNA * miRNA + kOff * C1 + a1 * C1 - kOn * C1 * miRNA + 2 * kOff * C2 + 2 * a2 * C2 - gm * miRNA
J2: -> C1; 2 * kOn * mRNA * miRNA - kOff * C1 - a1 * C1 - b1 * gm * C1 - kOn * C1 * miRNA + 2 * kOff * C2 + 2 * b2 * gm * C2
J3: -> C2; kOn * C1 * miRNA - 2 * kOff * C2 - a2 * C2 - 2 * b2 * gm * C2
J4: -> X_ZEB1; l_ZEB1 * mRNA - X_ZEB1

mRNA = 1.0; miRNA = 1.0; C1 = 0.0; C2 = 0.0; X_ZEB1 = 0.0
// MMI
a1 = 5; b1 = 0.65; a2 = 10; b2 = 0.4
gm = 1.0; kOff = 100; K = 5000; k_ZEB1 = 1
kOn := K * kOff
// Repression
mu = 0.07; r_I = 0.93; K_I_ZEB1 = 0.9; n_I_ZEB1 = 5; l_ZEB1 = 2

$F = 0
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

# tscn_domains = {'K': (0.1, 2), 'r': [0.9, 0.99], 'n': [0.5, 6.5], 'l': [1, 10], 'k': (1, 1), 'mu': (1, 1)}
# while True:
#     for p in ['K_I_ZEB1', 'l_ZEB1', 'k_ZEB1', 'mu', 'r_I']:
#         domain = tscn_domains[p.split('_')[0]]
#         if isinstance(domain, list):
#             value = np.random.uniform(low=domain[0], high=domain[1])
#         else:
#             mean, sd = domain
#             value = np.random.lognormal(np.log(mean) - (sd ** 2) / 2, sd)
#         if p[0] == 'n':
#             value = np.round(value)
#         runner[p] = value
#         if p[0] == 'l':
#             ini_combs[:, protein_col] = ini_combs[:, mrna_col] * value
#     attractors = multiattractor.findpointattractors(runner, ini_combs)
#     if len(attractors) == 3:
#         print(attractors)
#         print({p: runner[p] for p in runner.ps()})
#         input()

attractors = multiattractor.findpointattractors(runner, ini_combs)
print(multiattractor.pointattractorsdf(runner, attractors))

### MRNA TRANSCRIPTION RATE

bi_diagram = bifurcation.run_bifurcation(runner, 'k_ZEB1', (0.001, 20), ScanDirection='Negative')
regions = bifurcation.region_attractor_counts(bi_diagram)

plt.rcParams['font.size'] = 11
fig, axs = plt.subplots(ncols=2, sharex=True, figsize=(5, 3))
bifurcation.plot_split_bifurcation_diagram(runner, bi_diagram, ['X_ZEB1', 'miRNA'], axs, ['tab:blue', 'tab:green', 'tab:orange'])
for ax in axs:
    ax.set_yscale('log')
    ax.set_ylim(top=5, bottom=1e-6)
    ax.set_xlim(left=0.5, right=1)
    ax.set_xlabel('$k_R$', fontsize=14)
    ax.xaxis.set_minor_locator(mpltick.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(mpltick.LogLocator(base=10, numticks=10))
    ax.yaxis.set_minor_formatter(mpltick.NullFormatter())
axs[0].set_ylabel('Protein (AU)')
axs[1].set_ylabel('microRNA (AU)')
fig.tight_layout()
fig.savefig('images/bif_mmi2_mirep.svg')
print(regions)

### MICRORNA DEGRADATION RATE

### DEGRADATION RATE

runner['k_ZEB1'] = 0.9

bi_diagram = bifurcation.run_bifurcation(runner, 'gm', (0.001, 20), ScanDirection='Negative')
regions = bifurcation.region_attractor_counts(bi_diagram)

fig, axs = plt.subplots(ncols=2, sharex=True, figsize=(5, 3))
bifurcation.plot_split_bifurcation_diagram(runner, bi_diagram, ['X_ZEB1', 'miRNA'], axs, ['tab:blue', 'tab:green', 'tab:orange'])
for ax in axs:
    ax.set_yscale('log')
    ax.set_ylim(top=5, bottom=1e-6)
    ax.set_xlim(left=0.5, right=1.5)
    ax.set_xlabel(R'$\gamma$', fontsize=14)
    ax.xaxis.set_minor_locator(mpltick.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(mpltick.LogLocator(base=10, numticks=10))
    ax.yaxis.set_minor_formatter(mpltick.NullFormatter())
axs[0].set_ylabel('Protein (AU)')
axs[1].set_ylabel('microRNA (AU)')
fig.tight_layout()
fig.savefig('images/bif_mmi2_mirep_deg.svg')
print(regions)
