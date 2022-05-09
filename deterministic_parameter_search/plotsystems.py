import argparse
import matplotlib.pyplot as plt
import matplotlib.ticker as mpltick
import multiattractor
import numpy as np
import pandas as pd
import sobol
import tellurium as te
import warnings

parser = argparse.ArgumentParser()
parser.add_argument('model', type=str, help='model SB')
parser.add_argument('psets', type=str, help='input CSV')
parser.add_argument('output', type=str, help='output image')
parser.add_argument('yspecies', type=str, help='species name for Y axis')
parser.add_argument('xspecies', nargs='+', type=str, help='species names for X axis')
parser.add_argument('--base', '-b', type=str, help='pset base file')
parser.add_argument('--init-species', '-i', nargs='+', help='species names with nonzero initial conditions')
parser.add_argument('--min-attractors', '-a', type=int, help='minimum number of attractors to plot')
parser.add_argument('--max-attractors', '-A', type=int, help='maximum number of attractors to plot')
parser.add_argument('--default-attractors', type=int, help='attractor count to assume if missing')
parser.add_argument('--skip', '-s', type=int, help='number of systems to skip over')
parser.add_argument('--systems', '-n', type=int, default=10, help='maximum number of systems to show')
parser.add_argument('--ic-base', '-x', type=int, default=5, help='exponent base for initial conditions')
parser.add_argument('--verbose', '-v', nargs='?', const='df', help='print psets and attractors (df, dict, or jl)')
parser.add_argument('--publish', action='store_true', help='set font size for publication figure')
parser.add_argument('--size', nargs=2, type=float, help='figure width and height (inches)')
parser.add_argument('--names', nargs='+', help='override species names: yname xname...')
parser.add_argument('--sharex', action='store_true', help='link X axis ranges')
parser.add_argument('--ylinear', action='store_true', help='use a linear scale for the Y axis')
parser.add_argument('--xrange', action='append', nargs=3, type=float, help='set an X log-range: index min max')
args = parser.parse_args()

with open(args.model) as f:
    runner = te.loada(f.read())

if args.init_species:
    nonzero_ini_combs = np.power(args.ic_base, sobol.sample(len(args.init_species), 100, skip=99) * 6 - 4.5)
    ini_combs = np.zeros((nonzero_ini_combs.shape[0], len(runner.fs())))
    next_rna_col = 0
    for i, species in enumerate(runner.fs()):
        if species in args.init_species:
            ini_combs[:, i] = nonzero_ini_combs[:, next_rna_col]
            next_rna_col += 1
    update_protein_inicombs = lambda: None
else:
    ini_combs, update_protein_inicombs = multiattractor.makeinitialconditions(runner, base=args.ic_base)

if args.publish:
    plt.rcParams['font.size'] = 11
    labelprops = {'fontsize': 14}
else:
    labelprops = {}

if (',' in args.yspecies or ',' in args.xspecies[0]) and len(args.xspecies) == 1:
    fig, ax = plt.subplots(figsize=((4, 4) if args.size is None else tuple(args.size)))
    axs = [ax]
    genesets = {}
    for axis, arg in zip(['x', 'y'], [args.xspecies[0], args.yspecies]):
        label, genes = arg.split('=')
        geneset = {'label': label, 'genes': genes.split(',')}
        genesets[axis] = geneset
    xlabels = [genesets['x']['label']]
    ylabel = genesets['y']['label']
else:
    if args.size:
        figprops = {'figsize': tuple(args.size), 'sharey': True}
    else:
        figprops = {}
    if args.sharex:
        figprops['sharex'] = True
    fig, axs = plt.subplots(ncols=len(args.xspecies), **figprops)
    if len(args.xspecies) == 1:
        axs = [axs]
    xlabels = args.xspecies
    ylabel = args.yspecies
    genesets = None

psets = pd.read_csv(args.psets).to_dict('records')

if args.base:
    with open(args.base) as f:
        pset_base = eval(f.read())
else:
    pset_base = {}

pindex = 0
shown = 0
for pset_extension in psets:
    pset = {**pset_base, **pset_extension}
    if args.systems and shown >= args.systems:
        break
    pset_attractors = pset['attractors'] if 'attractors' in pset else args.default_attractors
    if (args.min_attractors and pset_attractors < args.min_attractors) or (args.max_attractors and pset_attractors > args.max_attractors):
        continue
    if args.skip and pindex < args.skip:
        pindex += 1
        continue
    for p, value in pset.items():
        if p in runner.ps():
            runner[p] = value
    update_protein_inicombs()
    attractors = multiattractor.findpointattractors(runner, ini_combs)
    if len(attractors) != pset_attractors:
        warnings.warn(f'pset {pindex} claimed {pset_attractors} attractors but has {len(attractors)} now')
        pindex += 1
        continue
    df = pd.DataFrame.from_records(attractors, columns=runner.fs())
    if genesets:
        def score(genes):
            total = None
            for g in genes:
                log = np.log(df[g])
                genescore = (log - log.min()) / (log.max() - log.min())
                if total is None:
                    total = genescore
                else:
                    total += genescore
            return total
        df[genesets['x']['label']] = score(genesets['x']['genes'])
        df[genesets['y']['label']] = score(genesets['y']['genes'])
        df = df.sort_values(genesets['y']['label'])
        axs[0].plot(df[genesets['x']['label']], df[genesets['y']['label']], 'o-')
    else:
        df = df.sort_values(args.yspecies)
        for ax, species in zip(axs, args.xspecies):
            ax.plot(df[species], df[args.yspecies], 'o-')
    if args.verbose:
        print(f'--- pset {pindex} (color {shown}) ---')
        if args.verbose == 'df':
            print(pset)
            print(df)
        elif args.verbose == 'dict':
            print(pset)
            print(df.to_dict('records'))
        else:
            print('Dict(' + ', '.join(f'"{k}" => {v}' for (k, v) in pset.items() if k in runner.ps()) + ')')
            print(attractors)
        print()
    pindex += 1
    shown += 1

if args.verbose == 'jl':
    print(runner.fs())

axs[0].set_ylabel(ylabel if args.names is None else args.names[0], **labelprops)
for ax, species in zip(axs, xlabels if args.names is None else args.names[1:]):
    ax.set_xlabel(species, **labelprops)
if genesets:
    axs[0].set_box_aspect(1)
else:
    for ax in axs:
        ax.set_xscale('log')
        if not args.ylinear:
            ax.set_yscale('log')
        if args.publish:
            ax.xaxis.set_minor_locator(mpltick.LogLocator(base=10, numticks=10))
            ax.xaxis.set_minor_formatter(mpltick.NullFormatter())
            if args.ylinear:
                ax.yaxis.set_major_locator(mpltick.MultipleLocator(0.5))
                ax.yaxis.set_minor_locator(mpltick.MultipleLocator(0.1))
            else:
                ax.yaxis.set_minor_locator(mpltick.LogLocator(base=10, numticks=10))
                ax.yaxis.set_minor_formatter(mpltick.NullFormatter())
    if args.xrange:
        for i, minlog, maxlog in args.xrange:
            axs[int(i)].set_xlim(left=(10 ** minlog), right=(10 ** maxlog))
fig.tight_layout()
fig.savefig(args.output)