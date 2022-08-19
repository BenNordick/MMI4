import argparse
import collections
import matplotlib.gridspec as mplgrid
import matplotlib.patches as mplpatch
import matplotlib.pyplot as plt
import matplotlib.ticker as mpltick
import numpy as np
import pandas as pd
import re

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='input sensitivity CSV file')
parser.add_argument('paramdef', type=str, help='parameter names/order file')
parser.add_argument('output', type=str, help='output image file')
parser.add_argument('--figsize', nargs=2, type=float, help='figure size in inches')
parser.add_argument('--scatter', nargs='?', type=float, const=1, help='plot values as points')
parser.add_argument('--median', action='store_true', help='mark medians and extrema')
parser.add_argument('--ignore', type=str, help='regex for parameter names to ignore')
args = parser.parse_args()

param_data = collections.defaultdict(list)
for row in pd.read_csv(args.input).to_dict('records'):
    for col in row:
        if not col.startswith('original_'):
            continue
        param = col[len('original_'):]
        if args.ignore and re.match(args.ignore, param):
            continue
        if not row.get(f'validated_{param}', True):
            continue
        param_data[param].append({'original': row[col], 'min': row[f'min_{param}'], 'max': row[f'max_{param}']})

param_names = {}
with open(args.paramdef) as f:
    for line in f:
        if line.strip() == '':
            continue
        p, name = line.split('=', maxsplit=1)
        param_names[p.strip()] = '$' + name.strip() + '$'
for p in param_data:
    if p not in param_names:
        param_names[p] = p

def uses_raw_span(param: str):
    return param.startswith(('f_', 'r_', 't_'))

main_params = {k: v for k, v in param_data.items() if not uses_raw_span(k)}
rawspan_params = {k: v for k, v in param_data.items() if uses_raw_span(k)}

width_ratios = [1, 2]
if rawspan_params:
    fig = plt.figure(figsize=args.figsize)
    gs = mplgrid.GridSpec(nrows=2, ncols=2, wspace=0, width_ratios=width_ratios, height_ratios=[len(main_params), len(rawspan_params)])
    ax_left = fig.add_subplot(gs[0, 0])
    ax_right = fig.add_subplot(gs[0, 1], sharey=ax_left)
    ax_rawspan = fig.add_subplot(gs[1, :])
else:
    fig, (ax_left, ax_right) = plt.subplots(1, 2, sharey=True, figsize=args.figsize, gridspec_kw={'width_ratios': width_ratios, 'wspace': 0})
    ax_rawspan = None

def plot_mark(ax, row, value, kwargs):
    ax.plot([value, value], [row - 1/3, row + 1/3], **kwargs)

CENTER_LIMIT = -2
main_row_names = []
rawspan_row_names = []
scatter_kwargs = {'color': 'k', 'alpha': args.scatter, 'linewidths': 0}
median_kwargs = {'color': 'tab:orange', 'lw': 2}
extrema_kwargs = {'color': 'k', 'lw': 1}
min_psets = len(next(iter(param_data.values())))
max_psets = 0
for p, display_name in param_names.items():
    psets = param_data[p]
    min_psets = min(min_psets, len(psets))
    max_psets = max(max_psets, len(psets))
    rect_kwargs = {'alpha': 1 / len(psets), 'fill': True, 'fc': 'tab:blue'}
    if p in main_params:
        row = len(main_row_names)
        lefts = []
        rights = []
        for pset in psets:
            left = np.log10(1 - pset['min'] / pset['original'])
            lefts.append(left)
            if left > CENTER_LIMIT:
                left_rect = mplpatch.Rectangle((CENTER_LIMIT, row - 1/3), left - CENTER_LIMIT, 2/3, **rect_kwargs)
                ax_left.add_patch(left_rect)
            right = np.log10((pset['max'] - pset['original']) / pset['original'])
            rights.append(right)
            if right > CENTER_LIMIT:
                right_rect = mplpatch.Rectangle((CENTER_LIMIT, row - 1/3), right - CENTER_LIMIT, 2/3, **rect_kwargs)
                ax_right.add_patch(right_rect)
        if args.scatter:
            ax_left.scatter(lefts, [row] * len(psets), **scatter_kwargs)
            ax_right.scatter(rights, [row] * len(psets), **scatter_kwargs)
        if args.median:
            plot_mark(ax_left, row, np.max(lefts), extrema_kwargs)
            plot_mark(ax_right, row, np.max(rights), extrema_kwargs)
            plot_mark(ax_left, row, np.median(lefts), median_kwargs)
            plot_mark(ax_right, row, np.median(rights), median_kwargs)
        main_row_names.append(display_name)
    else:
        row = len(rawspan_row_names)
        spans = []
        for pset in psets:
            span = min(pset['max'], 1) - max(pset['min'], 0)
            spans.append(span)
            rect = mplpatch.Rectangle((0, row - 1/3), span, 2/3, **rect_kwargs)
            ax_rawspan.add_patch(rect)
        if args.scatter:
            ax_rawspan.scatter(spans, [row] * len(psets), **scatter_kwargs)
        if args.median:
            plot_mark(ax_rawspan, row, np.max(spans), extrema_kwargs)
            plot_mark(ax_rawspan, row, np.median(spans), median_kwargs)
        rawspan_row_names.append(display_name)

def minor_tick_positions(powers):
    absvals = []
    for power in powers:
        absvals.extend(np.array([2, 3, 4, 5, 6, 7, 8, 9]) * (10 ** power))
    return np.log10(absvals)

ax_left.set_xlim(0, CENTER_LIMIT)
ax_left.set_xticks([0, -1, -2])
ax_left.set_xticklabels(['100%', '10%', '1%'])
ax_left.set_xticks(minor_tick_positions([-2, -1]), minor=True)
ax_left.set_xlabel('Proportion Subtracted')
ax_left.set_ylim(len(main_params) - 1/2, -1/2)
ax_left.set_yticks(range(len(main_row_names)))
ax_left.set_yticklabels(main_row_names)
ax_right.set_xlim(CENTER_LIMIT, -CENTER_LIMIT)
ax_right.set_xticks([-1, 0, 1, 2])
ax_right.set_xticklabels(['10%', '100%', R'$\times 10$', R'$\times 100$'])
ax_right.set_xticks(minor_tick_positions([-2, -1, 0, 1]), minor=True)
ax_right.tick_params(axis='y', left=False, labelleft=False)
ax_right.set_xlabel('Proportion Added')
if rawspan_params:
    ax_rawspan.set_xlabel(R'$\mathrm{Maximum} - \mathrm{Minimum}$')
    ax_rawspan.set_xlim(0, 1)
    ax_rawspan.set_ylim(len(rawspan_params) - 1/2, -1/2)
    ax_rawspan.xaxis.set_minor_locator(mpltick.MultipleLocator(0.05))
    ax_rawspan.set_yticks(range(len(rawspan_row_names)))
    ax_rawspan.set_yticklabels(rawspan_row_names)

print(f'Data from {min_psets} <= n <= {max_psets} psets')

fig.tight_layout()
fig.savefig(args.output)
