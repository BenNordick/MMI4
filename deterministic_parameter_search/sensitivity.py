import argparse
import multiattractor
import numpy as np
import pandas as pd
import sobol
import tellurium as te

parser = argparse.ArgumentParser()
parser.add_argument('model', type=str, help='input SB file')
parser.add_argument('psets', type=str, help='CSV file of psets to test')
parser.add_argument('output', type=str, help='output CSV')
parser.add_argument('--attractors', '-a', type=int, help='test psets with N attractors')
parser.add_argument('--base', '-b', type=str, help='pset base Python expression file')
parser.add_argument('--init-species', '-i', nargs='+', help='species names with nonzero initial conditions')
parser.add_argument('--ignore-params', '-g', nargs='+', help='parameter names to ignore')
parser.add_argument('--ic-base', '-x', type=int, default=5, help='exponent base for initial conditions')
parser.add_argument('--ic-count', '-c', type=int, default=100, help='number of initial condition points')
parser.add_argument('--ic-min-power', '-p', type=float, default=-4.5, help='minimum initial condition power')
parser.add_argument('--ic-max-power', '-P', type=float, default=1.5, help='maximum initial condition power')
parser.add_argument('--verbose', '-v', action='store_true', help='print ranges')
parser.add_argument('--debug', action='store_true', help='double-check boundaries')
parser.add_argument('--patch', type=str, help='redo failed double-checks or missing values in existing CSV')
parser.add_argument('--patch-missing', action='store_true', help='use normal integrator settings for patching')
args = parser.parse_args()

psets_df = pd.read_csv(args.psets)
if args.attractors is None:
    attractors = int(psets_df['attractors'].max())
else:
    attractors = args.attractors
print(f'Testing psets with {attractors} attractors')
if 'attractors' in psets_df.columns:
    psets = psets_df[psets_df['attractors'] == attractors].to_dict('records')
else:
    psets = psets_df.to_dict('records')

with open(args.model) as f:
    runner = te.loadAntimonyModel(f.read())

if args.base:
    with open(args.base) as f:
        pset_base = eval(f.read())
else:
    pset_base = {}

if args.init_species:
    ic_power_range = args.ic_max_power - args.ic_min_power
    nonzero_ini_combs = np.power(args.ic_base, sobol.sample(len(args.init_species), args.ic_count, skip=args.ic_count) * ic_power_range + args.ic_min_power)
    ini_combs = np.zeros((nonzero_ini_combs.shape[0], len(runner.fs())))
    next_rna_col = 0
    for i, species in enumerate(runner.fs()):
        if species in args.init_species:
            ini_combs[:, i] = nonzero_ini_combs[:, next_rna_col]
            next_rna_col += 1
    update_protein_inicombs = lambda: None
else:
    ini_combs, update_protein_inicombs = multiattractor.makeinitialconditions(runner, args.ic_count, 
        base=args.ic_base, minexp=args.ic_min_power, maxexp=args.ic_max_power)

known_parameters = set(runner.ps())
def unmodifiable_parameter(p):
    return any(p.startswith(bad) for bad in ['kOff', 'kOn'])
parameters_to_test = sorted(p for p in known_parameters if not unmodifiable_parameter(p))
if args.ignore_params:
    parameters_to_test = [p for p in parameters_to_test if p not in args.ignore_params]
print(f'Testing parameters {sorted(parameters_to_test)}')

def reset_parameters(pset):
    for p, value in {**pset_base, **pset}.items():
        if p in known_parameters:
            runner[p] = value

if args.patch and not args.patch_missing:
    fpa_kwargs = {'init_step': 0.25, 'min_step': 0.01, 'init_time': 500, 'max_time': 4000}
else:
    fpa_kwargs = {}

for pset in psets:
    reset_parameters(pset)
    attractor_points = multiattractor.findpointattractors(runner, ini_combs, **fpa_kwargs)
    if len(attractor_points) != attractors:
        raise RuntimeError(f'expected {attractors} attractors but found {len(attractor_points)}; pset: {pset}')
print(f'Verified attractor counts ({len(psets)} psets)')

def check_attractors():
    attractor_points = multiattractor.findpointattractors(runner, ini_combs, **fpa_kwargs)
    return len(attractor_points) == attractors

def binary_log_search(param, min_value, max_value, going_up):
    left = np.log2(min_value)
    right = np.log2(max_value)
    while right - left > 0.0001:
        half_range = (right - left) / 2
        runner[param] = 2 ** (left + half_range)
        if check_attractors() == going_up:
            left += half_range
        else:
            right -= half_range
    return 2 ** ((left + right) / 2)

def limit_search(param, initial_value, multiplier):
    prev = initial_value
    cur = runner[param] = initial_value * multiplier
    steps = 0
    while check_attractors():
        if steps > 10:
            return cur
        prev = cur
        cur = runner[param] = cur * multiplier
        steps += 1
    left_end = min(prev, cur)
    right_end = max(prev, cur)
    return binary_log_search(param, left_end, right_end, going_up=(multiplier > 1))

def double_check(param, original_value, left_end, right_end):
    runner[param] = min(left_end * 1.01, original_value)
    if not check_attractors():
        print(f'*** not stable above {param} left_end {left_end}')
        return False
    if left_end > original_value / (2 ** 10):
        runner[param] = left_end * 0.99
        if check_attractors():
            print(f'*** stable below {param} left_end {left_end}')
            return False
    runner[param] = max(right_end * 0.99, original_value)
    if not check_attractors():
        print(f'*** not stable below {param} right_end {right_end}')
        return False
    if right_end < original_value * (2 ** 10):
        runner[param] = right_end * 1.01
        if check_attractors():
            print(f'*** stable above {param} right_end {right_end}')
            return False
    return True

if args.patch:
    results = pd.read_csv(args.patch).to_dict('records')
else:
    results = []

for i, pset in enumerate(psets):
    print(f'{attractors}-stable pset #{i}')
    reset_parameters(pset)
    if args.patch:
        result = results[i]
        pset_params = [p for p in parameters_to_test if not result.get(f'validated_{p}', False)]
    else:
        result = {}
        pset_params = parameters_to_test
    for param in pset_params:
        original_value = runner[param]
        result[f'original_{param}'] = original_value
        left_end = limit_search(param, original_value, 1/2)
        result[f'min_{param}'] = left_end
        right_end = limit_search(param, original_value, 2)
        result[f'max_{param}'] = right_end
        if args.verbose:
            print(f'  {param}\t{left_end} < {original_value} < {right_end}')      
        if args.debug:
            result[f'validated_{param}'] = double_check(param, original_value, left_end, right_end)
        runner[param] = original_value
    if not args.patch:
        results.append(result)
    pd.DataFrame.from_dict(results).to_csv(args.output, index=False)
    if args.verbose:
        print()
